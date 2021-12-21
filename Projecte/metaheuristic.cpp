/***************************************************************************
    metaheuristic.cpp 
    (C) 2021 by C. Blum & M. Blesa
    
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Timer.h"
#include "Random.h"
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <vector>
#include <set>
#include <limits>
#include <iomanip>
#include <queue>
#include <math.h>

// global variables concerning the random number generator (in case needed)
time_t t;
Random* rnd;

// Data structures for the problem data
int n_of_nodes;
int n_of_arcs;
vector< set<int> > neighbors;

// string for keeping the name of the input file
string inputFile;

// computing time limit for each application of the metaheuristic
double time_limit = 3200.0;

// number of applications of the metaheuristic
int n_apps = 1;

// dummy parameters as examples for creating command line parameters 
// (see function read_parameters(...))
int n_population = 0;
int dummy_double_parameter = 0.0;


inline int stoi(string &s) {

  return atoi(s.c_str());
}

inline double stof(string &s) {

  return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

    int iarg = 1;
    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i") == 0) inputFile = argv[++iarg];
        // reading the computation time limit 
        // from the command line (if provided)
        else if (strcmp(argv[iarg],"-t") == 0) time_limit = atoi(argv[++iarg]); 
        // reading the number of applications of the metaheuristic 
        // from the command line (if provided)
        else if (strcmp(argv[iarg],"-n_apps") == 0) n_apps = atoi(argv[++iarg]); 
        // example for creating a command line parameter 
        // param1 -> integer value is stored in dummy_integer_parameter
        else if (strcmp(argv[iarg],"-n_population") == 0) {
            n_population = atoi(argv[++iarg]); 
        }
        // example for creating a command line parameter 
        // param2 -> double value is stored in dummy_double_parameter
        else if (strcmp(argv[iarg],"-param2") == 0) {
            dummy_double_parameter = atof(argv[++iarg]);
        }
        iarg++;
    }
}

/***************************
Genetic Algorithm parameters
***************************/
vector<vector<int>> population;

vector<int> dominant_set;
    
vector<int> valid_vertexs; // Saves the valid vertexs (if valid_vertexs[i] == total number of vertexs, then the individual i is PIDS)
vector<int> ds_sizes; // Saves the size of the dominant set of each individual
vector<int> fitness; // Saves the fitness for each individual of the current population

vector<int> candidate; // Saves the candidate individual
int candidate_i; // Indicates the index of the candidate in the current population
int candidate_size;
int last_candidate_size;
bool finish;

/************************
Genetic Algorithm methods
************************/

bool is_PIDS() {
    for (int n = 0; n < neighbors.size(); ++n) {
        float counter = 0;
        for (int m : neighbors[n]) {
            if (dominant_set[m] == 1) ++counter;
        }
        if (neighbors[n].size() != 0 and counter/neighbors[n].size() < 0.5) return false;
    }
    return true;
}


void initialize_population() {
    priority_queue<pair<int, int>> pq;
    for (int n = 0; n < n_of_nodes; ++n) pq.push(make_pair(-neighbors[n].size(), n));
    while (not pq.empty()) {
        pair<int, int> p = pq.top();

        priority_queue<pair<int, int>> adyacencies;
        int marked = 0;
        for (int m : neighbors[p.second]) {
            if (dominant_set[m] == 0) adyacencies.push(make_pair(neighbors[m].size(), m));
            else ++marked;
        }

        int counter = ceil(float(neighbors[p.second].size())/2) - marked;
        while (counter > 0) {
            dominant_set[adyacencies.top().second] = 1;
            --counter;
            adyacencies.pop();
        }
        pq.pop();
    }
    int size = 0;
    for (int i = 0; i < n_of_nodes; ++i) {
        if (dominant_set[i] == 1) ++size;
    }
    // Initialize the population
    srand(time(NULL));
    for (int i = 0; i < population.size(); ++i) {
        /*for (int j = 0; j < individual.size(); ++j) {
            individual[j] = 1; //rand()%2;
        }*/
        population[i] = dominant_set;
    }
    return;
}

void calculate_ds_sizes() {
    for (int i = 0; i < population.size(); ++i) {
        ds_sizes[i] = 0;
        for (int j = 0; j < population[i].size(); ++j) {
            if (population[i][j] == 1) ++ds_sizes[i];
        }
    }
    return;
}

void calculate_valid_vertexs() {
    for (int i = 0; i < population.size(); ++i) {
        // Iterate through each vertex n
        valid_vertexs[i] = 0;
        for (int n = 0; n < n_of_nodes; ++n) {
            int counter = 0;
            // Iterate through each adjacent vertex G[n][m]
            for (int m : neighbors[n]) {
                if (population[i][m] == 1) ++counter;
            }
            if (neighbors[n].size() == 0 or float(counter)/neighbors[n].size() >= 0.5) {
                ++valid_vertexs[i];
            }
        }
    }
    return;
}

void find_candidate() {
    int min;
    if (candidate_i == -1) min = n_of_nodes;
    else min = candidate_size;
    for (int i = 0; i < population.size(); ++i) {
        if (valid_vertexs[i] == n_of_nodes) { // is PIDS
            if (ds_sizes[i] < min) {
                min = ds_sizes[i];
                candidate = population[i];
                candidate_i = i;
                candidate_size = ds_sizes[i];
            }
        }
    }
    return;
}

void calculate_fitness() {
    for (int i = 0; i < population.size(); ++i) {
        if (candidate_i == -1) {
            fitness[i] = int((float(valid_vertexs[i])/n_of_nodes)*100);
        }
        else {
            if (valid_vertexs[i] == n_of_nodes) { // if it is PIDS
                fitness[i] = 100 - 10*(ds_sizes[i] - candidate_size); //int(100 - (float(100)/population[i].size() * (ds_sizes[i] - candidate_size)));
            }
            else { // if it is not PIDS
                fitness[i] = 0; //int(float(valid_vertexs[i])/n_of_nodes*100);
            }
        }
        //cout << " " << fitness[i];
    }
    bool found = false;
    int i = 0;
    while (not found and i < fitness.size()) {
        if (fitness[i] > 0) found = true;
        ++i;
    }
    if (not found) finish = true;
    return;
}

void calculate_fitness2() {
    for (int i = 0; i < population.size(); ++i) {
        fitness[i] = float(valid_vertexs[i])/n_of_nodes * 100;
        cout << " " << fitness[i];
    }
    cout << endl;
    return;
}

void selection() {
    vector<int> selected = vector<int>(n_population);
    int k = 0;
    while (k < n_population) {
        for (int i = 0; i < n_population; ++i) {
            int random = rand()%101;
            if (random <= fitness[i]) {
                selected[k] = i;
                ++k;
                if (k == n_population) break;
            }
        }
    }
    vector<vector<int>> new_population = vector<vector<int>>(n_population);
    for (int i = 0; i < n_population; ++i) {
        new_population[i] = population[selected[i]];
    }
    population = new_population;
    return;
}

void crossover(int probability) {
    for (int i = 0; i < n_population; i += 2) {
        if (i + 1 < n_population) {
            if (rand()%101 <= probability) {
                int random = rand()%n_of_nodes + 1;
                for (random; random < n_population; ++random) {
                    int aux = population[i][random];
                    population[i][random] = population[i+1][random];
                    population[i+1][random] = aux;
                }
            }
        }
    }
}

int mutate(int a) {
    if (a == 0) return 1;
    return 0;
}

void mutation(int probability) {
    for (int i = 0; i < n_population; ++i) {
        if (rand()%101 <= probability) {
            int j = rand()%n_of_nodes;
            population[i][j] = mutate(population[i][j]);
        }
    }
}

void print_individuals() {
    for (int i = 0; i < population.size(); ++i) {
        for (int j = 0; j < population[i].size(); ++j) {
            cout << population[i][j];
        }
        cout << endl;
    }
    return;
}

/**********
Main function
**********/

int main( int argc, char **argv ) {
    srand(time(NULL));

    read_parameters(argc,argv);
    
    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // initializing the random number generator. A random number 
    // between 0 and 1 is obtained with: double rnum = rnd->next();
    rnd = new Random((unsigned) time(&t));
    rnd->next();

    // vectors for storing the result and the computation time 
    // obtained by the <n_apps> applications of the metaheuristic
    vector<double> results(n_apps, std::numeric_limits<int>::max());
    vector<double> times(n_apps, 0.0);

    // opening the corresponding input file and reading the problem data
    ifstream indata;
    indata.open(inputFile.c_str());
    if (not indata) { // file couldn't be opened
        cout << "Error: file could not be opened" << endl;
    }

    indata >> n_of_nodes;
    indata >> n_of_arcs;
    neighbors = vector< set<int> >(n_of_nodes);
    int u, v;
    while (indata >> u >> v) {
        neighbors[u - 1].insert(v - 1);
        neighbors[v - 1].insert(u - 1);
    }
    indata.close();

    // main loop over all applications of the metaheuristic
    for (int na = 0; na < n_apps; ++na) {

        // the computation time starts now
         Timer timer;

        // Example for requesting the elapsed computation time at any moment: 
        //double ct = timer.elapsed_time(Timer::VIRTUAL);

        cout << "start application " << na + 1 << endl;

        // HERE GOES YOUR METAHEURISTIC

        // For implementing the metaheuristic you probably want to take profit 
        // from the greedy heuristic and/or the local search method that you 
        // already developed.
        //
        // Whenever the best found solution is improved, first take the 
        // computation time as explained above. Say you store it in variable ct.
        // Then, write the following to the screen: 
        // cout << "value " << <value of the new best found solution>;
        // cout << "\ttime " << ct << endl;
        //
        // Store the value of the new best found solution in vector results: 
        // results[na] = <value of the new best found solution>;
        //
        // And store the current computation time (that is, the time measured 
        // at that moment and stored in variable "ct") in vector times: 
        // times[na] = ct;
        //
        // Stop the execution of the metaheuristic 
        // once the time limit "time_limit" is reached.
        
        // Selection of couples
        // initialization of the genetic algorithm parameters
        population = vector<vector<int>>(n_population);
        valid_vertexs = vector<int>(n_population, 0);
        ds_sizes = vector<int>(n_population, 0);
        
        dominant_set = vector<int> (n_of_nodes, 0);

        // Applies greedy once to get an individual, each individual of the population will be this one
        initialize_population();

        // Select any individual from the initial population
        candidate = population[0];
        candidate_i = 0;
        candidate_size = 0;
        for (int i = 0; i < candidate.size(); ++i) if (candidate[i] == 1) ++candidate_size;
        last_candidate_size = candidate_size;

        // All individuals starts being the same one, and it is PIDS so fitness will be scored as 100
        fitness = vector<int>(n_population, 100);

        finish = false;

        while (not finish) {
            double ct = timer.elapsed_time(Timer::VIRTUAL);
            selection();
            crossover(100);
            mutation(100);
            calculate_ds_sizes();
            calculate_valid_vertexs();
            find_candidate();
            if (candidate_i >= 0 && candidate_size < last_candidate_size) {
                last_candidate_size = candidate_size;
                cout << "PIDS found, size = " << candidate_size;
                cout << ", time: " << ct << " seconds" << endl;
                //for (int i = 0; i < candidate.size(); ++i) cout << candidate[i];
                //cout << endl;
                results[na] = candidate_size;
            }
            calculate_fitness();
            times[na] = ct;
            if (ct > time_limit) finish = true;
        }
        cout << "Dominant set:" << endl;
        for (int i = 0; i < candidate.size(); ++i) {
            if (candidate[i] == 1) cout << i << " ";
        }
        cout << endl;
        cout << "PIDS size = " << candidate_size << endl;
        cout << "Is it PIDS? " << (is_PIDS() ? "Yes" : "No") << endl;

        cout << "end application " << na + 1 << endl;
    }
    // calculating the average of the results and computation times, 
    // and their standard deviations, and write them to the screen
    double r_mean = 0.0;
    int r_best = std::numeric_limits<int>::max();
    double t_mean = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        r_mean = r_mean + results[i];
        if (int(results[i]) < r_best) r_best = int(results[i]);
        t_mean = t_mean + times[i];
    }
    r_mean = r_mean/double(results.size());
    t_mean = t_mean/double(times.size());
    double rsd = 0.0;
    double tsd = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        rsd = rsd + pow(results[i]-r_mean,2.0);
        tsd = tsd + pow(times[i]-t_mean,2.0);
    }
    rsd = rsd/double(results.size());
    if (rsd > 0.0) {
        rsd = sqrt(rsd);
    }
    tsd = tsd/double(results.size());
    if (tsd > 0.0) {
        tsd = sqrt(tsd);
    }
    // printing statistical information
    cout << r_best << "\t" << r_mean << "\t" << rsd << "\t";
    cout << t_mean << "\t" << tsd << endl;
}

