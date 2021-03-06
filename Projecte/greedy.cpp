/***************************************************************************
    greedy.cpp 
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
vector<int> dominant_set;

// string for keeping the name of the input file
string inputFile;

// dummy parameters as examples for creating command line parameters 
// see function read_parameters(...)
int dummy_integer_parameter = 0;
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
        if (strcmp(argv[iarg],"-i")==0) inputFile = argv[++iarg];
        
        // example for creating a command line parameter param1 
        //-> integer value is stored in dummy_integer_parameter
        else if (strcmp(argv[iarg],"-param1")==0) 
            dummy_integer_parameter = atoi(argv[++iarg]); 
            
        // example for creating a command line parameter param2 
        //-> double value is stored in dummy_double_parameter
        else if (strcmp(argv[iarg],"-param2")==0) 
            dummy_double_parameter = atof(argv[++iarg]);  
            
        iarg++;
    }
}

bool is_PIDS() {
    for (int n = 0; n < neighbors.size(); ++n) {
        float counter = 0;
        for (int m : neighbors[n]) {
            if (dominant_set[m] == 1) ++counter;
        }
        if (neighbors[n].size() != 0 and counter/neighbors[n].size() >= 0.5) return false;
    }
    return true;
}


/************
Main function
*************/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);
    
    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // initializing the random number generator. 
    // A random number in (0,1) is obtained with: double rnum = rnd->next();
    rnd = new Random((unsigned) time(&t));
    rnd->next();

    // variables for storing the result and the computation time 
    // obtained by the greedy heuristic
    double results = std::numeric_limits<int>::max();
    double time = 0.0;

    // opening the corresponding input file and reading the problem data
    ifstream indata;
    indata.open(inputFile.c_str());
    if(not indata) { // file couldn't be opened
        cout << "Error: file could not be opened" << endl;
    }

    indata >> n_of_nodes;
    indata >> n_of_arcs;
    neighbors = vector< set<int> >(n_of_nodes);
    dominant_set = vector<int> (n_of_nodes, 0);
    int u, v;
    while(indata >> u >> v) {
        neighbors[u - 1].insert(v - 1);
        neighbors[v - 1].insert(u - 1);
    }
    indata.close();

    // the computation time starts now
    Timer timer;

    // Example for requesting the elapsed computation time at any moment: 
    // double ct = timer.elapsed_time(Timer::VIRTUAL);

    // HERE GOES YOUR GREEDY HEURISTIC
    // When finished with generating a solution, first take the computation 
    // time as explained above. Say you store it in variable ct.
    // Then write the following to the screen: 
    // cout << "value " << <value of your solution> << "\ttime " << ct << endl;
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

    double ct = timer.elapsed_time(Timer::VIRTUAL);
    int size = 0;
    cout << "Dominant set:" << endl;
    for (int i = 0; i < n_of_nodes; ++i) {
        if (dominant_set[i] == 1) {
            ++size;
            cout << i << " ";
        }
    }
    cout << endl;
    cout << "Is it PIDS? " << (is_PIDS() ? "Yes" : "No") << endl;
    cout << "PIDS size = " << size << endl;
    cout << "Finished in " << ct << " seconds" << endl;
}

