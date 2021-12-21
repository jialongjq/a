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
string resultFile;

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

        else if (strcmp(argv[iarg],"-r")==0) resultFile = argv[++iarg];

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
        if (neighbors[n].size() != 0 and counter/neighbors[n].size() < 0.5) return false;
    }
    return true;
}

bool is_minimal() {
  	for (int n = 0; n < n_of_nodes; ++n) {
		if (dominant_set[n] == 1) {
            bool pids = true;
			dominant_set[n] = 0;
			for (int m : neighbors[n]) {
				int dominant_neighbors_of_m = 0;
				for (int l : neighbors[m]) if (dominant_set[l] == 1) ++dominant_neighbors_of_m;
				if (neighbors[m].size() != 0 and float(dominant_neighbors_of_m)/neighbors[m].size() < 0.5) pids = false;
			}
			if (pids) return false;
			dominant_set[n] = 1;
		}
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

    Timer timer;
    double ct1 = timer.elapsed_time(Timer::VIRTUAL);


    // opening the corresponding input file and reading the problem data
    ifstream indata;
    indata.open(inputFile.c_str());
    if(not indata) { // file couldn't be opened
        cout << "Error: file could not be opened" << endl;
    }

    indata >> n_of_nodes;
    indata >> n_of_arcs;
    neighbors = vector<set<int>>(n_of_nodes);
    int u, v;
    while(indata >> u >> v) {
        neighbors[u - 1].insert(v - 1);
        neighbors[v - 1].insert(u - 1);
    }
    indata.close();

    ifstream inresult;
    inresult.open(resultFile.c_str());
    if(not inresult) { // file couldn't be opened
        cout << "Error: file could not be opened" << endl;
    }
    int x;
    dominant_set = vector<int>(n_of_nodes, 0);
    while (inresult >> x) {
      dominant_set[x] = 1;
    }
    //for (int i = 0; i < dominant_set.size(); ++i) cout << " " << dominant_set[i];

    if (is_PIDS()) {
      cout << "La solució donada es dominador d'influencia positiva";
      if (is_minimal()) cout << " i minimal." << endl;
      else cout << " pero no minimal." << endl;
    }
    else cout << "La solució donada no és mínim conjunt dominador d'influència positiva."<< endl;
    // Example for requesting the elapsed computation time at any moment:
    double ct2 = timer.elapsed_time(Timer::VIRTUAL);
    cout << "\ttime " << ct2-ct1 << endl;

}
