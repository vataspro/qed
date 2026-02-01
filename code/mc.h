/*
    4D scalar field lattice
*/
#pragma once
#include <iostream>
#include <string>
#include <iterator>
#include <fstream>
#include <vector>
#include <cmath>
#include <array>
#include "ranlxs.h"
#include "lattice.h"

using namespace std;

class MCRunner {
  public:
    U1& field;

    long int acc;
    long int totiter;
    double eps;

    double S_glob();
    double S_loc(int x, int mu);
    void step(int x, int mu);
    void sweep(int local_rep=1);

    void overlx_step(int x, int mu); 
    void overlx_sweep();

    // Constructor
    MCRunner(U1& field_, double eps_);
};


