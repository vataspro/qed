/*
    Code for a 4D lattice containing a U1 field.


    Because the lattice is 4D, the field also has Nd=4
    components.

*/
#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <array>
#include "ranlxs.h"
#include <complex>

using namespace std;

constexpr int Nd = 4;

int randint(int N);

class U1 {

  public:
    vector<array<double, Nd>> U; // Gauge field angle theta
    int N; // Lattice dimension length
    int V; // Lattice volume
    double beta;

    U1(int N_, double beta_); // Constructor of N^4 lattice field

    int index(int x, int y, int z, int t); // index of point

    double plaquette(int x, int mu, int nu); // local plaquette
    double plaquette(); // average plaquette

    complex<double> staple(int x, int mu);
    double S_loc(int x, int mu);
    vector<array<int, Nd>> fwd;
    vector<array<int, Nd>> bwd;

    double monopole_density(int x, int mu, int nu);
    double monopole_density();
    double magnetic_flux(int x, int mu, int nu, int rho);


    //int fwd(int nu); // Forward neighbour in nu-direction 
    //int bwd(int nu); // Backward --

};

