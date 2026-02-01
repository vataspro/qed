#include "mc.h"
#include <iostream>

/*
    Lattice Constructor
*/
MCRunner::MCRunner(U1& field_, double eps_) : acc(0), totiter(0), eps(eps_), field(field_) {}


// Metropolis Step
void MCRunner::step(int x, int mu) {

    double p;
    ranlxs(&p, 1);
    
    double dx;
    ranlxs(&dx, 1);
    dx = eps * (dx - 0.5);

    double dS;
    dS = field.S_loc(x, mu);

    
    double temp1;
    temp1 = field.U[x][mu];
    //cout << temp1 << endl;
    
    double old = field.U[x][mu];
    field.U[x][mu] += dx;

    field.U[x][mu] = fmod(field.U[x][mu] + 2*M_PI, (2 * M_PI));
    dS = field.S_loc(x, mu) - dS; // Snew - Sold

    // pacc = exp(-dS)
    if (exp(-dS) < p) { //if reject
        field.U[x][mu] = old; //-= dx;

    } else { acc++;
//        cout << eps << " " << dx << " " << temp1 << " " << field.U[x][mu] << endl;
    }

    totiter++;

}

void MCRunner::sweep(int local_rep) {
    for (int x=0; x<pow(field.N, 4); x++) 
    for (int mu=0; mu<4; mu++)
    for (int loc_i=0; loc_i<local_rep; loc_i++) {
            step(x, mu);
        }
}

void MCRunner::overlx_step(int x, int mu) {

    double theta_staple = arg(field.staple(x, mu));
    field.U[x][mu] = 2. * theta_staple - field.U[x][mu];

    // ensure [0, 2*pi)
    // Max angle is 8*pi
    field.U[x][mu] = fmod(field.U[x][mu] + 8*M_PI, 2*M_PI);
}

void MCRunner::overlx_sweep() {

    int x, mu;
    for (int i=0; i<field.V; i++) {
        x = randint(field.V);
        mu = randint(Nd);
        overlx_step(x, mu);
    }
}


