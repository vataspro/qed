#include "lattice.h"

#include <iostream>

int randint(int N) {
    double x;
    ranlxs(&x, 1);
    return static_cast<int>(N*x);
}

U1::U1(int N_, double beta_) : U(pow(N_, Nd)), fwd(pow(N_, Nd)), bwd(pow(N_, Nd)), N(N_), beta(beta_), V(pow(N_, Nd)) {

    // Bad hot start, as ranlxs is U[0, 1] rather than [0, 2pi]
    for (auto& v : U) {
        ranlxs(v.data(), v.size()); 
        for (auto& theta : v) {
            theta *= 2.0 * M_PI;
        }
    }

    /*
    for (int x = 0; x < V; x++)
    for (int mu = 0; mu < Nd; mu++)
    if (U[x][mu] != 0.0) {
        cout << "BAD LINK " << x << " " << mu
                  << " = " << U[x][mu] << endl;
        abort();
    }*/

    // Get forward neighbours
    for (int x=0; x<N; x++)
    for (int y=0; y<N; y++) 
    for (int z=0; z<N; z++) 
    for (int t=0; t<N; t++) {

        fwd[index(x, y, z, t)][0] = index((x+1)%N, y, z, t);
        fwd[index(x, y, z, t)][1] = index(x, (y+1)%N, z, t);
        fwd[index(x, y, z, t)][2] = index(x, y, (z+1)%N, t);
        fwd[index(x, y, z, t)][3] = index(x, y, z, (t+1)%N);

        bwd[index(x, y, z, t)][0] = index((x-1+N)%N, y, z, t);
        bwd[index(x, y, z, t)][1] = index(x, (y-1+N)%N, z, t);
        bwd[index(x, y, z, t)][2] = index(x, y, (z-1+N)%N, t);
        bwd[index(x, y, z, t)][3] = index(x, y, z, (t-1+N)%N);

    }
}


complex<double> U1::staple(int x, int mu) {

    complex<double> sum = 0;
    
    for (int nu = 0; nu < Nd; nu++) {
        if (nu == mu) continue;

        int x_minus_nu = bwd[x][nu];

        // forward plaquette
        sum += polar(1.0,
            U[x][nu]
          + U[fwd[x][nu]][mu]
          - U[fwd[x][mu]][nu]
        );

        // backward plaquette
        sum += polar(1.0,
          - U[x_minus_nu][nu]
          + U[x_minus_nu][mu]
          + U[fwd[x_minus_nu][mu]][nu]
        );
    }

    return sum;
}
            
double U1::S_loc(int x, int mu) {

    double S = 0;
    for (int nu = 0; nu < Nd; nu++) {
        if (nu == mu) continue;

        // forward plaquette
        double Umn = U[x][mu] + U[fwd[x][mu]][nu] - U[fwd[x][nu]][mu] - U[x][nu];
        S += cos(Umn);

        // backward plaquette
        int xmnu = bwd[x][nu];
        Umn = U[xmnu][mu] + U[fwd[xmnu][mu]][nu] - U[fwd[xmnu][nu]][mu] - U[xmnu][nu];
        S += cos(Umn);
    }
    return - beta * S;
}


double U1::plaquette(int x, int mu, int nu) {
    return U[x][mu] + U[fwd[x][mu]][nu] - U[fwd[x][nu]][mu] -     U[x][nu];
}
double U1::plaquette() {

    double sum = 0;
    for (int x=0; x<V; x++) {
        // six plaquettes per point in 4d
        // one for each pair 01, 02, 03, 12, 13, 23
        for (int mu=0; mu<3; mu++)
        for (int nu=mu+1; nu<4; nu++) {

            double Umn = plaquette(x, mu, nu);//U[x][mu] + U[fwd[x][mu]][nu] - U[fwd[x][nu]][mu] - U[x][nu];
            sum += cos(Umn);

        }
    }

    return sum / (6 * V);
}

double U1::monopole_density(int x, int mu, int nu) {
    double plaq = plaquette(x, mu, nu);
    //M = fmod(M + 4 * M_PI, 2 * M_PI);
    int M=0;
    while (plaq > M_PI) {
        ++M;
        plaq -= 2* M_PI;
    }
    while (plaq < -M_PI) {
        --M;
        plaq += 2*M_PI;
    }

    return static_cast<double>(M);
}


double U1::magnetic_flux(int x, int mu, int nu, int rho) {
    int x_plus_rho = fwd[x][rho];
    return monopole_density(x_plus_rho, mu, nu) - monopole_density(x, mu, nu);
}

double U1::monopole_density() {

    double M = 0;
    for (int x=0; x<V; x++) {
        // 0 component
        M += abs( + magnetic_flux(x, 1, 2, 3) 
                  - magnetic_flux(x, 1, 3, 2)
                  + magnetic_flux(x, 2, 3, 1));
        /*
        // 1 component
        M += abs( + magnetic_flux(x, 0, 2, 3) 
                  - magnetic_flux(x, 0, 3, 2)
                  + magnetic_flux(x, 2, 3, 0));
        // 2 component
        M += abs( - magnetic_flux(x, 0, 1, 3) 
                  + magnetic_flux(x, 0, 3, 1)
                  - magnetic_flux(x, 1, 3, 0));
        // 3 component
        M += abs( + magnetic_flux(x, 0, 1, 2) 
                  - magnetic_flux(x, 0, 2, 1)
                  + magnetic_flux(x, 1, 2, 0));
        */
    }

        return M / (4. * V) * N;
};


int U1::index(int x, int y, int z, int t) {
    return N*(N*(N*x + y) + z) + t;
}
