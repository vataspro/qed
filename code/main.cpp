#include <iostream>
#include <iomanip>
#include "lattice.h"
#include "mc.h"

using namespace std;

int main(int argc, char* argv[]) {


    if (argc != 2) {
        cerr << "Usage " << argv[0] << " beta" << endl;
        exit(1);
    }

    const int N = 4;
    double beta = atof(argv[1]);
    double eps=0.1;

    U1 lattice(N, beta);
    MCRunner runner(lattice, eps);

    for (int j=0; j<10000; j++) {
        runner.sweep();
        runner.overlx_sweep();
    }

    for (int k=0; k<200; k++) {
        for (int i=0; i<200; i++) {
        runner.sweep();
        runner.overlx_sweep();
        //cout << static_cast<double>(runner.acc) / (lattice.V * 4 * static_cast<double>(100 + 100*k + i)) << endl;
        }
        cout << setprecision(16) << lattice.plaquette() << " ";
        cout << setprecision(16) << lattice.monopole_density() << endl;
    }

    return 0;
}
