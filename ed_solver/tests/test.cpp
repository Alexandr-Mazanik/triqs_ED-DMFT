#include "solver.h"

int main() {
	double U = 6, t = 1, Beta = 5, h = 0.01, mu = 0;
	int Nbath = 3;
	int Nw = 200;

	Solver solver(Beta, U, t, mu, h, Nbath, Nw);

	complex* Delta_up   = new complex[Nw];
    complex* Delta_down = new complex[Nw];

    for(int w = 0; w < Nw; w++) {
        Delta_up[w]   = complex(0.0, -0.5);  
        Delta_down[w] = complex(0.0, -0.5);
    }

    std::cout << "Starting ED solve...\n";
    double lnZ = solver.init(Delta_up, Delta_down, Nbath, h, mu);

    std::cout << "\nlnZ_imp = " << lnZ << "\n";
    std::cout << "G_up(iw_0) = " << solver.gu[0] << "\n";
    std::cout << "Sigma_up(iw_0) = " << solver.sigmau[0] << "\n";

    delete[] Delta_up;
    delete[] Delta_down;

    return 0;
}
