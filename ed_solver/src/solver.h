#pragma once

#include "helpers.h"

class Solver
{   
    constexpr int n(int j, int l) const {
        return (j>>l) & 1;
    }

    int n_bath_max, NFock, n_bath;

    double ** H, ** Psi, ** cuPsi,** cdPsi,  *E;

    long double Z;

    double * t2u, * t2d, *eu, *ed; 
    double adjust_bath_levels(complex * Delta_w, double * t2_bath, double * e_bath);
    double adjust_bath_couplings(complex * Delta_w, double * t2_bath,  double * e_bath);

    int * nud;
    
public:
    double Beta;
    int Nbath_max;
    int Nw, NwED;

    double tolerance;
    Solver(double Beta, int Nbath_max, int Nw, int NwED=20);
    complex * Iw;
    double * Iw2;

    complex * gu, * gd;
    complex * sigmau, * sigmad;

    double init(double U, complex * Delta_up, complex *Delta_down, int Nbath, double h_loc=0, double mu_loc=0);
    
    ~Solver();
};
