#pragma once

#include "helpers.h"

class SolverLanczos
{   
    constexpr int n(int j, int l) const {
        return (j>>l) & 1;
    }

    constexpr int fermion_sign(int j, int l) const {
        return (__builtin_popcount(j & ((1 << l) - 1)) & 1) ? -1 : 1;
    }

    long double Z;
    int n_bath_max = -1, NFock = 0;
    int n_bath;

    VectorXcd Iw;
    VectorXd  Iw2;

    VectorXi nud;
    SparseMatrixXd H;
    MatrixXd Psi;
    VectorXd E;

    VectorXd t2u, t2d, eu, ed; 
    double adjust_bath_levels(const VectorXcd& Delta_w, VectorXd& t2_bath, VectorXd& e_bath);
    double adjust_bath_couplings(const VectorXcd& Delta_w, VectorXd& t2_bath, VectorXd& e_bath);

    int get_block_idx(const VectorXd& psi) const;
    
public:
    double tolerance;
    double Beta;
    int Nbath_max;
    int Nw, NwED;
    
    SolverLanczos(double Beta, int Nbath_max, int Nw, int NwED=20);

    VectorXcd gu, gd;
    VectorXcd sigmau, sigmad;

    double init(double U, const VectorXcd& Delta_up, const VectorXcd& Delta_down, 
        int Nbath, int Nm, int Nkr, double h_loc=0, double mu_loc=0);
};
