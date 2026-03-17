import numpy as np
from dataclasses import dataclass, field

from triqs.gf import *

from ed_solver import EDSolverCore
from ed_solver import EDLSolverCore


@dataclass
class Bath:
    eps_up: list[float] = field(default_factory=list)
    eps_down: list[float] = field(default_factory=list)
    t2_up: list[float] = field(default_factory=list)
    t2_down: list[float] = field(default_factory=list)


class Solver:
    def __init__(self, beta, nbath_max=3, n_iw=1025, n_wED=20, method='Jacobi'):
        if method == 'Jacobi':
            self.cpp_solver = EDSolverCore(beta, nbath_max, n_iw, n_wED)
        elif method == 'Lanczos':
            self.cpp_solver = EDLSolverCore(beta, nbath_max, n_iw, n_wED)
        else:
            raise RuntimeError(f"Invalid method: {method} - expected 'Jacobi' or 'Lanczos'.")

        self.beta = beta
        self.n_iw = n_iw       
        self.lnZ = None

        self.nbath_max = nbath_max
        self.bath_parameters = Bath()
        
        self.method = method
        self.warn_seen = False

        self.iw_mesh = MeshImFreq(beta=beta, S='Fermion', n_max=n_iw)
        self.Sigma_iw = BlockGf(mesh=self.iw_mesh, gf_struct=[('up', 1), ('down', 1)])
        self.Sigma_iw.zero()
        self.G_iw = self.Sigma_iw.copy()
        self.Delta_iw = self.Sigma_iw.copy()
        
        # default initial guess for Delta
        self.set_initial_guess()
    
    def extend_to_negative_freq(self, gf_positive):
        gf_positive.data[0 : self.n_iw] = np.conj(gf_positive.data[2 * self.n_iw : self.n_iw - 1 : -1])

    def set_initial_guess(self, bath=None):
        if bath == None:
            iw_arr = np.array(list(self.iw_mesh.values()))
            iw2_arr = -np.array(np.abs(iw_arr) ** 2, dtype=float)

            for spin in ['up', 'down']:
                self.Delta_iw[spin].data[self.n_iw : 2 * self.n_iw, 0, 0] = -4j / np.sqrt(-iw2_arr[self.n_iw : 2 * self.n_iw] + 1)
                self.extend_to_negative_freq(self.Delta_iw[spin])
        else: 
            self.Delta_iw.zero()
            for l in range(self.nbath_max):
                self.Delta_iw['up'] << self.Delta_iw['up'] + bath.t2_up[l] * inverse(iOmega_n - bath.eps_up[l])
                self.Delta_iw['down'] << self.Delta_iw['down'] + bath.t2_down[l] * inverse(iOmega_n - bath.eps_down[l])

    def solve(self, U, nbath, h, mu, n_m=None, n_kr=None):
        
        if self.method == 'Jacobi':
            self.lnZ = self.cpp_solver.solve(U,
                                        self.Delta_iw['up'].data[self.n_iw : 2 * self.n_iw].squeeze(),
                                        self.Delta_iw['down'].data[self.n_iw : 2 * self.n_iw].squeeze(),
                                        nbath, h, mu)
        elif self.method == 'Lanczos':
            if ((n_m == None or n_kr == None) and not self.warn_seen):
                n_m = 10
                n_kr = 10
                print("Warning: n_m (number of eigenstates) and n_kr (Krylov subspace size) not set.\n"
                      f"Using defaults: n_m = {n_m}, n_kr = {n_kr}.");
                self.warn_seen = True;

            self.lnZ = self.cpp_solver.solve(U,
                                        self.Delta_iw['up'].data[self.n_iw : 2 * self.n_iw].squeeze(),
                                        self.Delta_iw['down'].data[self.n_iw : 2 * self.n_iw].squeeze(),
                                        nbath, n_m, n_kr, h, mu)
        
        self.G_iw['up'].data[self.n_iw : 2 * self.n_iw, 0, 0] = self.cpp_solver.G_up[:]
        self.G_iw['down'].data[self.n_iw : 2 * self.n_iw, 0, 0] = self.cpp_solver.G_down[:]
        self.Sigma_iw['up'].data[self.n_iw : 2 * self.n_iw, 0, 0] = self.cpp_solver.Sigma_up[:]
        self.Sigma_iw['down'].data[self.n_iw : 2 * self.n_iw, 0, 0] = self.cpp_solver.Sigma_down[:]
        
        for spin in ['up', 'down']:
            for gf_positive in [self.G_iw[spin], self.Sigma_iw[spin]]:
                self.extend_to_negative_freq(gf_positive)

        self.bath_parameters.eps_up = self.cpp_solver.eu 
        self.bath_parameters.eps_down = self.cpp_solver.ed 
        self.bath_parameters.t2_up = self.cpp_solver.t2u 
        self.bath_parameters.t2_down = self.cpp_solver.t2d 

    def double_occupancy(self):
        return self.cpp_solver.double_occupancy()

    def mean_occupancy(self):
        return self.cpp_solver.mean_occupancy()
                