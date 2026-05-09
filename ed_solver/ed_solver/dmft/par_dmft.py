from dataclasses import dataclass, field
from pathlib import Path
import warnings

import numpy as np 
from triqs.gf import *

from ed_solver import Bath, Solver
from .dmft_bar import bar_fmt

from h5 import *
from tqdm import tqdm


@dataclass
class DataPoint:
    T: float
    U: float
    t: float = 0.25
    mu: float = 0.0
    h: float = 0.0

    DO: float | None = None      # double occupancy
    n: float | None = None       # mean occupancy 
    n_up: float | None = None    # mean ⟨n↑⟩ 
    n_down: float | None = None  # mean ⟨n↓⟩
    F: float | None = None       # lattice free energy

    converged: bool = False


class ParDMFT:
    def __init__(self, solver, k_mesh, eps):
        self.solver = solver
        self.k_mesh = k_mesh
        self.epsilon = eps

        self.G_loc_iw = self.solver.Sigma_iw.copy()

        self.data = None

    def solve(self, n_bath, alpha, n_loops, atol=1e-7, rtol=1e-5):
        if self.data is None:
            raise RuntimeError(
                "DMFT has not been configured yet.\n"
                "Call .set_point(DataPoint(...)) first.\n"
            )

        converged = False
        diff_up_prev = 0
        diff_down_prev = 0

        bar = tqdm(range(n_loops), 
                    desc=f"Performing DMFT for T = {self.data.T:.4f}, U = {self.data.U:.4f} ...",
                    leave=False)

        for it in bar:

            self.solver.solve(U=self.data.U, nbath=n_bath, h=self.data.h, mu=self.data.mu)
            
            self.G_loc_iw.zero()
            for spin in ['up', 'down']:
                for k in self.k_mesh.values():
                    self.G_loc_iw[spin] += inverse(inverse(self.solver.G_iw[spin]) + 
                                              self.solver.Delta_iw[spin] - self.epsilon(k, self.data.t, spin))
                self.G_loc_iw[spin] /= len(self.k_mesh)
            
            diff_up   = np.max(np.abs(self.G_loc_iw['up'].data - self.solver.G_iw['up'].data))
            diff_down = np.max(np.abs(self.G_loc_iw['down'].data - self.solver.G_iw['down'].data))
            bar.set_postfix(
                du=bar_fmt(diff_up, diff_up_prev),
                dd=bar_fmt(diff_down, diff_down_prev)
            )
            diff_up_prev = diff_up
            diff_down_prev = diff_down

            converged = (np.allclose(self.G_loc_iw['up'].data, self.solver.G_iw['up'].data, atol=atol, rtol=rtol) and
                            np.allclose(self.G_loc_iw['down'].data, self.solver.G_iw['down'].data, atol=atol, rtol=rtol)) 

            if converged:
                self.data.converged = True
                break

            for spin in ['up', 'down']:
                self.solver.Delta_iw[spin] += alpha * (inverse(self.solver.G_iw[spin]) - inverse(self.G_loc_iw[spin]))  

        if not converged:
            diff_up   = np.max(np.abs(self.G_loc_iw['up'].data   - self.solver.G_iw['up'].data))
            diff_down = np.max(np.abs(self.G_loc_iw['down'].data - self.solver.G_iw['down'].data))

            msg = (
                f"\nWarning: DMFT cycle did not converge after {n_loops} iterations.\n"
                f"  T = {self.data.T:.4f}, U = {self.data.U:.4f}\n"
                f"  Last difference spin 'up'   (max |G_loc - G_imp|): {diff_up:.2e}\n"
                f"  Last difference spin 'down' (max |G_loc - G_imp|): {diff_down:.2e}\n"
            )

            warnings.warn(msg, RuntimeWarning, stacklevel=2)

    def set_point(self, data: DataPoint):
        if not isinstance(data, DataPoint):
            raise TypeError("data must be a 'DataPoint' instance\n")

        if abs(self.solver.beta - 1 / data.T) > 1e-10:
            raise ValueError(
                f"Inconsistent inverse temperature (beta):\n"
                f"  Solver: beta = {self.solver.beta:.4g}\n"
                f"  DMFT: 1/T = {1/data.T:.4g}\n"
            )

        self.data = data

    def initial_guess(self, filename, init_point: DataPoint, direction="undirected"):
        if not isinstance(init_point, DataPoint):
            raise TypeError("initial point must be a 'DataPoint' instance\n")

        TUstr = f'T{init_point.T:.4f}U{init_point.U:.4f}'
        filename = Path(filename).resolve()

        try:
            with HDFArchive(str(filename), 'r') as arch:
                if direction not in arch:
                    raise RuntimeError(f"The file '{filename}' does not contain data for direction '{direction}'.")
                else:
                    subgroup = arch[direction]

                if (round(init_point.T, 4), round(init_point.U, 4), True) not in subgroup['points']:
                    raise RuntimeError(
                        f"Initial point with T = {init_point.T:.4f} and U = {init_point.U:.4f} "
                        f"not found in file '{filename}' with direction '{direction}'.\n"
                    ) 
                else:
                    group = subgroup[TUstr]

                eu = group['eps_up']; ed = group['eps_down']
                t2u = group['t2_up']; t2d = group['t2_down']
                bath = Bath(eu, ed, t2u, t2d)
                mu = group['mu']

                self.solver.set_initial_guess(bath)
                self.data = DataPoint(T=init_point.T, U=init_point.U, mu=mu)

        except:
            raise FileNotFoundError(f"File not found: '{filename}'\n")

    def calculate_observables(self, double_occupancy=False, 
                                    mean_occupancy=False,
                                    n_up=False, n_down=False):
        if double_occupancy:
            self.data.DO = self.solver.double_occupancy()
        if mean_occupancy:
            self.data.n = self.solver.mean_occupancy()
        if n_up:
            self.data.n_up = self.solver.n_up()
        if n_down:
            self.data.n_down = self.solver.n_down()

    def export_solver_state(self, filename, direction="undirected"):
        TUstr = f'T{self.data.T:.4f}U{self.data.U:.4f}'
        filename = Path(filename).resolve()

        with HDFArchive(str(filename), 'a') as arch:
            if direction not in arch:
                arch.create_group(direction)
                arch[direction]['points'] = []

            subgroup = arch[direction]   

            current_points = subgroup['points']
            new_point = (np.round(self.data.T, 4), np.round(self.data.U, 4), self.data.converged)
            
            if not any(np.isclose(new_point[0], p[0]) and np.isclose(new_point[1], p[1]) 
                      for p in current_points):
                
                current_points.append(new_point)
                subgroup['points'] = current_points 

            if TUstr not in subgroup:
                subgroup.create_group(TUstr)

            group = subgroup[TUstr]
            group['G-iw']     = self.solver.G_iw
            group['G_loc-iw'] = self.G_loc_iw
            group['Sigma-iw'] = self.solver.Sigma_iw
            group['Delta-iw'] = self.solver.Delta_iw

            group['mu'] = self.data.mu
            group['h'] = self.data.h

            group['eps_up'] = self.solver.bath_parameters.eps_up
            group['eps_down'] = self.solver.bath_parameters.eps_down
            group['t2_up'] = self.solver.bath_parameters.t2_up
            group['t2_down'] = self.solver.bath_parameters.t2_down
            
            group['converged'] = self.data.converged

            group['double_occupancy'] = self.data.DO
            group['mean_occupancy'] = self.data.n
            group['n_up'] = self.data.n_up
            group['n_down'] = self.data.n_down
