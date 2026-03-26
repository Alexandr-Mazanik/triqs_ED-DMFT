from dataclasses import dataclass, field
from pathlib import Path
import warnings

import numpy as np 
from triqs.gf import *

from ed_solver import Bath, Solver
from .par_dmft import DataPoint
from .dmft_bar import bar_fmt

from h5 import *
from tqdm import tqdm


@dataclass
class MagneticDataPoint(DataPoint):
    B: float = 0.0
    q: int = 1


class MagneticParDMFT:
    def __init__(self, solver, k_mesh, eps):
        self.solver = solver
        self.k_mesh = k_mesh
        self.epsilon = eps

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

        G_loc_iw = self.solver.Sigma_iw.copy()
        G_magnetic_iw = BlockGf(mesh=self.solver.Sigma_iw.mesh, gf_struct=[('up', self.data.q), ('down', self.data.q)])

        with tqdm(range(n_loops),
                  desc=f"DMFT for T = {self.data.T:.4f}, U = {self.data.U:.4f}, B = {self.data.B:.3f}, q = {self.data.q}",
                  leave=False) as bar:

            for it in bar:

                self.solver.solve(U=self.data.U, nbath=n_bath, h=self.data.h, mu=self.data.mu)

                G_magnetic_iw.zero()
                for spin in ['up', 'down']:
                    M = Gf(mesh=self.solver.Sigma_iw.mesh, target_shape=(self.data.q, self.data.q))
                    m = M.copy()

                    m.data[:, np.arange(self.data.q), np.arange(self.data.q)] = (inverse(self.solver.G_iw[spin]) + 
                                                                                  self.solver.Delta_iw[spin]).data[:, 0, 0, None]
                    for k in self.k_mesh.values():
                        M.data[:] = m.data[:] - self.epsilon(k, self.data.t, self.data.B, self.data.q, spin)
                        G_magnetic_iw[spin] += inverse(M)

                    G_magnetic_iw[spin] /= len(self.k_mesh)

                    G_loc_iw[spin].data[:, 0, 0] = np.trace(G_magnetic_iw[spin].data, axis1=1, axis2=2) / self.data.q

                diff_up   = np.max(np.abs(G_loc_iw['up'].data - self.solver.G_iw['up'].data))
                diff_down = np.max(np.abs(G_loc_iw['down'].data - self.solver.G_iw['down'].data))
                bar.set_postfix(
                    du=bar_fmt(diff_up, diff_up_prev),
                    dd=bar_fmt(diff_down, diff_down_prev)
                )
                diff_up_prev = diff_up
                diff_down_prev = diff_down

                converged = (np.allclose(G_loc_iw['up'].data, self.solver.G_iw['up'].data, atol=atol, rtol=rtol) and
                            np.allclose(G_loc_iw['down'].data, self.solver.G_iw['down'].data, atol=atol, rtol=rtol))  

                if converged:
                    self.data.converged = True
                    break

                for spin in ['up', 'down']:
                    self.solver.Delta_iw[spin] += alpha * (inverse(self.solver.G_iw[spin]) - inverse(G_loc_iw[spin]))  

        if not converged:
            diff_up   = np.max(np.abs(G_loc_iw['up'].data   - self.solver.G_iw['up'].data))
            diff_down = np.max(np.abs(G_loc_iw['down'].data - self.solver.G_iw['down'].data))

            msg = (
                f"\nWarning: DMFT cycle did not converge after {n_loops} iterations.\n"
                f"  T = {self.data.T:.4f}, U = {self.data.U:.4f}\n"
                f"  B = {self.data.B:.3f}, q = {self.data.q}\n"
                f"  Last difference spin 'up'   (max |G_loc - G_imp|): {diff_up:.2e}\n"
                f"  Last difference spin 'down' (max |G_loc - G_imp|): {diff_down:.2e}\n"
            )

            warnings.warn(msg, RuntimeWarning, stacklevel=2)

    def set_point(self, data: MagneticDataPoint):
        if not isinstance(data, MagneticDataPoint):
            raise TypeError("data must be a 'MagneticDataPoint' instance\n")

        if abs(self.solver.beta - 1 / data.T) > 1e-10:
            raise ValueError(
                f"Inconsistent inverse temperature (beta):\n"
                f"  Solver: beta = {self.solver.beta:.4g}\n"
                f"  DMFT: 1/T = {1/data.T:.4g}\n"
            )

        self.data = data

    def initial_guess(self, filename, init_point: DataPoint | MagneticDataPoint, direction="undirected"):

        if isinstance(init_point, DataPoint):
            point_str = f'T{init_point.T:.4f}U{init_point.U:.4f}'
        elif isinstance(init_point, MagneticDataPoint):
            point_str = f'B{self.data.B:.4f}'
        else:
            raise TypeError("initial point must be a 'DataPoint' or 'MagneticDataPoint' instance\n")

        filename = Path(filename).resolve()

        try:
            with HDFArchive(str(filename), 'r') as arch:
                if direction not in arch:
                    raise RuntimeError(f"The file '{filename}' does not contain data for direction '{direction}'.")
                else:
                    subgroup = arch[direction]

                if isinstance(init_point, DataPoint):
                    if (round(init_point.T, 4), round(init_point.U, 4), True) not in subgroup['points']:
                        raise RuntimeError(
                            f"Initial point with T = {init_point.T:.4f} and U = {init_point.U:.4f} "
                            f"not found in file '{filename}' with direction '{direction}'.\n"
                        ) 
                    else:
                        group = subgroup[point_str]
                else:
                    if (np.round(init_point.B, 4), True) not in subgroup['points']:
                        raise RuntimeError(
                            f"Initial point with B = {init_point.B:.4f} "
                            f"not found in file '{filename}' with direction '{direction}'.\n"
                        ) 
                    else:
                        group = subgroup[point_str]

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
        Bstr = f'B{self.data.B:.4f}'
        filename = Path(filename).resolve()

        with HDFArchive(str(filename), 'a') as arch:
            if direction not in arch:
                arch.create_group(direction)
                arch[direction]['points'] = []

            subgroup = arch[direction]   

            current_points = subgroup['points']
            new_point = (np.round(self.data.B, 4), self.data.converged)
            
            if not any(np.isclose(new_point[0], p[0]) for p in current_points):
                
                current_points.append(new_point)
                subgroup['points'] = current_points 

            if Bstr not in subgroup:
                subgroup.create_group(Bstr)

            group = subgroup[Bstr]
            group['G-iw']     = self.solver.G_iw
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
