import numpy as np

from triqs.gf import *
import ed_solver as ed
from h5 import *

from tqdm import tqdm


def main():
    t = 1.0                              # hopping
    beta = 10.0                          # 1/T, where T is the temperature in units of t
    mu = 0                               # half-filling
    h = 0.01                             # paramagnetic field   
    U_list = np.arange(0.0, 10.0, 1.0)   # list of interaction parameters 

    n_loops_max = 200                    # maximum number of DMFT iterations

    n_bath = 3                           # number of bath levels used in ED
    alpha = 0.2                          # hybridization update parameter


    # Construct the impurity solver
    S = ed.Solver(beta=beta, nbath_max=n_bath, method='Jacobi')
    delta_old = S.Delta_iw.copy()

    for U in tqdm(U_list, 
                  desc="Scan over U...", 
                  position=0, leave=True):

        # First guess for Delta
        S.G_iw << SemiCircular(2*t)
        S.Delta_iw << t**2 * S.G_iw
        delta_old << S.G_iw

        # DMFT loop with self-consistency
        for i in tqdm(range(n_loops_max), 
                      desc=f"Performing DMFT calculation for U = {U:.2f}", 
                      position=1, leave=False):

            # Solve the impurity problem
            S.solve(U=U, nbath=n_bath, h=h, mu=mu)
            
            # self-consistency
            S.Delta_iw += alpha * (t**2 * S.G_iw - delta_old)

            # convergence check
            converged = np.allclose(delta_old['up'].data, S.Delta_iw['up'].data, rtol=1e-5)

            delta_old << S.Delta_iw

            if converged or i == n_loops_max - 1:
                export_state(S, U)
            
            if converged:
                break
            else:
                if i == n_loops_max - 1:
                    print(f"\nReached max iterations ({n_loops_max}) without convergence (U = {U:.2f})\n")


def export_state(S, U):
    Ustr = f'U{U:.2f}'

    with HDFArchive("data/ed_bethe.h5", 'a') as arch:
        if Ustr not in arch:
            arch.create_group(Ustr)
        group = arch[Ustr]

        group['G-res']     = S.G_iw
        group['Sigma-res'] = S.Sigma_iw
        group['Delta-res'] = S.Delta_iw

        group['double_occupancy'] = S.double_occupancy()


if __name__ == "__main__":
    main()
