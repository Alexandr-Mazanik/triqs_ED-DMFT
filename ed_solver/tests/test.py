import numpy as np
import ed_solver


def main():
	Beta = 5 
	U = 6 
	h = 0.01 
	mu = 0
	Nbath = 3
	Nw = 200

	solver = ed_solver.EDSolverCore(Beta, Nbath, Nw)

	Delta_up = np.full(solver.Nw, -0.5j, dtype=complex)
	Delta_down = np.full(solver.Nw, -0.5j, dtype=complex)

	print("Starting ED solve...\n")
	lnZ = solver.solve(U, Delta_up, Delta_down, Nbath, h, mu)

	print(f"lnZ_imp = {np.round(lnZ, 4)}")
	print(f"G_up(iw_0) = {np.round(solver.G_up[0], 4)}")
	print(f"Sigma_up(iw_0) = {np.round(solver.Sigma_up[0], 4)}")


if __name__ == "__main__":
	main()
