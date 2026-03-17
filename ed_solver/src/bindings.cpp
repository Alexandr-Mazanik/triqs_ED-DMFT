#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include "solver.h"
#include "solverLanczos.h"

namespace py = pybind11;
using complex = std::complex<double>;

PYBIND11_MODULE(_ed_solver, m, py::mod_gil_not_used()) {
  m.doc() = "Exact Diagonalization impurity solver for DMFT";

  py::class_<Solver>(m, "EDSolverCore")
    .def(py::init<double, int, int, int>(),
       py::arg("Beta"), py::arg("Nbath_max"), py::arg("Nw"), py::arg("NwED") = 20,
       "Create the ED solver instance.")

    .def("solve", [](Solver& self, 
                     double U,
             	       const py::array_t<complex>& Delta_up,
                     const py::array_t<complex>& Delta_down,
                     int Nbath, double h_loc, double mu_loc) {
      
      	if (Delta_up.ndim() != 1 ||  Delta_up.size() != self.Nw)
            throw py::value_error("Delta_up must be a 1D array of size Nw");
        
        if (Delta_down.ndim() != 1 || Delta_down.size() != self.Nw) 
            throw py::value_error("Delta_down must be a 1D array of size Nw");

          return self.init(U,
                           const_cast<complex * >(Delta_up.data()),
                   		     const_cast<complex * >(Delta_down.data()),
                           Nbath, h_loc, mu_loc);

      },
      
      py::arg("U"), py::arg("Delta_up"), py::arg("Delta_down"),
      py::arg("N_bath"), py::arg("h_loc"), py::arg("mu_loc"),
      "Solve the impurity problem for given hybridization functions.\n"
      "U: interaction term \n"
      "Delta_up/down: complex arrays of size Nw (positive frequencies only)\n"
      "Nbath: number of bath sites to use (≤ Nbath_max)\n"
      "h_loc: local AF field\n"
      "mu_loc: chemical potential")

    .def("double_occupancy", &Solver::double_occupancy,
      "Returns the double occupancy ⟨n↑ n↓⟩")

    .def("mean_occupancy", &Solver::mean_occupancy,
      "Returns the mean occupancy ⟨n⟩")

    .def_property_readonly("eu", [](const Solver& s) {
      return py::array_t<double>(s.Nbath_max, s.eu);
    }, "Bath on-site energy levels for spin up (shape: Nbath_max,)")

    .def_property_readonly("ed", [](const Solver& s) {
      return py::array_t<double>(s.Nbath_max, s.ed);
    }, "Bath on-site energy levels for spin down (shape: Nbath_max,)")

    .def_property_readonly("t2u", [](const Solver& s) {
      return py::array_t<double>(s.Nbath_max, s.t2u);
    }, "Coupling strengths between impurity and bath sites for spin up (shape: Nbath_max,)")

    .def_property_readonly("t2d", [](const Solver& s) {
      return py::array_t<double>(s.Nbath_max, s.t2d);
    }, "Coupling strengths between impurity and bath sites for spin down (shape: Nbath_max,)")

    .def_property_readonly("G_up", [](const Solver& s) {
      return py::array_t<std::complex<double>>(s.Nw, s.gu);
    }, "Impurity Green's function (spin up) on positive Matsubara frequencies (shape: (Nw,))")

    .def_property_readonly("G_down", [](const Solver& s) {
      return py::array_t<std::complex<double>>(s.Nw, s.gd);
    }, "Impurity Green's function (spin G_down) on positive Matsubara frequencies (shape: (Nw,))")

    .def_property_readonly("Sigma_up", [](const Solver& s) {
      return py::array_t<std::complex<double>>(s.Nw, s.sigmau);
    }, "Impurity self-energy (spin up) on positive Matsubara frequencies (shape: (Nw,))")

    .def_property_readonly("Sigma_down", [](const Solver& s) {
      return py::array_t<std::complex<double>>(s.Nw, s.sigmad);
    }, "Impurity self-energy (spin down) on positive Matsubara frequencies (shape: (Nw,))")

    .def_readonly("Nw", &Solver::Nw)
    .def_readonly("Beta", &Solver::Beta)
    .def_readonly("Nbath_max", &Solver::Nbath_max);


  py::class_<SolverLanczos>(m, "EDLSolverCore")
    .def(py::init<double, int, int, int>(),
      py::arg("Beta"), py::arg("Nbath_max"), py::arg("Nw"), py::arg("NwED") = 20,
      "Create the ED solver (Lanczos algorithm) instance.")

    .def("solve", 
      &SolverLanczos::init,
      py::arg("U"), py::arg("Delta_up"), py::arg("Delta_down"),
      py::arg("N_bath"), py::arg("Nm"), py::arg("Nkr"), 
      py::arg("h_loc") = 0.0, py::arg("mu_loc") = 0.0,
      "Solve the impurity problem for given hybridization functions.\n"
      "U: interaction term \n"
      "Delta_up/down: complex arrays of size Nw (positive frequencies only)\n"
      "Nbath: number of bath sites to use (≤ Nbath_max)\n"
      "Nm: number of lowest-energy eigenstates to compute\n"
      "Nkr: size of the Krylov subspace\n"
      "h_loc: local AF field\n"
      "mu_loc: chemical potential")

    .def_property_readonly("G_up",
      [](const SolverLanczos& s) { return s.gu; },
      "Impurity Green's function (spin up) on positive Matsubara frequencies (shape: (Nw,))")

    .def_property_readonly("G_down",
      [](const SolverLanczos& s) { return s.gd; },
      "Impurity Green's function (spin down) on positive Matsubara frequencies (shape: (Nw,))")

    .def_property_readonly("Sigma_up",
      [](const SolverLanczos& s) { return s.sigmau; },
      "Impurity self-energy (spin up) on positive Matsubara frequencies (shape: (Nw,))")

    .def_property_readonly("Sigma_down",
      [](const SolverLanczos& s) { return s.sigmad; },
      "Impurity self-energy (spin down) on positive Matsubara frequencies (shape: (Nw,))")

    .def_readonly("Nw", &SolverLanczos::Nw)
    .def_readonly("Beta", &SolverLanczos::Beta)
    .def_readonly("Nbath_max", &SolverLanczos::Nbath_max);
}
