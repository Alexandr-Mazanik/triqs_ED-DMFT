from ._ed_solver import EDSolverCore
from ._ed_solver import EDLSolverCore

import warnings

try:
    from .triqs_wrapper import Solver
    _triqs_available = True
except ImportError:
    _triqs_available = False
    Solver = None

    if not hasattr(warnings, '_ed_solver_triqs_warning_shown'):
        warnings.warn(
            "TRIQS is not installed. The 'Solver' class (TRIQS-compatible wrapper) is unavailable.\n Only the original EDSolverCore class can be used.",
            UserWarning
        )
        warnings._ed_solver_triqs_warning_shown = True

__version__ = "0.1.0"
__all__ = ["EDSolverCore", "EDLSolverCore", "Solver"]
