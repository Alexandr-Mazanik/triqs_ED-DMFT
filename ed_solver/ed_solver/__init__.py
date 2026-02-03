from ._ed_solver import EDSolverCore
from ._ed_solver import EDLSolverCore

import warnings

try:
    from .triqs_wrapper import Solver, Bath
    from .dmft import ParDMFT
    _triqs_available = True
except ImportError:
    _triqs_available = False
    Solver = None

    if not hasattr(warnings, '_ed_solver_triqs_warning_shown'):
        msg = (
                "TRIQS is not installed.\n"
                "The 'Solver' class (TRIQS-compatible wrapper) is unavailable.\n"
                "The DMFT implementations are unavailable.\n"
                "Only the original EDSolverCore and EDLSolverCore classes can be used.\n"
        )

        warnings.warn(msg, UserWarning, stacklevel=2)
        warnings._ed_solver_triqs_warning_shown = True

__version__ = "0.1.0"
__all__ = [
    "EDSolverCore", 
    "EDLSolverCore", 
    "Solver",
    "Bath",
    "ParDMFT"
]
