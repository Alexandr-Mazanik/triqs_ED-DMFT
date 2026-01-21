import os
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension
from glob import glob

conda_prefix = os.environ.get("CONDA_PREFIX")
if conda_prefix is None:
    raise RuntimeError(
        "This package must be built inside an activated conda/mamba environment "
        "with 'eigen' and 'spectra-cpp' installed via conda-forge."
    )

include_dirs = [
    os.path.join(conda_prefix, "include"),         
    os.path.join(conda_prefix, "include", "eigen3"),
]

ext_modules = [
    Pybind11Extension(
        "ed_solver._ed_solver",
        sorted(glob("src/*.cpp")),
        cxx_std=17,
        include_dirs=include_dirs,
        extra_compile_args=["-Ofast", "-Wall", "-march=native"],
    ),
]

setup(
    name="ed_solver",
    version="0.1.0",
    description="Exact Diagonalization impurity solver for DMFT",
    ext_modules=ext_modules,
    packages=["ed_solver"],       
    package_dir={"": "."},   
    zip_safe=False,
)
