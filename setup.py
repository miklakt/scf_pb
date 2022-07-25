from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

__version__ = "0.0.1"

ext_modules = [
    Pybind11Extension("scf_pb",
        ["scf_pb/src/pybind.cpp"],
        define_macros = [('VERSION_INFO', __version__)],
        ),
]

setup(
    name="scf_pb",
    version=__version__,
    author="Mikhail Laktionov",
    author_email="miklakt@gmail.com",
    #url="https://github.com/pybind/python_example",
    description="A test project using pybind11",
    long_description="",
    ext_modules=ext_modules,
    #extras_require={"test": "pytest"},
    cmdclass={"build_ext": build_ext},
    #zip_safe=False,
    python_requires=">=3.6",
)