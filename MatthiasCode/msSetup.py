"""
Setup Cython Code for Forward Backward and Viterbi Path
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

from distutils.core import setup
from Cython.Build import cythonize

# python3 msSetup.py build_ext --inplace
# cythonize -a -i yourmod.pyx

setup(
    ext_modules = cythonize("msFun.pyx")
)
