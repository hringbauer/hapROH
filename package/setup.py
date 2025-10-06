from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

###
# We distribute the source only,
# to avoid trouble with platform


### Code for Cython / C Extension
USE_CYTHON = True 
ext = '.pyx' if USE_CYTHON else '.c'
extensions = [Extension("hapsburg.cfunc", ["hapsburg/cfunc" + ext], include_dirs=[np.get_include()])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

###
setup(
    ext_modules=extensions,
)
