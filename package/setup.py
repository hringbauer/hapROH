from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Distutils import build_ext

###
# We distribute the source only,
# to avoid trouble with platform

### MODIFY IF YOU WANT TO USE CYTHON FOR BUILDING:
USE_CYTHON = True   

# Code for Cython / C Extension
ext = '.pyx' if USE_CYTHON else '.c'
extensions = [Extension("hapsburg.cfunc", ["hapsburg/cfunc" + ext])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)
  
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="hapROH",
    version="0.52",  # a means alpha 0.3a4
    packages=find_packages(),
    ext_modules=extensions,
    python_requires='>=3.6',
    install_requires=['numpy', 'pandas', 'scipy', 'h5py', 'psutil', 'numdifftools', 'cython', 'matplotlib', 'pysam'],
)
