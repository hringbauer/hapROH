from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Distutils import build_ext

###
# I distribute the source only,
# to avoid trouble with platform
# specific wheels. Pls build yourself!

### MODIFY IF YOU WANT TO USE CYTHON FOR BUILDING!
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
    name="hapsburg",
    version="0.1a4",  # a means alpha
    author="Harald Ringbauer",
    author_email="harald_ringbauer@hms.harvard.edu",
    description="Calling long ROH in 1240k ancient DNA data using modern reference panel",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hringbauer/HAPSBURG",
    packages=find_packages(),
    #cmdclass=cmdclass,
    ext_modules=extensions,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)