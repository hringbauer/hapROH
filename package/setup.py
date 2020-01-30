from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Distutils import build_ext

###
# I distribute the source only,
# to avoid trouble with platform
# specific wheels. Pls build yourself!

USE_CYTHON = True   

# Code for Cython / C Extension

### First attempt
ext = '.pyx' if USE_CYTHON else '.c'
extensions = [Extension("cfunc", ["hapsburg/cfunc" + ext])]
if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

### 2nd attempt see https://stackoverflow.com/questions/4505747/how-should-i-structure-a-python-package-that-contains-cython-code
#cmdclass = {}
#extensions = []
#if USE_CYTHON:
#    extensions += [
#        Extension("cfunc", ["cfunc.pyx"]),
#    ]
#    cmdclass.update({'build_ext': build_ext})
  
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="hapsburg",
    version="0.1a3",  # a means alpha
    author="Harald Ringbauer",
    author_email="harald_ringbauer@hms.harvard.edu",
    description="Identifction of long ROH in 1240k ancient DNA data using modern reference panel",
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
