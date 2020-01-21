from setuptools import setup, find_packages
from setuptools.extension import Extension


# Eventually change Name to setup_pypi.py
# Documentation for packaging at
# https://packaging.python.org/tutorials/packaging-projects/
# Documentation for C extensino (here CYython):
# For version numbers:
# (https://www.python.org/dev/peps/pep-0440/)
###

USE_CYTHON = False   # Whether to use CYthon for extension

# Code for Cython / C Extension
ext = '.pyx' if USE_CYTHON else '.c'
extensions = [Extension("cfunc", ["cfunc" + ext])]
if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="hapsburg",
    version="0.1a1",  # a means alpha
    author="Harald Ringbauer",
    author_email="hringbauer@uchicago.edu",
    description="Identifction of long ROH in 1240k ancient DNA data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hringbauer/HAPSBURG",
    packages=setuptools.find_packages(),
    ext_modules=extensions
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
