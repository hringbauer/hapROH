Installation
===============

You can install the package using the Package manager pip:

``python3 -m pip install hapROH``

If you have it already installed and want to upgrade to a newer hapROH version you can use:

``python3 -m pip install hapROH --upgrade``

The package distributes source code. The setup.py contains information that should automatically install the package.
For customized installations, find more info in the section c Extension (TODO: figure out how to make a link).

c Extension
************

For performance reasons, the heavy lifting of the algorithms is coded into c methods (cfunc.c). The package is set up so that this "extension" is built during installation. This is done automatically from cfunc.pyx via the package cython (as CYTHON=True in setup.py by default). You can also set CYTHON=False, then the extension is directly compiled from cfunc.c (experimental, not tested on all platforms).


Dependencies
************

The basic dependencies of the package are kept minimal. They are sufficient for the core functions of the algorithms (numpy, pandas, scipy & h5py). If you want to use extended analysis and plotting functionality, there are extra Python packages that you need to install manually (e.g. via pip or conda).

    If you want to use the advanced plotting functionality, you need ``matplotlib``.

    For plotting of maps, you need ``basemap`` (warning: installing can be tricky on some architectures).

    If you want to use the effective population size fitting functionality from ROH output, you require the package ``statsmodels``.

