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

For performance reasons, the heavy lifting of the algorithms is coded into c methods (cfunc.c). This "extension" is built via cython from cfunc.pyx This should be done automatically via the package cython (as CYTHON=True in setup.py by default).
You can also set CYTHON=False, then the extension is compiled from cfunc.c directly (experimental, not tested on all platforms).


Dependencies
************

The basic dependencies are kept minimal. They are sufficient for the core functions of ROH calling (numpy, pandas, scipy & h5py). If you want to use extended analysis and plotting functionality, there are extra Python packages that you need to install (e.g. via pip or conda).

    If you want to use the advanced plotting functionality, you need ``matplotlib``.

    For plotting of maps, you need ``basemap`` (warning: installing can be tricky on some architectures).

    If you want to use the effective population size fitting functionality from ROH output, you require the package ``statsmodels``.

