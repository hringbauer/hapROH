Installation
===============

You can install the package using the Package manager pip:

cpython3 -m pip install hapROH``

If you have it already installed and want to upgrade to a newer hapROH version you can use:

``python3 -m pip install hapROH --upgrade``

The package distributes source code. The setup.py contains information that should automatically install the package.
For customized installations, find more info in the section c Extension (TODO: figure out how to make a link).

c Extension
************

For performance reasons, the heavy lifting of the algorithms is coded into c methods (cfunc.c). The package is set up so that this "extension" is built during installation. This is done automatically from cfunc.pyx via the package cython (as CYTHON=True in setup.py by default). You can also set CYTHON=False, then the extension is directly compiled from cfunc.c (experimental, not tested on all platforms).


Dependencies
************

The basic dependencies of the package are kept minimal. They are sufficient for the core functions of the algorithms. When ``hapROH`` is installed, the following dependent Python packages should be automatically installed without any action on your part: 
* ``numpy`` for calculaions with numerical arrays at C speed 
* ``pandas`` for handling databases and tables at C speed 
* ``scipy`` for statistical operations at C speed
* ``h5py`` for handling hdf5, a file format with partial I/O
* ``psutil`` for process monitoring
* ``numdifftools`` for automatic numerical differentiation


If you want to use extended analysis and plotting functionality, there are extra Python packages that you need to install manually (e.g. via pip or conda).

    If you want to use the advanced plotting functionality, you need ``matplotlib``.

    For plotting of maps, you need ``basemap`` (warning: installing can be tricky on some architectures).

    If you want to use the effective population size fitting functionality from ROH output, you require the package ``statsmodels``.
    
Upgrading    
************
If you installed a stable ``hapROH`` release via pip and wish to upgrade to the latest stable release, you can do so by adding ``--upgrade``:

``pip install --upgrade hapROH``
