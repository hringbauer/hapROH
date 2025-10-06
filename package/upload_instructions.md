# Instructions on how to upload a package to PYPI
This file is based on update and upload instructions from `https://packaging.python.org/tutorials/packaging-projects/`
October 2025: Switched from setuptools to pyproject.toml setup. 

I distribute the source - and no wheels (for now). Because the Cython extension makes it not pure Python, building the wheel is not trivial.

## First steps
- On Leipzig HPC Cluster (primary source to create package since 2025):
Activate Python environment and go to the hapROH package folder:

pyenvhpc312
cd /mnt/archgen/users/hringbauer/git/hapROH/package/

%%% LEGACY as of 2025
%- On Leipzig Archgen Cluster (virtualenv does not work for installation as of 2025...)
%cd /mnt/archgen/users/hringbauer/git/hapROH/package/

%%% LEGACY as of 2024 (Leipzig is now the primary development server)
%- On Chicago cluster:  
%cd /project2/jnovembre/hringbauer/HAPSBURG/package
%module load python

%- On Harvard O2:
%envpython37

### For a local test install (e.g., for unit tests), run from ./package/ folder:
pip3 install ./

Add flag  `--user` if not in a virtual environment.

### Run Tests of Expected Behavior
- Use `../Notebooks/Tests/hapROH_test_leipzig.ipynb` to run tests of the expected behavior of hapROH
- Use `../Notebooks/Tests/hapCon_test_leipzig.ipynb` to run tests of hapCon.

# Updated in 2025 to pyproject.toml architecture
### Create the Source Package 
Update version in `./pyproject.toml` to next version number and update `./change_log.md`

### Clean up prior packages
Delete previous `./dist/*`:

rm ./dist/*

### Build package (only source no wheel)
python3 -m build --sdist

%%% Legacy (<2025)
%### Create the Source Package 
%Update version in `./setup.py` to next version number and update `./change_log.md`

%### Update setuptools. 
%Delete previous `./dist/*` (alternatively, be specific below what to upload):  

%rm ./dist/*

%### Run the setup file:
%python3 setup.py sdist

### Upload to PyPi
### For full PyPi server
python3 -m twine upload dist/* 

### [Alternatively] Upload to the test server (for testing)
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/* 

%## To test whether extensions builds
%python3 setup.py build_ext --inplace

### Test install
Follow the instructions on the PyPI site of `hapROH`.

Re-install to check whether the new hapROH version installs properly:

python3 -m pip install --upgrade --no-deps --force-reinstall hapROH

Add the `--user` flag on Python settings where packages need to go to the user. 

# Ressources for further reading
### for packaging: 
https://packaging.python.org/tutorials/packaging-projects/

### for version numbers:
https://www.python.org/dev/peps/pep-0440/
