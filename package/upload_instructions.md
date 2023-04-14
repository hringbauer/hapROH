# Instructions how to upload package to PYPI

This file is based on update and upload instructions from https://packaging.python.org/tutorials/packaging-projects/

## First steps

- On Leipzig Archgen Cluster:
cd /mnt/archgen/users/hringbauer/git/hapROH/package/

%%% LEGACY now (Leipzig is primary development server)
%- On Chicago cluster:  
%cd /project2/jnovembre/hringbauer/HAPSBURG/package
%module load python

- On Harvard O2:
envpython37

### For a local test install (e.g. for unit tests), run from ./package/ folder:
pip3 install --user ./


### Run Tests of Expected Behavior
- Use `../Notebooks/Tests/hapROH_test_leipzig.ipynb` to run tests of expected behavior of hapROH
- Use `../Notebooks/Tests/hapCon_test_leipzig.ipynb` to run tests of hapCon.

### Create the Source Package 
Update version in `./setup.py` to next version number and update `./change_log.md`

### Update setuptools. 
Delete previous `./dist/*` (alternatively be specific below what to upload):  

rm ./dist/*

### Run the setup file:
python3 setup.py sdist

### Upload to the Sources (copy into shell, to interactively do it!)
### For full PyPi server (user name `hringbauer`)
python3 -m twine upload dist/* 
### Alternatively: Upload on test server (for testing)
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/* 


## To test whether extensions builds
python3 setup.py build_ext --inplace

# Further Documentation 
### To install via pip:
Follow instructions on pypi site of `hapROH`.

`python3 -m pip install --user --upgrade --no-deps --force-reinstall hapROH`

### for packaging: 
https://packaging.python.org/tutorials/packaging-projects/

### For C extension (here cython):

### for version numbers:
https://www.python.org/dev/peps/pep-0440/
