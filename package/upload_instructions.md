# Instructions how to upload package to PYPI

Summary of update and upload instructions from https://packaging.python.org/tutorials/packaging-projects/

### Go to base folder
cd /project2/jnovembre/hringbauer/HAPSBURG/package

- On Chicago cluster:  
module load python

- Or on O2:
envpython37

### For local test install (e.g. for unit tests), run from ./package/ folder:
pip3 install --user ./


### Run Tests of Expected Behavior
Use `/Notebooks/Tests/roh_caller_test_chicago.ipynb` to run tests of expected behavior of hapROH
Use `/Notebooks/Tests/hapCon_test_chicago.ipynb` to run tests of hapCon.

### Create the Source Package 
Update version in setup.py to next version number

### Update setuptools. 
Delete previous ./dist/* (alternatively be specific below what to upload):  

rm ./dist/*

Run the setup file:
python3 setup.py sdist

### Upload to the Sources (copy into shell, to interactively do it!)
### For full PyPi server
python3 -m twine upload dist/* 
### Alternatively: Upload on test server (for testing)
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/* 


## To test whether extensions builds
python3 setup.py build_ext --inplace

# Further Documentation 
### To install via pip:
Follow instructions on pypi site of `hapROH`.

### for packaging: 
https://packaging.python.org/tutorials/packaging-projects/

### For C extension (here cython):

### for version numbers:
https://www.python.org/dev/peps/pep-0440/
