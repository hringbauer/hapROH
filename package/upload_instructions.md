# Instructions how to upload package to PYPI

Summary of update and upload instructions from https://packaging.python.org/tutorials/packaging-projects/

### Go to base folder
cd /project2/jnovembre/hringbauer/HAPSBURG/package

On Chicago cluster:  
module load python/3.7.0

### Create the Source Package 
Update version in setup.py to next version number (if wanted)

### Update setuptools. 
If needed, delete previous dist/* (or be specific below what to upload)

python3 setup.py sdist

### Upload to the Sources (copy into shell, to interactively do it!)
### 1) For test server
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/* 
### 2) For full PyPi server
python3 -m twine upload dist/* 

## To test whether extensions builds
python3 setup.py build_ext --inplace

# Further Documentation 
### To install via pip:
Look at instructions on pypi site

### for packaging: 
https://packaging.python.org/tutorials/packaging-projects/

### For C extension (here cython):

### for version numbers:
https://www.python.org/dev/peps/pep-0440/
