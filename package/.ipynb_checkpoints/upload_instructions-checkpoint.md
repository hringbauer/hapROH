# Instructions how to upload package to PYPI

Summaring up update and upload instructions from https://packaging.python.org/tutorials/packaging-projects/

### Go to base folder
cd /project2/jnovembre/hringbauer/HAPSBURG/package

On Chicago cluster:  
module load python/3.7.0

### Create the Source Package 
Update version in setup.py to next version number (if wanted)

### (update setuptools if needed, delete previous dist/* or be specific below)
python3 setup.py sdist

### Upload to the Sources (copy into shell, to interactively do it!)
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

## To test whether extensions builds
python3 setup.py build_ext --inplace

# Further Documentation 
### for packaging: 
https://packaging.python.org/tutorials/packaging-projects/

### For C extension (here cython):

### for version numbers:
https://www.python.org/dev/peps/pep-0440/
