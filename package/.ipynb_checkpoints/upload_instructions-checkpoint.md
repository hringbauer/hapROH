# Instructions how to upload package to PYPI

Summaring the following:
Follow update and upload instructions on https://packaging.python.org/tutorials/packaging-projects/

### Go to base folder
cd /project2/jnovembre/hringbauer/HAPSBURG/package

On Chicago cluster:  
module load python/3.7.0

### Create the Source Package (update setuptools if needed)
python3 setup.py sdist

### Upload to the Sources (copy into shell, to interactively do it!)
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

# Further Documentation 
### for packaging: 
https://packaging.python.org/tutorials/packaging-projects/

### For C extension (here cython):

### for version numbers:
https://www.python.org/dev/peps/pep-0440/