# Step-by-step guide for building mesas from source on a Mac M1 chip
# [Check Most Up-to-date Notion Page Here](https://www.notion.so/Step-by-step-guide-for-building-mesas-from-source-on-a-Mac-M1-chip-c9804e43bb06481790d6bf714923b47b?pvs=4)

# Step 1: Download newest mesas.py from Github

```bash
# get the newest GitHub repo
git clone https://github.com/charman2/mesas.git
# get into mesas folder
cd <your_root_folder>/mesas 

```

---

# Step 2: Config a virtual environment for mesas.py

## 1. Create a conda env

```bash
# make sure to load conda module first if on a HPC
conda create -n <your-mesas-source-env-name>
conda activate <your-mesas-source-env-name>
```

## 2. Setup a fortran compile

### a) Routine

```bash
mkdir mesas/sas/cdflib90/_build
conda install fortran-compiler cmake -c conda-forge
cmake -S mesas/sas/cdflib90 -B mesas/sas/cdflib90/_build
cmake --build mesas/sas/cdflib90/_build
```

### b) If you are updating mesas.py/rebuilding fortran compiler, ALWAYS run the following command line and then repeat a)

```bash
# Otherwise ignore this part
rm -rf mesas/sas/cdflib90/_build
# followed by the 4 lines in a)
```

### c) If failed at the 3rd line of a): `cmake -S mesas/sas/cdflib90 -B mesas/sas/cdflib90/_build`, with a problem in Fortran compile

```bash
# just like b), you need to remove old build
rm -rf mesas/sas/cdflib90/_build
mkdir mesas/sas/cdflib90/_build
# update gfortran and gcc
conda uninstall gfortran
conda install -c conda-forge gfortran
brew reinstall gcc # note you need to install homebrew first
```

Now, proceed to a) to continue the routine.  Note, I always open another terminal to continue the routine process just to be safe. Also, check your current working directory if running in new terminal. 

## 3.  Build mesas.py as a Python package

### In case you did not config python for the conda virtual env, i.e. having this issue: `command not found: pip`:

```bash
conda install python 
```

### Install package (make sure you are in `<your_root_folder>/mesas`):

```bash
pip install .
```

To debug and get verbose out, try: `pip install -v .` 

---

# Step 3: Check if https://github.com/charman2/mesas is running properly

## a) Make sure you are not in the folder where you downloaded https://github.com/charman2/mesas

```bash
# you can simply do this if you have been strictly following my guide so far
cd ..
```

## b) Get into Python

```bash
# simply type python
python
# If you are in python, the command line will look like >>>
# to exit python, type exit()
```

## c) In python, run:

```python
import mesas
# If this failed, you did not successfully install mesas
from mesas.sas.model import Model
# If this failed, two possibilities:
# 1. you are in ./mesas/, see Step 3 a)
# 2. your fortran compiler failed, see Step 2
```