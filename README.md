# Scripts and Data for the "Nanobody Article"

This repository contains scripts and data for the article entitled 
    _Investigation and Assessment of a Machine Learning and 
    Molecular Dynamics Assisted anti-SARS-CoV2 Nanobody Design 
    Workflow_.

## Requirements

These scripts were run using:

- Python 3.11.7
- NumPy 1.26.3
- Pandas 2.2.0 
- BioPython 1.84.dev0
- Matplotlib 3.8.2
- SciPy 1.12.0
- Scikit-learn 1.4.0

## Uses

The contents are the following:
- `run_data_processors` contains basic scripts that convert the ML webserver inputs
 and outputs in order to be compatible with each other. The number at the beginning 
 of the filenames indicate their ideal execution order. Put every script into the same
 directory with the corresponding files to be processed and run them without any arguments.
- `data` contains the ML run results, as well as the processed versions of each ML input
 and output file.
- `pdb_files` contain the corresponding structure files for the ML algorithms. The MD
 simulations have been started from the FoldX mutated versions of these files.
- `main.py` will do all ML-related data analysis using the files from `data`. 