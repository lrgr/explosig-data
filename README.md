[![Build Status](https://travis-ci.org/lrgr/explosig-data.svg?branch=master)](https://travis-ci.org/lrgr/explosig-data)
[![PyPI](https://img.shields.io/pypi/v/explosig-data)](https://pypi.org/project/explosig-data/)

## ExploSig Data

Helpers for processing mutation data into standard formats originally developed for the [ExploSig](https://github.com/lrgr/explosig) family of tools.

- [Documentation](https://lrgr.github.io/explosig-data/)

### Installation

```sh
pip install explosig-data
```

### Example 

With raw SSM/MAF file from ICGC or TCGA:

```python
>>> import explosig_data as ed

>>> # Step 1: Process into the ExploSig "standard format":
>>> data_container = ed.standardize_ICGC_ssm_file('path/to/ssm.tsv') # if ICGC
>>> data_container = ed.standardize_TCGA_maf_file('path/to/maf.tsv') # if TCGA

>>> # Step 2: Process further
>>> data_container.extend_df().to_counts_df('SBS_96', ed.categories.SBS_96_category_list())

>>> # Step 3: Access any processed dataframe of interest:
>>> ssm_df = data_container.ssm_df
>>> extended_df = data_container.extended_df
>>> counts_df = data_container.counts_dfs['SBS_96']


>>> # Alternatively, use without the chaining API:
>>> ssm_df = ed.standardize_ICGC_ssm_file('path/to/ssm.tsv', wrap=False) # if ICGC
>>> ssm_df = ed.standardize_TCGA_maf_file('path/to/maf.tsv', wrap=False) # if TCGA
>>> extended_df = ed.extend_ssm_df(ssm_df)
>>> counts_df = ed.counts_from_extended_ssm_df(
        extended_df, 
        category_colname='SBS_96', 
        category_values=ed.categories.SBS_96_category_list()
    )
```

With data already in the ExploSig "standard format":

```python
>>> import explosig_data as ed
>>> import pandas as pd

>>> # Step 0: Load the data into a dataframe, for example by reading from a TSV file.
>>> ssm_df = pd.read_csv('path/to/standard.tsv', sep='\t')

>>> # Step 1: Wrap the dataframe using the container class to allow use of the chainable functions.
>>> data_container = ed.SimpleSomaticMutationContainer(ssm_df)

>>> # Now see step 2 above (or the alternative steps above).
```


### Development

Build and install from the current directory.

```sh
python setup.py sdist bdist_wheel && pip install .
```