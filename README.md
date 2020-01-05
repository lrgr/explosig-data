[![Build Status](https://travis-ci.org/lrgr/explosig-data.svg?branch=master)](https://travis-ci.org/lrgr/explosig-data)
[![PyPI](https://img.shields.io/pypi/v/explosig-data)](https://pypi.org/project/explosig-data/)

## ExploSig Data

Helpers for processing mutation data into standard formats originally developed for the [ExploSig](https://github.com/lrgr/explosig) family of tools.

- [Documentation](https://lrgr.github.io/explosig-data/)

### Installation

```sh
pip install explosig-data
```

### Example Usage

```python
>>> import explosig_data as ed

>>> # With chaining
>>> container = (ed
        .standardize_ICGC_ssm_file('path/to/ssm.tsv')
        .extend_df()
        .to_counts_df('SBS_96', ed.categories.SBS_96_category_list())
    )
>>> counts_df = container.counts_dfs['SBS_96']


>>> # Without chaining
>>> ssm_df = ed.standardize_ICGC_ssm_file('path/to/ssm.tsv', wrap=False)
>>> # or
>>> ssm_df = ed.standardize_TCGA_maf_file('path/to/maf.tsv', wrap=False)
>>> extended_df = ed.extend_ssm_df(ssm_df)
>>> counts_df = ed.counts_from_extended_ssm_df(
        extended_df, 
        category_colname='SBS_96', 
        category_values=ed.categories.SBS_96_category_list()
    )
```


### Development

Build and install from the current directory.

```sh
python setup.py sdist bdist_wheel && pip install .
```