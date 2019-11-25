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
>>> ed.extend_ssm_df(ed.convert_to_ssm_df_from_ICGC_file('~/Desktop/simple_somatic_mutation.open.ALL-US.tsv.gz'))
```


### Development

Build and install from the current directory.

```sh
python setup.py sdist bdist_wheel && pip install .
```