import os
import logging

from .constants import *
from .utils import clean_ssm_df
from .ssm_extended import extend_ssm_df
from .ssm_counts import counts_from_extended_ssm_df
from .data_source_ICGC import standardize_ICGC_ssm_file
from .data_source_TCGA import standardize_TCGA_maf_file


def _setup():
    os.makedirs(os.path.expanduser(os.path.join('~', '.explosig')), exist_ok=True)
    os.makedirs(os.path.expanduser(os.path.join('~', '.explosig', 'genes')), exist_ok=True)
    os.makedirs(os.path.expanduser(os.path.join('~', '.explosig', 'genomes')), exist_ok=True)

_setup()


    

