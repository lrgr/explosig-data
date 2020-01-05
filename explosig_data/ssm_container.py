import pandas as pd
import logging

from .ssm_extended import extend_ssm_df
from .ssm_counts import counts_from_extended_ssm_df

class SimpleSomaticMutationContainer(object):

    def __init__(self, ssm_df):
        self.ssm_df = ssm_df
        self.extended_df = None
        self.counts_dfs = {}
    
    def extend_df(self, **kwargs):
        self.extended_df = extend_ssm_df(self.ssm_df, **kwargs)
        return self
    
    def to_counts_df(self, category_colname, category_values, **kwargs):
        self.counts_dfs[category_colname] = counts_from_extended_ssm_df(self.extended_df, category_colname, category_values, **kwargs)
        return self