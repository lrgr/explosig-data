import pandas as pd
import logging

from .constants import *
from .i_o import get_logger, get_df_drop_message


# Helper functions
def convert_with_map(row, index, convert_map):
  try:
    return convert_map[row[index]]
  except KeyError:
    return NAN_VAL

def clean_ssm_df(df):
    """Perform the final stage of standardization of a simple somatic mutation dataframe.
    
    Parameters
    ----------
    df : `pd.DataFrame`
        A simple somatic mutation dataframe that contains all of the expected columns.
    
    Returns
    -------
    `pd.DataFrame`
        The dataframe with typed columns, sorted rows, and filtered rows (filtered if NaN/invalid chromosome, NaN start pos, or NaN end pos).
    """
    # Drop mutations with NaN chromosome
    filtered_df = df.dropna(subset=[COLNAME.CHR.value])
    logging.debug(get_df_drop_message(COLNAME.CHR.value, "NaN value", df, filtered_df))
    df = filtered_df

    # Drop mutations with NaN start position
    filtered_df = df.dropna(subset=[COLNAME.POS_START.value])
    logging.debug(get_df_drop_message(COLNAME.POS_START.value, "NaN value", df, filtered_df))
    df = filtered_df

    # Drop mutations with NaN end position
    filtered_df = df.dropna(subset=[COLNAME.POS_END.value])
    logging.debug(get_df_drop_message(COLNAME.POS_END.value, "NaN value", df, filtered_df))
    df = filtered_df
    
    # Drop mutations with invalid chromosome
    filtered_df = df.loc[df[COLNAME.CHR.value].isin(CHROMOSOMES)]
    logging.debug(get_df_drop_message(COLNAME.CHR.value, "invalid value", df, filtered_df))
    df = filtered_df

    # Ensure correct types before sorting
    df[COLNAME.CHR.value] = df[COLNAME.CHR.value].apply(str) # make sure everything is a string
    df[COLNAME.POS_START.value] = df[COLNAME.POS_START.value].astype(int)
    df[COLNAME.POS_END.value] = df[COLNAME.POS_END.value].astype(int)

    # Sort the mutations by sample and then genomic location
    df[COLNAME.CHR.value] = pd.Categorical(df[COLNAME.CHR.value], CHROMOSOMES, ordered=True)
    df = df.sort_values([COLNAME.PATIENT.value, COLNAME.SAMPLE.value, COLNAME.CHR.value, COLNAME.POS_START.value])

    # Restrict to the standard set of columns
    return df[SSM_COLUMNS]


