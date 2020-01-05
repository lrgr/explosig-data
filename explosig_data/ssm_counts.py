import logging
import pandas as pd

from .constants import *
from .categories import *
from .i_o import get_logger, get_df_drop_message


def counts_from_extended_ssm_df(extended_df, category_colname, category_values,
                                sparse_output=False, console_verbosity=logging.DEBUG):
    """Construct a count matrix dataframe from a simple somatic mutation dataframe that has already been "extended".
    
    Parameters
    ----------
    extended_df : `pd.DataFrame`
        An extended simple somatic mutation dataframe (e.g. produced by the `extend_ssm_df` function).
    category_colname : `str`
        The category column name.
    category_values : `list`
        A list of all possible values for the category column. These will become the column names of the output dataframe.
    sparse_output : `bool`, optional
        Whether the returned dataframe will be in a sparse format, by default `False`
    console_verbosity : `int`, optional
        Logging verbosity enum value, by default `logging.DEBUG`
    
    Returns
    -------
    `pd.DataFrame`
        A mutation count dataframe. If not sparse, index is sample IDs, columns are category values, cells are count values.
    
    Raises
    ------
    `ValueError`
        Raises error if expected columns are missing from the input dataframe.
    """
    
    ssm_df = extended_df
    categories = category_values

    dtypes = {
        COLNAME.PATIENT.value: str,
        COLNAME.SAMPLE.value: str,
        COLNAME.CHR.value: str,
        COLNAME.POS_START.value: object,
        COLNAME.POS_END.value: object,
        COLNAME.TSTRAND.value: str,
        COLNAME.REF.value: str,
        COLNAME.VAR.value: str,
        COLNAME.FPRIME.value: str,
        COLNAME.TPRIME.value: str,
        COLNAME.MUT_TYPE.value: str,
        category_colname: str
    }

    # Check that the input df contains the expected columns.
    missing_cols = list(set(dtypes.keys()) - set(ssm_df.columns.values))
    if len(missing_cols) == 1 and missing_cols[0] == category_colname:
        raise ValueError("Input dataframe is missing the category column.")
    elif len(missing_cols) > 1:
        raise ValueError("Input dataframe is missing too many columns.")

    # Filter out mutations with categories not in our lists
    filtered_df = ssm_df.loc[ssm_df[category_colname].isin(set(categories))].copy()
    logging.debug(get_df_drop_message(category_colname, "invalid value", ssm_df, filtered_df))
    ssm_df = filtered_df

    # remove na's
    filtered_df = ssm_df.dropna(subset=[COLNAME.VAR.value])
    logging.debug(get_df_drop_message(COLNAME.VAR.value, "NaN value", ssm_df, filtered_df))
    ssm_df = filtered_df

    filtered_df = ssm_df.dropna(subset=[COLNAME.REF.value])
    logging.debug(get_df_drop_message(COLNAME.REF.value, "NaN value", ssm_df, filtered_df))
    ssm_df = filtered_df

    groups = ssm_df.groupby([COLNAME.SAMPLE.value, category_colname])

    counts_df = groups.size().reset_index(name='counts')

    # TODO: figure out how to factor in transcription strand column. and donor column.
    #       (easy with sparse output format (just add to the groupby),
    #           but for matrix-style output format need to look into using pandas MultiIndex columns)

    if sparse_output:
        return counts_df
    else:
        # Initialize sample x context matrix to hold mutation context counts
        counts_matrix_df = pd.DataFrame(data=0, columns=categories, index=list(ssm_df[COLNAME.SAMPLE.value].unique()))
        # Pivot to get a (possibly) incomplete matrix
        partial_matrix_df = counts_df.pivot(index=COLNAME.SAMPLE.value, columns=category_colname, values='counts')
        # Expand back out and fill empty cells with zeros
        counts_matrix_df = counts_matrix_df.add(partial_matrix_df, fill_value=0)
        # Select columns using the category list array (pandas add function does reorder the columns for some reason)
        counts_matrix_df = counts_matrix_df[categories]
        return counts_matrix_df