import logging
import pandas as pd


from .constants import *
from .categories import *
from .i_o import get_logger, get_df_drop_message
from .genomes import get_human_genomes_dict
from .genes import get_human_genes_dict

# Add columns containing five prime and three prime flanking base pairs.
def add_flanking_columns(df, genomes):

    # Calculate number of flanking base pairs to add
    def get_flanking_size(row):
        return 6 * max( len(row[COLNAME.REF.value]), len(row[COLNAME.VAR.value]) )

    logging.info("Adding 5' flanking base column...")

    # if these are single-base substitution mutations, easy to pass in reference base during this step
    # to perform additional genome lookup assertions that can catch position indexing differences
    df[COLNAME.FPRIME.value] = df.apply(lambda row: genomes[row[COLNAME.ASSEMBLY.value]].lflank(
        chr_name=str(row[COLNAME.CHR.value]),
        pos=int(row[COLNAME.POS_START.value]),
        gstrand=row[COLNAME.GSTRAND.value],
        size=get_flanking_size(row),
        reference_base=(row[COLNAME.REF.value] if row[COLNAME.MUT_TYPE.value] == MUT_TYPE_VAL.SBS.value else None)), axis='columns')
    
    logging.info("Adding 3' flanking base column...")

    df[COLNAME.TPRIME.value] = df.apply(lambda row: genomes[row[COLNAME.ASSEMBLY.value]].rflank(
        chr_name=str(row[COLNAME.CHR.value]),
        pos=int(row[COLNAME.POS_END.value]),
        gstrand=row[COLNAME.GSTRAND.value],
        size=get_flanking_size(row),
        reference_base=(row[COLNAME.REF.value] if row[COLNAME.MUT_TYPE.value] == MUT_TYPE_VAL.SBS.value else None)), axis='columns')

    return df

# Add a category column for the given category name and list functions.
def add_mutation_category_column(df, category_functions):
    # Add category column
    for category_name, (category_name_func, mut_types) in category_functions.items():
        logging.info("Adding category {colname} column...".format(colname=category_name))
        df[category_name] = df.apply(lambda row: (category_name_func(row) if row[COLNAME.MUT_TYPE.value] in mut_types else NAN_VAL), axis='columns')

    return df

# Add a column specifying whether the mutation is on the transcribed or non-transcribed strand.
def add_transcription_strand_column(df, genes):
    def determine_tstrand(row):
        try:
            tstrand = genes[row[COLNAME.ASSEMBLY.value]].strand(
                chr_name=row[COLNAME.CHR.value], 
                pos=row[COLNAME.POS_START.value], 
                gstrand=row[COLNAME.GSTRAND.value]
            )
        except ValueError:
            tstrand = NAN_VAL
        return tstrand

    df[COLNAME.TSTRAND.value] = df.apply(determine_tstrand, axis='columns')

    return df

def extend_ssm_df(ssm_df, category_functions=None, genomes=None, genes=None, 
                    console_verbosity=logging.DEBUG):
    """Extend a standardized simple somatic mutation dataframe by adding the following columns: flanking bases, transcription strand, mutation category.
    
    Parameters
    ----------
    ssm_df : `pd.DataFrame`
        An already-standardized simple somatic mutation dataframe.
    category_functions : `dict`, optional
        Dictionary mapping category column names to tuples: (category_name_func, `list` of applicable mutation type enum values).
    genomes : `dict`, optional
        Dictionary mapping genome assembly enum values to Genome objects.
    genes : `dict`, optional
        Dictionary mapping genome assembly enum values to GeneLookup objects.
    
    Returns
    -------
    pd.DataFrame
        [description]
    """

    get_logger(console_verbosity=console_verbosity)

    if category_functions == None:
        category_functions = {
            'INDEL_Alexandrov2018_83': (INDEL_Alexandrov2018_83_category_name, [MUT_TYPE_VAL.INS.value, MUT_TYPE_VAL.DEL.value]),
            'DBS_78': (DBS_78_category_name, [MUT_TYPE_VAL.DBS.value]),
            'SBS_96': (SBS_96_category_name, [MUT_TYPE_VAL.SBS.value]),
        }
    
    if genomes == None:
        genomes = get_human_genomes_dict()
    
    if genes == None:
        genes = get_human_genes_dict()
    
    ssm_df = add_flanking_columns(ssm_df, genomes)
    ssm_df = add_transcription_strand_column(ssm_df, genes)
    ssm_df = add_mutation_category_column(ssm_df, category_functions)

    #logging.info('Adding distance to previous mutation column')
    #ssm_df = add_dist_to_prev_mut_column(ssm_df)

    #logging.info('Adding rolling mean column')
    #ssm_df = add_rolling_mean_column(ssm_df)

    return ssm_df



    
    