
import logging
import pandas as pd

from .constants import *
from .utils import clean_ssm_df, convert_with_map
from .i_o import get_logger, get_df_drop_message
from .ssm_container import SimpleSomaticMutationContainer

col_dtypes = {
    'icgc_mutation_id': str,
    'icgc_donor_id': str,
    'icgc_sample_id': str,
    'chromosome': str,
    'chromosome_start': object,
    'chromosome_end': object,
    'chromosome_strand': str,
    'assembly_version': str,
    'mutation_type': str,
    'reference_genome_allele': str,
    'mutated_to_allele': str,
    'total_read_count': object,
    'mutant_allele_read_count': object,
    'sequencing_strategy': str
}
col_renames = {
    "icgc_donor_id": COLNAME.PATIENT.value,
    "icgc_sample_id": COLNAME.SAMPLE.value,
    "chromosome": COLNAME.CHR.value,
    "chromosome_start": COLNAME.POS_START.value,
    "chromosome_end": COLNAME.POS_END.value,
    "reference_genome_allele": COLNAME.REF.value,
    "mutated_to_allele": COLNAME.VAR.value,
    "mutation_type": COLNAME.MUT_TYPE.value,
    "sequencing_strategy": COLNAME.SEQ_TYPE.value,
    "assembly_version": COLNAME.ASSEMBLY.value,
    "chromosome_strand": COLNAME.GSTRAND.value
}

def standardize_ICGC_ssm_file(input_ssm_file, wrap=True, filter_by_seq_type=None, 
                                        cancer_type='unknown', provenance='unknown', cohort='unknown', 
                                        col_dtypes=col_dtypes, col_renames=col_renames,
                                        console_verbosity=logging.DEBUG):
    """Convert to explosig simple somatic mutation ("standard") format from the ICGC simple somatic mutation format.
    
    Parameters
    ----------
    input_ssm_file : `str`
        Path to the ICGC simple somatic mutation file.
    wrap : `bool`, optional
        Whether to wrap the return value for chaining, by default `True`
    filter_by_seq_type : `str` or `list`, optional
        A sequencing type or list of sequencing types by which to filter, by default None
    cancer_type : `str`, optional
        Value to fill the Cancer Type column, by default 'unknown'
    provenance : `str`, optional
        Value to fill the Provenance column, by default 'unknown'
    cohort : `str`, optional
        Value to fill the Cohort column, by default 'unknown'
    col_dtypes : `dict`, optional
        Dictionary mapping input column names to data types.
    col_renames : `dict`, optional
        Dictionary mapping input column names to standard column name constants.
    console_verbosity : `int`, optional
        Logging verbosity, by default `logging.DEBUG`
    
    Returns
    -------
    `pd.DataFrame`
        The simple somatic mutation dataframe in a standardized format. This dataframe can be passed to the `extend...` functions.
    """
    get_logger(console_verbosity=console_verbosity)

    
    ssm_df = pd.read_csv(input_ssm_file, sep='\t', usecols=col_dtypes.keys(), dtype=col_dtypes)
    logging.debug("Input df has %d rows" % ssm_df.shape[0])

    # Standardize column names
    ssm_df = ssm_df.rename(columns=col_renames)

    # Mapped sequencing types
    seq_type_map = {
        'RNA-Seq': SEQ_TYPE_VAL.RNASEQ.value,
        'WXS': SEQ_TYPE_VAL.WXS.value,
        'WGS': SEQ_TYPE_VAL.WGS.value
    }
    ssm_df[COLNAME.SEQ_TYPE.value] = ssm_df.apply(lambda row: convert_with_map(row, COLNAME.SEQ_TYPE.value, seq_type_map), axis='columns')

    logging.debug("Standardized sequencing types resulting in %d WGS, %d WXS, %d RNA-Seq, %d NaN rows" % (
        ssm_df.loc[ssm_df[COLNAME.SEQ_TYPE.value] == SEQ_TYPE_VAL.WGS.value].shape[0],
        ssm_df.loc[ssm_df[COLNAME.SEQ_TYPE.value] == SEQ_TYPE_VAL.WXS.value].shape[0],
        ssm_df.loc[ssm_df[COLNAME.SEQ_TYPE.value] == SEQ_TYPE_VAL.RNASEQ.value].shape[0],
        ssm_df.loc[ssm_df[COLNAME.SEQ_TYPE.value] == NAN_VAL].shape[0]
    ))

    if filter_by_seq_type != None:
        if type(filter_by_seq_type) == str:
            ssm_df = ssm_df.loc[ssm_df[COLNAME.SEQ_TYPE.value] == filter_by_seq_type]
        elif type(filter_by_seq_type) == list:
            ssm_df = ssm_df.loc[ssm_df[COLNAME.SEQ_TYPE.value].isin(filter_by_seq_type)]
        logging.debug("After restricting to sequencing type %s, df has %d rows" % (str(filter_by_seq_type), ssm_df.shape[0]))

    ssm_df.sort_values(by=[COLNAME.PATIENT.value, COLNAME.SAMPLE.value, COLNAME.POS_START.value, "total_read_count"], na_position='last', inplace=True)
    
    # In ICGC ssm files, identical mutations often have multiple rows because there is a different row for each gene consequence.
    # May also have multiple rows for the same mutation if the sample had both WXS and WGS sequencing, for example.
    ssm_df.drop_duplicates(subset=["icgc_mutation_id", COLNAME.PATIENT.value, COLNAME.SAMPLE.value, COLNAME.SEQ_TYPE.value], keep='first', inplace=True)

    logging.debug("After dropping rows with duplicate mutation ID, patient ID, sample ID, and sequencing type, df has %d rows" % ssm_df.shape[0])
    
    ssm_df[COLNAME.CANCER_TYPE.value], ssm_df[COLNAME.PROVENANCE.value], ssm_df[COLNAME.COHORT.value] = cancer_type, provenance, cohort

    # TODO: update this indel logic
    def convert_mut_type(row):
        if row[COLNAME.MUT_TYPE.value] == 'single base substitution':
            return MUT_TYPE_VAL.SBS.value
        elif len(row[COLNAME.REF.value]) == 2 and len(row[COLNAME.VAR.value]) == 2:
            return MUT_TYPE_VAL.DBS.value
        elif len(row[COLNAME.REF.value]) == 1 and row[COLNAME.REF.value] == '-':
            return MUT_TYPE_VAL.INS.value
        elif len(row[COLNAME.VAR.value]) == 1 and row[COLNAME.VAR.value] == '-':
            return MUT_TYPE_VAL.DEL.value
        else:
            return NAN_VAL
    
    ssm_df[COLNAME.MUT_TYPE.value] = ssm_df.apply(convert_mut_type, axis='columns')

    logging.debug("Assigned mutation types resulting in %d SBS, %d DBS, %d INS, %d DEL, %d NaN" % (
        ssm_df.loc[ssm_df[COLNAME.MUT_TYPE.value] == MUT_TYPE_VAL.SBS.value].shape[0],
        ssm_df.loc[ssm_df[COLNAME.MUT_TYPE.value] == MUT_TYPE_VAL.DBS.value].shape[0],
        ssm_df.loc[ssm_df[COLNAME.MUT_TYPE.value] == MUT_TYPE_VAL.INS.value].shape[0],
        ssm_df.loc[ssm_df[COLNAME.MUT_TYPE.value] == MUT_TYPE_VAL.DEL.value].shape[0],
        ssm_df.loc[ssm_df[COLNAME.MUT_TYPE.value] == NAN_VAL].shape[0]
    ))

    assembly_map = {
        'GRCh37': ASSEMBLY_VAL.HG19.value,
        'GRCh38': ASSEMBLY_VAL.HG38.value
    }
    ssm_df[COLNAME.ASSEMBLY.value] = ssm_df.apply(lambda row: convert_with_map(row, COLNAME.ASSEMBLY.value, assembly_map), axis='columns')

    gstrand_map = {
        '1': GSTRAND_VAL.PLUS.value
    }
    ssm_df[COLNAME.GSTRAND.value] = ssm_df.apply(lambda row: convert_with_map(row, COLNAME.GSTRAND.value, gstrand_map), axis='columns')

    ssm_df = clean_ssm_df(ssm_df)
    
    if wrap:
        return SimpleSomaticMutationContainer(ssm_df)
    else:
        return ssm_df
