import logging
import pandas as pd

from .constants import *
from .utils import clean_ssm_df, convert_with_map
from .i_o import get_logger, get_df_drop_message
from .ssm_container import SimpleSomaticMutationContainer

col_dtypes = {
    "Tumor_Sample_Barcode": str,
    "Matched_Norm_Sample_Barcode": str,
    "Chromosome": str,
    "Start_Position": object,
    "End_Position": object,
    "Variant_Type": str,
    "Strand": str,
    "STRAND": str,
    "Reference_Allele": str,
    "Tumor_Seq_Allele1": str,
    "Tumor_Seq_Allele2": str,
    "Match_Norm_Seq_Allele1": str,
    "Match_Norm_Seq_Allele2": str,
    "FILTER": str,
    "CONTEXT": str,
    "NCBI_Build": str,
    "Hugo_Symbol": str,
    "Variant_Classification": str
}

col_renames = {
    "Tumor_Sample_Barcode": COLNAME.SAMPLE.value,
    "Chromosome": COLNAME.CHR.value,
    "Start_Position": COLNAME.POS_START.value,
    "End_Position": COLNAME.POS_END.value,
    "Variant_Type": COLNAME.MUT_TYPE.value,
    "STRAND": COLNAME.TSTRAND.value,
    "Strand": COLNAME.GSTRAND.value,
    "Reference_Allele": COLNAME.REF.value,
    "Tumor_Seq_Allele2": COLNAME.VAR.value,
    "NCBI_Build": COLNAME.ASSEMBLY.value,
    "Variant_Classification": COLNAME.MUT_CLASS.value,
    "Hugo_Symbol": COLNAME.GENE_SYMBOL.value
}

def standardize_TCGA_maf_file(input_maf_file, wrap=True,
                                        cancer_type='unknown', provenance='unknown', cohort='unknown', 
                                        col_dtypes=col_dtypes, col_renames=col_renames,
                                        console_verbosity=logging.DEBUG):
    """Convert to explosig simple somatic mutation ("standard") format from the TCGA PanCanAtlas MAF format.
    
    Parameters
    ----------
    input_maf_file : `str`
        Path to a TCGA PanCanAtlas MAF file.
    wrap : `bool`, optional
        Whether to wrap the return value for chaining, by default `True`
    cancer_type : `str`, optional
        The cancer type value on which to filter.
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

    maf_df = pd.read_csv(input_maf_file, sep="\t", usecols=col_dtypes.keys(), dtype=col_dtypes)
    logging.debug("Input df has %d rows" % maf_df.shape[0])
    
    maf_df = maf_df.rename(columns=col_renames)
    # set sequencing strategy to be whole exome sequencing (WXS)
    # not part of the MAF but WR manually verified via the mc3 paper (Ellrot et al 2018)
    maf_df[COLNAME.SEQ_TYPE.value] = SEQ_TYPE_VAL.WXS.value
    # set patient to be first 12 characters of sample
    maf_df[COLNAME.PATIENT.value] = [barcode[0:12] for barcode in maf_df[COLNAME.SAMPLE.value]]
    
    # remove mutations where Filter column contains 'nonpreferredpair' or 'oxog' or 'StrandBias'
    maf_df = maf_df.loc[~maf_df["FILTER"].str.contains('StrandBias|oxog|nonpreferredpair')]

    logging.debug("After removing mutations where FILTER column contains 'nonpreferredpair' or 'oxog' or 'StrandBias', df has %d rows" % maf_df.shape[0])

    # set cohort and provenance
    maf_df[COLNAME.COHORT.value] = cohort
    maf_df[COLNAME.PROVENANCE.value] = provenance
    maf_df[COLNAME.CANCER_TYPE.value] = cancer_type

    # TODO: update this indel logic
    def convert_mut_type(row):
        # required because len(NaN) returns a TypeError "object of type 'float' has no len()"
        if row.isnull()[[COLNAME.REF.value, COLNAME.VAR.value]].any():
            return NAN_VAL
        if row[COLNAME.MUT_TYPE.value] == 'SNP':
            return MUT_TYPE_VAL.SBS.value
        elif len(row[COLNAME.REF.value]) == 2 and len(row[COLNAME.VAR.value]) == 2:
            return MUT_TYPE_VAL.DBS.value
        elif len(row[COLNAME.REF.value]) == 1 and row[COLNAME.REF.value] == '-':
            return MUT_TYPE_VAL.INS.value
        elif len(row[COLNAME.VAR.value]) == 1 and row[COLNAME.VAR.value] == '-':
            return MUT_TYPE_VAL.DEL.value
        else:
            return NAN_VAL
    
    maf_df[COLNAME.MUT_TYPE.value] = maf_df.apply(convert_mut_type, axis=1)

    logging.debug("Assigned mutation types resulting in %d SBS, %d DBS, %d INS, %d DEL, %d NaN" % (
        maf_df.loc[maf_df[COLNAME.MUT_TYPE.value] == MUT_TYPE_VAL.SBS.value].shape[0],
        maf_df.loc[maf_df[COLNAME.MUT_TYPE.value] == MUT_TYPE_VAL.DBS.value].shape[0],
        maf_df.loc[maf_df[COLNAME.MUT_TYPE.value] == MUT_TYPE_VAL.INS.value].shape[0],
        maf_df.loc[maf_df[COLNAME.MUT_TYPE.value] == MUT_TYPE_VAL.DEL.value].shape[0],
        maf_df.loc[maf_df[COLNAME.MUT_TYPE.value] == NAN_VAL].shape[0]
    ))

    maf_df = clean_ssm_df(maf_df)
    
    if wrap:
        return SimpleSomaticMutationContainer(maf_df)
    else:
        return maf_df