import os
import sys
import snakemake as snakemake_api
import tempfile
import yaml
import numpy as np
import pandas as pd
import logging

from .constants import *
from .i_o import get_logger, get_df_drop_message
from .categories import *

def _setup():
    print("Setting up explosig-data")

    os.makedirs(os.path.expanduser(os.path.join('~', '.explosig')), exist_ok=True)
    os.makedirs(os.path.expanduser(os.path.join('~', '.explosig', 'genomes')), exist_ok=True)

    
    # Download human genome files
    config = {
        "output": {
            "hg19": os.path.expanduser(os.path.join('~', '.explosig', "genomes", "hg19.fa")),
            "hg38": os.path.expanduser(os.path.join('~', '.explosig', "genomes", "hg38.fa"))
        }
    }

    # Since snakemake() function can only handle "flat" dicts using the direct config= parameter,
    # need to write the config dict to a temporary file and instead pass in to configfile=
    with tempfile.NamedTemporaryFile(mode='w') as temp:
        yaml.dump(config, temp, default_flow_style=False)
        snakefile = os.path.join(os.path.dirname(__file__), 'snakefiles', 'genomes', 'human.smk')
        snakemake_api.snakemake(snakefile=snakefile, configfiles=[temp.name])

_setup()

def convert_with_map(row, index, convert_map):
  try:
    return convert_map[row[index]]
  except KeyError:
    return NAN_VAL

def clean_ssm_df(df):
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
    df.sort_values([COLNAME.PATIENT.value, COLNAME.SAMPLE.value, COLNAME.CHR.value, COLNAME.POS_START.value], inplace=True)

    # Restrict to the standard set of columns
    return df[SSM_COLUMNS]


def convert_to_ssm_df_from_ICGC_file(input_file, filter_by_seq_type=None, 
                                        cancer_type='unknown', provenance='unknown', cohort='unknown', 
                                        console_verbosity=logging.DEBUG):
    get_logger(console_verbosity=console_verbosity)

    dtypes = {
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
    ssm_df = pd.read_csv(input_file, sep='\t', usecols=dtypes.keys(), dtype=dtypes)
    logging.debug("Input df has %d rows" % ssm_df.shape[0])

    # Standardize column names
    ssm_df.rename(columns={
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
    }, inplace=True)

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
    return ssm_df

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
def add_transcription_strand_column(df, transcripts):
    def determine_tstrand(row):
        try:
            tstrand = transcripts[row[COLNAME.ASSEMBLY.value]].strand(
                chr_name=row[COLNAME.CHR.value], 
                pos=row[COLNAME.POS_START.value], 
                gstrand=row[COLNAME.GSTRAND.value]
            )
        except ValueError:
            tstrand = NAN_VAL
        return tstrand

    df[COLNAME.TSTRAND.value] = df.apply(determine_tstrand, axis='columns')

    return df

def extend_ssm_df(ssm_df, category_functions=None, genomes=None, transcripts=None):
    if category_functions == None:
        category_functions = {
            'INDEL_Alexandrov2018_83': (INDEL_Alexandrov2018_83_category_name, [MUT_TYPE_VAL.INS.value, MUT_TYPE_VAL.DEL.value]),
            'DBS_78': (DBS_78_category_name, [MUT_TYPE_VAL.DBS.value]),
            'SBS_96': (SBS_96_category_name, [MUT_TYPE_VAL.SBS.value]),
        }
    
    if genomes == None:
        genomes = {
            ASSEMBLY_VAL.HG19.value: None,
            ASSEMBLY_VAL.HG38.value: None
        }
    
    if transcripts == None:
        transcripts = {
            ASSEMBLY_VAL.HG19.value: None,
            ASSEMBLY_VAL.HG38.value: None
        }
    
    ssm_df = add_flanking_columns(ssm_df, genomes)
    ssm_df = add_transcription_strand_column(ssm_df, transcripts)
    ssm_df = add_mutation_category_column(ssm_df, category_functions)

    #logging.info('Adding distance to previous mutation column')
    #ssm_df = add_dist_to_prev_mut_column(ssm_df)

    #logging.info('Adding rolling mean column')
    #ssm_df = add_rolling_mean_column(ssm_df)

    return ssm_df



    
    
    

