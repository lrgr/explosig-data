import os
from enum import Enum

EXPLOSIG_DATA_DIR = os.path.expanduser(os.path.join('~', '.explosig'))

BASES       = [ 'A', 'C', 'G', 'T' ]
CHROMOSOMES = list(map(str, list(range(1, 23)))) + ['X', 'Y']
PURINES     = set('AG')
BASE_PAIR   = dict(A='T', C='G', G='C', T='A')

NAN_VAL     = 'NaN'

class CONSORTIA(Enum):
    ICGC = "ICGC"

class COLNAME(Enum):
    # Column names for mutation tables
    PATIENT = 'Patient'
    SAMPLE = 'Sample'
    ALT_SAMPLE = 'Alternative Sample Name'
    CANCER_TYPE = 'Cancer Type'
    PROVENANCE = 'Provenance'
    COHORT = 'Cohort'
    CHR = 'Chromosome'
    POS_START = 'Start Position'
    POS_END = 'End Position'
    REF = 'Reference Sequence'
    VAR = 'Variant Sequence'
    GSTRAND = 'Genomic Strand' # Watson-Crick Strand
    SEQ_TYPE = 'Sequencing Strategy'
    MUT_TYPE = 'Mutation Type'
    MUT_CLASS = 'Mutation Classification'
    GENE_SYMBOL = 'Gene Symbol'
    ASSEMBLY = 'Assembly Version'

    # Supplemental column names
    FPRIME = "5' Flanking Bases"
    TPRIME = "3' Flanking Bases"
    TSTRAND = 'Transcriptional Strand'
    MUT_DIST = 'Distance to Previous Mutation'
    NEAREST_MUT = 'Distance to Nearest Mutation'
    MUT_DIST_ROLLING_MEAN = 'Rolling Mean of 6 Mutation Distances'


SSM_COLUMNS = [
    COLNAME.PATIENT.value, 
    COLNAME.SAMPLE.value, 
    COLNAME.CANCER_TYPE.value, 
    COLNAME.PROVENANCE.value, 
    COLNAME.COHORT.value,
    COLNAME.CHR.value, 
    COLNAME.POS_START.value, 
    COLNAME.POS_END.value, 
    COLNAME.REF.value, 
    COLNAME.VAR.value, 
    COLNAME.GSTRAND.value, 
    COLNAME.SEQ_TYPE.value, 
    COLNAME.MUT_TYPE.value, 
    COLNAME.ASSEMBLY.value
]

class TSTRAND_VAL(Enum):
    PLUS = '+'
    MINUS = '-'

class GSTRAND_VAL(Enum):
    PLUS = '+'
    MINUS = '-'

class SEQ_TYPE_VAL(Enum):
    WXS = 'WXS'
    WGS = 'WGS'
    RNASEQ = 'RNA-Seq'
    MRXS = 'MRE-Seq' # multi region exome sequencing

class MUT_TYPE_VAL(Enum):
    SBS = 'SBS'
    INS = 'INS'
    DEL = 'DEL'
    DBS = 'DBS'
    TBS = 'TBS'

class ASSEMBLY_VAL(Enum):
    HG19 = 'GRCh37'
    HG38 = 'GRCh38'