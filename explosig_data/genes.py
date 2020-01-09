import os
import bisect
import pandas as pd

import snakemake as snakemake_api
import tempfile
import yaml

from .constants import *

def row_matches(row, pos):
    if row[COLNAME.POS_START.value] <= pos and row[COLNAME.POS_END.value] >= pos:
        return True
    return False

def search_left(i, pos, chr_pos_list, chr_df):
    found_rows = []
    while i >= 0 and i < len(chr_pos_list):
        potential_row = chr_df.iloc[i]
        if row_matches(potential_row, pos):
            found_rows.append(potential_row)
        i -= 1
        if potential_row[COLNAME.POS_END.value] < pos:
            break
    return found_rows

def search_right(i, pos, chr_pos_list, chr_df):
    found_rows = []
    while i >= 0 and i < len(chr_pos_list):
        potential_row = chr_df.iloc[i]
        if row_matches(potential_row, pos):
            found_rows.append(potential_row)
        i += 1
        if potential_row[COLNAME.POS_START.value] > pos:
            break
    return found_rows
            
class GeneLookup:
    def __init__(self, transcripts_filepath):
        colnames = ['gene', COLNAME.CHR.value, COLNAME.TSTRAND.value, COLNAME.POS_START.value, COLNAME.POS_END.value]
        df = pd.read_csv(transcripts_filepath, sep='\t', usecols=[0, 2, 3, 4, 5], dtype={0:str, 2:str, 3:str, 4:int, 5:int}, header=None, names=colnames)
        # Drop the UCSC-style 'chr' prefix
        df[COLNAME.CHR.value] = df.apply(lambda row: row[COLNAME.CHR.value][3:], axis='columns')
        # Drop the non-standard chromosomes
        df = df.loc[df[COLNAME.CHR.value].isin(CHROMOSOMES)]
        # Sort so we can do binary search
        df = df.sort_values(by=[COLNAME.CHR.value, COLNAME.POS_END.value, COLNAME.POS_START.value])
        # Filter the dataframe so it can be accessed by assembly and chromosome
        self.df = {}
        self.positions = {}
        for chromosome in CHROMOSOMES:
            self.df[chromosome] = df.loc[df[COLNAME.CHR.value] == chromosome]
            # Store the position column as a list so that python's bisect can be used
            self.positions[chromosome] = self.df[chromosome][COLNAME.POS_END.value].tolist()
    
    def strand(self, chr_name, pos, gstrand):
        chr_name = str(chr_name)
        pos = int(pos)
        assert (chr_name in CHROMOSOMES)
        assert (gstrand == GSTRAND_VAL.PLUS.value) # TODO update position when GSTRAND is not plus

        # Restrict df to current chromosome
        chr_df = self.df[chr_name]
        chr_pos_list = self.positions[chr_name]
        
        found_rows = []
        # Search
        i = bisect.bisect_left(chr_pos_list, pos)
        if i and i >= 0:
            found_rows = (
                search_left(i, pos, chr_pos_list, chr_df) + 
                search_right(i, pos, chr_pos_list, chr_df)
            )
        else:
            raise ValueError("No transcript matches found")
        
        if len(found_rows) > 0:
            # Check for conflicts
            strands = set()
            for found_row in found_rows:
                strands.add(found_row[COLNAME.TSTRAND.value])
            if strands == {TSTRAND_VAL.PLUS.value, TSTRAND_VAL.MINUS.value}:
                return ('%s,%s' % (TSTRAND_VAL.PLUS.value, TSTRAND_VAL.MINUS.value))
            elif len(strands) == 1:
                return list(strands)[0]
            else:
                raise ValueError("Not sure what is going on")
        else:
            raise ValueError("No transcript matches found.")

def download_human_genes():
    config = {
        "output": {
            "hg19": os.path.join(EXPLOSIG_DATA_DIR, "genes", "refFlat19.txt"),
            "hg38": os.path.join(EXPLOSIG_DATA_DIR, "genes", "refFlat38.txt")
        }
    }

    # Since snakemake() function can only handle "flat" dicts using the direct config= parameter,
    # need to write the config dict to a temporary file and instead pass in to configfile=
    with tempfile.NamedTemporaryFile(mode='w') as temp:
        yaml.dump(config, temp, default_flow_style=False)
        snakefile = os.path.join(os.path.dirname(__file__), 'snakefiles', 'genes', 'human.smk')
        snakemake_api.snakemake(snakefile=snakefile, configfiles=[temp.name])

def get_human_genes_dict():
    download_human_genes()
    return {
        ASSEMBLY_VAL.HG19.value: GeneLookup(os.path.join(EXPLOSIG_DATA_DIR, "genes", "refFlat19.txt")),
        ASSEMBLY_VAL.HG38.value: GeneLookup(os.path.join(EXPLOSIG_DATA_DIR, "genes", "refFlat38.txt"))
    }