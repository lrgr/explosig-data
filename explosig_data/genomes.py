import os
import logging
from abc import abstractmethod
from Bio import SeqIO
from Bio.Seq import reverse_complement
import twobitreader

import snakemake as snakemake_api
import tempfile
import yaml

from .constants import *

class Genome:
    @abstractmethod
    def __init__(self, genome_filepath):
        raise NotImplementedError
    
    @abstractmethod
    def seq(self, chr_name, start, end, gstrand):
        raise NotImplementedError
    
    @abstractmethod
    def base(self, chr_name, pos, gstrand):
        raise NotImplementedError
    
    def lflank(self, chr_name, pos, gstrand, size, reference_base=None):
        if reference_base != None:
            # optional verification (convenient if computing for SBS mutation)
            reference_base_from_genome = self.base(chr_name=chr_name, pos=pos, gstrand=gstrand)
            assert(reference_base_from_genome == 'N' or reference_base_from_genome == reference_base)
        return self.seq(chr_name=chr_name, start=pos-size-1, end=pos-1, gstrand=gstrand)
    
    def rflank(self, chr_name, pos, gstrand, size, reference_base=None):
        if reference_base != None:
            # optional verification (convenient if computing for SBS mutation)
            reference_base_from_genome = self.base(chr_name=chr_name, pos=pos, gstrand=gstrand)
            assert(reference_base_from_genome == 'N' or reference_base_from_genome == reference_base)
        return self.seq(chr_name=chr_name, start=pos, end=pos+size, gstrand=gstrand)

class FastaGenome(Genome):
    def __init__(self, genome_filepath):
        logging.debug('Loading genome...')

        with open(genome_filepath, "r") as IN:
            self.genome = SeqIO.to_dict(SeqIO.parse(IN, "fasta"))
        logging.debug('Loading genome complete')
        
    def seq(self, chr_name, start, end, gstrand):
        assert (gstrand == GSTRAND_VAL.PLUS.value) # TODO update position when GSTRAND is not plus
        return str(self.genome[chr_name][start:end].seq)

    
    def base(self, chr_name, pos, gstrand):
        assert (gstrand == GSTRAND_VAL.PLUS.value) # TODO update position when GSTRAND is not plus
        return str(self.genome[chr_name][pos-1])


class TwoBitGenome(Genome):
    def __init__(self, genome_filepath):
        logging.debug('Loading genome...')
        
        self.genome = twobitreader.TwoBitFile(genome_filepath)
        logging.info('Loading genome complete')
        
    def seq(self, chr_name, start, end, gstrand):
        assert (gstrand == GSTRAND_VAL.PLUS.value) # TODO update position when GSTRAND is not plus
        return self.genome[chr_name][start:end].upper()
    
    def base(self, chr_name, pos, gstrand):
        assert (gstrand == GSTRAND_VAL.PLUS.value) # TODO update position when GSTRAND is not plus
        return str(self.genome[chr_name][pos-1]).upper()

def download_human_genomes():
    config = {
        "output": {
            "hg19": os.path.join(EXPLOSIG_DATA_DIR, "genomes", "hg19.fa"),
            "hg38": os.path.join(EXPLOSIG_DATA_DIR, "genomes", "hg38.fa")
        }
    }

    # Since snakemake() function can only handle "flat" dicts using the direct config= parameter,
    # need to write the config dict to a temporary file and instead pass in to configfile=
    with tempfile.NamedTemporaryFile(mode='w') as temp:
        yaml.dump(config, temp, default_flow_style=False)
        snakefile = os.path.join(os.path.dirname(__file__), 'snakefiles', 'genomes', 'human.smk')
        snakemake_api.snakemake(snakefile=snakefile, configfiles=[temp.name])

def get_human_genomes_dict():
    download_human_genomes()
    return {
        ASSEMBLY_VAL.HG19.value: FastaGenome(os.path.join(EXPLOSIG_DATA_DIR, "genomes", "hg19.fa")),
        ASSEMBLY_VAL.HG38.value: FastaGenome(os.path.join(EXPLOSIG_DATA_DIR, "genomes", "hg38.fa"))
    }