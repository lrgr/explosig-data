from urllib.request import urlretrieve
import gzip
import shutil
from math import floor

HG19_REFFLAT_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz'
HG38_REFFLAT_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz'

def print_download_progress(x, y, z):
    if floor(x*y*100/z) - floor((x-1)*y*100/z) >= 1:
        print("Download progress: {}%".format(floor(x*y*100/z)))

# Rules
rule genes_human_all:
  input:
    config['output']['hg19'],
    config['output']['hg38'],

rule genes_human_hg19_extract:
    input:
        config['output']['hg19'] + '.gz'
    output:
        config['output']['hg19']
    run:
        with gzip.open(config['output']['hg19'] + '.gz', 'rb') as f_in:
            with open(config['output']['hg19'], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

rule genes_human_hg38_extract:
    input:
        config['output']['hg38'] + '.gz'
    output:
        config['output']['hg38']
    run:
        with gzip.open(config['output']['hg38'] + '.gz', 'rb') as f_in:
            with open(config['output']['hg38'], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

rule genes_human_hg19_download:
    output:
        config['output']['hg19'] + '.gz'
    run:
        urlretrieve(
            HG19_REFFLAT_URL, 
            config['output']['hg19'] + '.gz', 
            reporthook=print_download_progress
        )

rule genes_human_hg38_download:
    output:
        config['output']['hg38'] + '.gz'
    run:
        urlretrieve(
            HG19_REFFLAT_URL, 
            config['output']['hg38'] + '.gz', 
            reporthook=print_download_progress
        )