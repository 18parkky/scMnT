import os
import time
import argparse
import pandas as pd

start_time = time.time()

script_description = 'Given allele table, split it into smaller allele tables by flag and save to disk'

parser = argparse.ArgumentParser(description=script_description)
# Required parameters
parser.add_argument('-at', '--dir_allele_table',   help="Directory of input MS allele table", required=True)
# Optional parameters
parser.add_argument('-PATH', '--PATH_out',          help='PATH of the output files (default: current directory)', required=False, type=str, default=os.getcwd() )

args = vars(parser.parse_args())

dir_MS_allele_table = args["dir_allele_table"]
PATH_out            = args["PATH_out"]

df = pd.read_csv(dir_MS_allele_table, sep='\t')

filename    = os.path.splitext( os.path.basename( dir_MS_allele_table ))[0]

for flag, edf in df.groupby("flag"):
    edf.to_csv(f"{PATH_out}/{filename}.flag_{flag}.tsv", sep='\t', index=None) 

print(f"Finished split_MSAT_byFlag.py\t(Total time taken: {round(time.time() - start_time, 2)} seconds)")