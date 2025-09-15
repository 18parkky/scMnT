import os
import time
import logging
import argparse 

import scipy.stats
import numpy as np
import pandas as pd
import scanpy as sc

import NanoMnT.nanomnt_utility as nanomnt_utility

def preprocessAlleleTable( AlleleTable, min_STR_length=10, max_STR_length=24, include_GC_repeats=False ):
    
    ### 1. Filter out low-quality flankings (e.g., indels within flankings)        
    col_flanking_quality = list()
    for tup2 in AlleleTable.itertuples():
        bf = f'{tup2.left_flanking_seq}{tup2.right_flanking_seq}'
        if '*' in bf:
            col_flanking_quality.append( 'Poor' )
        elif bf.upper() != bf:
            col_flanking_quality.append( 'Poor' )
        else:
            col_flanking_quality.append( 'Good' )
            
    AlleleTable['flanking_quality'] = col_flanking_quality
    AlleleTable = AlleleTable[(AlleleTable['flanking_quality']=='Good')].copy()
    
    ### 2. Filter out G/C repeats
    if include_GC_repeats == False:
        AlleleTable = AlleleTable[(AlleleTable['repeat_unit'].isin(['A', 'T']))].copy()
    
    ### 3. Filter out reads without CB or UMI
    AlleleTable.dropna(inplace=True,)
    AlleleTable = AlleleTable[(AlleleTable['reference_STR_allele']<=max_STR_length) & 
                              (AlleleTable['reference_STR_allele']>=min_STR_length)].copy()
    
    AlleleTable['diff'] = AlleleTable['read_STR_allele'] - AlleleTable['reference_STR_allele']

    return AlleleTable

def runscMSI_score( PATH_SCANPY_ADATA, PATH_STR_ALLELE_TABLE, filename, alleleTable_filename, DIR_OUT ):
    
    # Load data
    adata = sc.read_h5ad(PATH_SCANPY_ADATA)
    
    AlleleTable = pd.read_csv(PATH_STR_ALLELE_TABLE, sep='\t')
    AlleleTable = preprocessAlleleTable(AlleleTable)
    AlleleTable.to_csv(f'{DIR_OUT}/{alleleTable_filename}.preprocessed.tsv', sep='\t', index=False)

    # Create Identifier->MSProfile dictionary
    dict_Identifier_To_MSprofile = dict()

    for Identifier, edf in AlleleTable.groupby("Identifier"):
        edf_o = edf['diff'].dropna()
        if len(edf_o) > 0:
            dict_Identifier_To_MSprofile[Identifier] = [ np.mean(edf_o), np.std(edf_o), len(edf_o), -1 * np.mean(edf_o) * np.std(edf_o), ]
            
    for Identifier in adata.obs['Identifier']:
        try: dict_Identifier_To_MSprofile[Identifier]
        except KeyError: dict_Identifier_To_MSprofile[Identifier]=[0, 0, 0, 0]
    
    # [2] Overlay microsatellite information to Scanpy object
    adata.obs['AvgSTRDiff'] = [ dict_Identifier_To_MSprofile[Identifier][0] for Identifier in adata.obs['Identifier'] ]
    adata.obs['StdSTRDiff'] = [ dict_Identifier_To_MSprofile[Identifier][1] for Identifier in adata.obs['Identifier'] ]
    adata.obs['NumSTRLoci'] = [ dict_Identifier_To_MSprofile[Identifier][2] for Identifier in adata.obs['Identifier'] ]
    adata.obs['MSI_score']  = [ dict_Identifier_To_MSprofile[Identifier][3] for Identifier in adata.obs['Identifier'] ]
    adata.write(f'{DIR_OUT}/{filename}.NanoMnT.h5ad')
    
    return 
    
def main():
    start_time  = time.time()

    github_link = "https://github.com/18parkky/scMnT"
    script_description = f"[scMnT] Given a Scanpy object (h5ad) and AlleleTable (tsv), label MSI score to the Scanpy object and save to disk"
    parser = argparse.ArgumentParser(description=script_description)

    # Required parameters
    parser.add_argument('-a', '--PATH_SCANPY_ADATA',      
                        help="PATH to the single-cell Scanpy object with *unique identifier for each cell stored at adata.obs['Identifier']*", 
                        required=True,
                        )
    parser.add_argument('-at', '--PATH_STR_ALLELE_TABLE',      
                        help="PATH to the Allele Table (generated using getAlleleTable) *unique identifier for each cell stored at df['Identifier']*", 
                        required=True,
                        )
    parser.add_argument('-out', '--DIR_OUT',          
                        help='Directory to write output files (default: current directory)', 
                        required=False, 
                        type=str, 
                        default=os.getcwd() 
                        )
    
    args = vars(parser.parse_args())
    
    PATH_SCANPY_ADATA       = args["PATH_SCANPY_ADATA"]
    PATH_STR_ALLELE_TABLE   = args["PATH_STR_ALLELE_TABLE"]
    DIR_OUT                 = args["DIR_OUT"]
    
    filename                = os.path.splitext(os.path.basename(PATH_SCANPY_ADATA))[0]
    alleleTable_filename    = os.path.splitext(os.path.basename(PATH_STR_ALLELE_TABLE))[0]
    
    nanomnt_utility.checkAndCreate(DIR_OUT)
    PATH_log = f"{DIR_OUT}/nanomnt.scMSI_score.{filename}.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'{k}\t:\t{v}')
        
    runscMSI_score( PATH_SCANPY_ADATA, PATH_STR_ALLELE_TABLE, filename, alleleTable_filename, DIR_OUT )
    logging.info(f"Finished scMSI-score.py\t(Total time taken: {nanomnt_utility.getElapsedTime(start_time)} seconds)")

        
if __name__ == "__main__":
    main()
