import os
import time
import logging
import argparse 
import numpy as np
import pandas as pd
import scanpy as sc
import scmnt_utility

def preprocessAlleleTable( AlleleTable, min_STR_length, max_STR_length, include_GC_repeats=False ):
    
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
    AlleleTable = AlleleTable[(AlleleTable['reference_STR_allele']<=max_STR_length) & (AlleleTable['reference_STR_allele']>=min_STR_length)].copy()
    
    AlleleTable['diff'] = AlleleTable['read_STR_allele'] - AlleleTable['reference_STR_allele']

    return AlleleTable

def runscMSI_score( adata, PATH_Allele_Table, minimum_MS_length, maximum_MS_length, filename, DIR_out ):
    
    AlleleTable = pd.read_csv(PATH_Allele_Table, sep='\t')
    AlleleTable = preprocessAlleleTable(AlleleTable, minimum_MS_length, maximum_MS_length)
    AlleleTable.to_csv(f'{DIR_out}/{filename}.AlleleTable.preprocessed.tsv', sep='\t', index=False)

    if 'CB' not in list(adata.obs.columns):
        adata.obs['CB'] = [ idx.split('-')[0] for idx in adata.obs.index ]
        
    # [1] Create CB->MSProfile dictionary
    dict_CB_to_MSprofile = dict()

    for CB, edf in AlleleTable.groupby("CB"):
        edf_o = edf['diff'].dropna()
        if len(edf_o) > 0:
            dict_CB_to_MSprofile[CB] = [ np.mean(edf_o), np.std(edf_o), len(edf_o), -1 * np.mean(edf_o) * np.std(edf_o), ]

    # [2] Overlay microsatellite information to Scanpy object
    MS_columns = ['AvgSTRDiff', 'StdSTRDiff', 'NumSTRLoci', 'MSI_score']
    # Retrieve MS profiles for each cell, defaulting to [0, 0, 0, 0] if CB not found
    MS_profile = [dict_CB_to_MSprofile.get(cb, [0, 0, 0, 0]) for cb in adata.obs['CB']]
    for i, col in enumerate(MS_columns):
        adata.obs[col] = [vals[i] for vals in MS_profile]
    adata.write(f'{DIR_out}/{filename}.scMnT.h5ad')

    return 
    
def main():
    start_time  = time.time()

    github_link = "https://github.com/18parkky/scMnT"
    script_description = f"[scMnT] Given a cell type-annotated Scanpy object and its AlleleTable (generated using getAlleleTable.py), label MSI score to the Scanpy object. For details, see GitHub: {github_link}"
    parser = argparse.ArgumentParser(description=script_description)

    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    
    # Required parameters
    required.add_argument('-a', '--PATH_Scanpy_obj',      
                        help="Path to Scanpy object (.h5ad) with cell types in .obs['CellType']. Tumor cells must be labeled 'Tumor'.",
                        required=True,
                        )
    
    required.add_argument('-at', '--PATH_Allele_Table',      
                        help="PATH to the Allele Table (generated using getAlleleTable.py)", 
                        required=True,
                        )
    
    required.add_argument('-n', '--filename',      
                        help="Name of the output files", 
                        required=True,
                        )
    
    # Optional parameters
    optional.add_argument('-min_l', '--minimum_MS_length',          
                        help='Minimum length of the MS loci to assess (default: 10)', 
                        required=False, 
                        type=int, 
                        default=10 
                        )

    optional.add_argument('-max_l', '--maximum_MS_length',          
                        help='Maximum length of the MS loci to assess (default: 24)', 
                        required=False, 
                        type=int, 
                        default=24 
                        )
        
    optional.add_argument('-out', '--DIR_OUT',          
                        help='Directory to write output files (default: current directory)', 
                        required=False, 
                        type=str, 
                        default=os.getcwd() 
                        )
    
    args = vars(parser.parse_args())
    
    PATH_Scanpy_obj     = args['PATH_Scanpy_obj']
    PATH_Allele_Table   = args['PATH_Allele_Table']
    filename            = args['filename']
    minimum_MS_length   = args['minimum_MS_length']
    maximum_MS_length   = args['maximum_MS_length']
    DIR_out             = args['DIR_out']
    
    scmnt_utility.checkAndCreate(DIR_out)
    PATH_log = f"{DIR_out}/scMnT.scMSI_score.{filename}.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'{k}\t:\t{v}')

    # Load data
    adata = sc.read_h5ad(PATH_Scanpy_obj)
    if 'CellType' not in list(adata.obs.columns):
        logging.error(f"Cell type information must be available at .obs['CellType']!")
        raise ValueError
        
    runscMSI_score( adata, PATH_Allele_Table, minimum_MS_length, maximum_MS_length, filename, DIR_out )
    logging.info(f"Finished scMnT_score.py\t(Total time taken: {scmnt_utility.getElapsedTime(start_time)} seconds)")

        
if __name__ == "__main__":
    main()
