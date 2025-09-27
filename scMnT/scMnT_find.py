import os
import time
import logging
import argparse 

import scipy.stats
import numpy as np
import pandas as pd
import scanpy as sc

import scmnt_utility

def cliffs_delta(x, y):
    nx = len(x)
    ny = len(y)
    n_greater = sum(i > j for i in x for j in y)
    n_less = sum(i < j for i in x for j in y)
    delta = (n_greater - n_less) / (nx * ny)
    return delta

def runscMSI_find( adata, mininum_loci, filename, DIR_out ):
    
    # Load data
    adata_obs = adata.obs[(adata.obs['NumSTRLoci']>=mininum_loci)].copy()

    NormalCellTypes = [ CellType for CellType in set(adata_obs['CellType']) if CellType[:5] != 'Tumor' ]
    logging.info(f'Normal cell types:\t{", ".join(NormalCellTypes)}')
    NormalMSIscores = adata_obs[adata_obs['CellType'].isin( NormalCellTypes )]['MSI_score']

    TestResults = list()
    for Tumor_CellType, edf in adata_obs[~(adata_obs['CellType'].isin(NormalCellTypes))].groupby('CellType', observed=True):
        MSI_scores = edf['MSI_score']
        
        sample_n = min([len(NormalMSIscores), len(MSI_scores)])
        
        stat, pval = scipy.stats.ks_2samp(MSI_scores.sample(sample_n, random_state=42), NormalMSIscores.sample(sample_n, random_state=42))
        delta = cliffs_delta(MSI_scores, NormalMSIscores)

        TestResults.append( [Tumor_CellType, pval, delta, sample_n] )
        
    TestResults = pd.DataFrame(TestResults, columns=['CT', 'pval', 'delta', 'n_cells'])
    TestResults['logPval'] = [ -np.log(pval+10**-10) for pval in TestResults['pval'] ]
    
    TestResults.to_csv(f'{DIR_out}/{filename}.findMSI_results.tsv', sep='\t', index=False)

def main():
    start_time  = time.time()

    github_link = "https://github.com/18parkky/scMnT"
    script_description = f"[scMnT] Given a MSI score-labeled Scanpy object (output of scMSI-score.py), find MSI cells"
    parser = argparse.ArgumentParser(description=script_description)

    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    
    # Required parameters
    required.add_argument('-a', '--PATH_Scanpy_obj',      
                        help="Path to Scanpy object (.h5ad) processed by scMnT_score.py",
                        required=True,
                        )
    
    required.add_argument('-n', '--filename',      
                        help="Name of the output files", 
                        required=True,
                        )
    
    optional.add_argument('-m', '--mininum_loci',     
                        help="Minimum number of MS loci per cell required for MSI/MSS identification (default: 10). Cells with fewer loci are excluded.",
                        required=False, 
                        type=int,
                        default=10,
                        )
    
    optional.add_argument('-out', '--DIR_out',          
                        help='Directory to write output files (default: current directory)', 
                        required=False, 
                        type=str, 
                        default=os.getcwd() 
                        )
    
    args = vars(parser.parse_args())
    
    PATH_Scanpy_obj = args['PATH_Scanpy_obj']
    filename        = args['filename']
    mininum_loci    = args['mininum_loci']
    DIR_out         = args["DIR_out"]
        
    scmnt_utility.checkAndCreate(DIR_out)
    PATH_log = f"{DIR_out}/scMnT.scMSI_find.{filename}.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'{k}\t:\t{v}')

    # Load data
    adata = sc.read_h5ad(PATH_Scanpy_obj)
    
    runscMSI_find( adata, mininum_loci, filename, DIR_out )
    logging.info(f"Finished scMnT-score.py\t(Total time taken: {scmnt_utility.getElapsedTime(start_time)} seconds)")

        
if __name__ == "__main__":
    main()
