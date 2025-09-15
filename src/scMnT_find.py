import os
import time
import logging
import argparse 

import scipy.stats
import numpy as np
import pandas as pd
import scanpy as sc

import NanoMnT.nanomnt_utility as nanomnt_utility

def cliffs_delta(x, y):
    nx = len(x)
    ny = len(y)
    n_greater = sum(i > j for i in x for j in y)
    n_less = sum(i < j for i in x for j in y)
    delta = (n_greater - n_less) / (nx * ny)
    return delta

def runscMSI_find( PATH_NANOMNT_SCANPY_ADATA, NORMAL_CELLTYPES, MINIMUM_LOCI, filename, DIR_OUT ):
    
    # Load data
    adata = sc.read_h5ad(PATH_NANOMNT_SCANPY_ADATA)
    df = adata.obs[(adata.obs['NumSTRLoci']>=MINIMUM_LOCI)].copy()

    NormalMSIscores = df[df['CellType'].isin(NORMAL_CELLTYPES)]['MSI_score']

    TestResults = list()
    for CellType, edf in df[~(df['CellType'].isin(NORMAL_CELLTYPES))].groupby('CellType', observed=True):
        MSI_scores = edf['MSI_score']
        
        sample_n = min([len(NormalMSIscores), len(MSI_scores)])
        
        stat, pval = scipy.stats.ks_2samp(MSI_scores.sample(sample_n, random_state=42), NormalMSIscores.sample(sample_n, random_state=42))
        delta = cliffs_delta(MSI_scores, NormalMSIscores)

        TestResults.append( [CellType, pval, delta, sample_n] )
        
    TestResults = pd.DataFrame(TestResults, columns=['CT', 'pval', 'delta', 'n_cells'])
    TestResults['logPval'] = [ -np.log(pval+10**-10) for pval in TestResults['pval'] ]
    
    TestResults.to_csv(f'{DIR_OUT}/{filename}.findMSI_results.tsv', sep='\t', index=False)

    # pval_threshold  = 0.01 
    # delta_threshold = 0.4

    # MSI_celltypes = list()

    # TestResults_significant = TestResults[(TestResults['pval']<=pval_threshold) & 
    #                                     ((TestResults['delta']>=delta_threshold) | (TestResults['delta']<=-delta_threshold))]

    # for tup in TestResults_significant.itertuples():
    #     MSI_celltypes.append( tup.CT )

    # TestResults_significant.reset_index(inplace=True, drop=True)
    # TestResults_significant

def main():
    start_time  = time.time()

    github_link = "https://github.com/18parkky/NanoMnT"
    script_description = f"Given a MSI score-labeled Scanpy object (output of scMSI-score.py), find MSI cells"
    parser = argparse.ArgumentParser(description=script_description)

    # Required parameters
    parser.add_argument('-a', '--PATH_NANOMNT_SCANPY_ADATA',      
                        help="PATH to the single-cell Scanpy object with MSI score labeled (output of scMSI-score.py) *Cell group (e.g., cell type) must be pre-computed and accessible at adata.obs['CellGroup']", 
                        required=True,
                        )

    parser.add_argument('-n', '--NORMAL_CELLTYPES',      
                        help="Comma separated list of normal cell types (cells that are surely MSS, such as monocytes) to use as reference (e.g., monocyte,T,fibroblast)", 
                        required=True, type=str,
                        )
    
    parser.add_argument('-m', '--MINIMUM_LOCI',     
                        help="Minimum number of MS loci required per cell to be included in the analysis (default: 10). Cells with fewer loci will be excluded.", 
                        required=False, 
                        type=int,
                        default=10,
                        )
    
    parser.add_argument('-out', '--DIR_OUT',          
                        help='Directory to write output files (default: current directory)', 
                        required=False, 
                        type=str, 
                        default=os.getcwd() 
                        )
    
    args = vars(parser.parse_args())
    
    PATH_NANOMNT_SCANPY_ADATA   = args["PATH_NANOMNT_SCANPY_ADATA"]
    NORMAL_CELLTYPES            = [CT.strip() for CT in args["NORMAL_CELLTYPES"].split(',')]
    MINIMUM_LOCI                = args['MINIMUM_LOCI']
    DIR_OUT                     = args["DIR_OUT"]
    
    filename = os.path.splitext(os.path.basename(PATH_NANOMNT_SCANPY_ADATA))[0]
    
    nanomnt_utility.checkAndCreate(DIR_OUT)
    PATH_log = f"{DIR_OUT}/nanomnt.scMSI_find.{filename}.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'{k}\t:\t{v}')
        
    runscMSI_find( PATH_NANOMNT_SCANPY_ADATA, NORMAL_CELLTYPES, MINIMUM_LOCI, filename, DIR_OUT )
    logging.info(f"Finished scMSI-score.py\t(Total time taken: {nanomnt_utility.getElapsedTime(start_time)} seconds)")

        
if __name__ == "__main__":
    main()