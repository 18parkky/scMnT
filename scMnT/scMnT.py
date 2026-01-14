import os
import time
import logging
import argparse

def main():
    start_time  = time.time()

    github_link = "https://github.com/18parkky/scMnT"
    script_description = f"[scMnT] Using (1) a CB/UB-tagged BAM, (2) a cell type-annotated Scanpy object, and (3) a reference genome FASTA, this script generates a read-level microsatellite information table and scores cells by MSI. For details, see GitHub: {github_link}"
    parser = argparse.ArgumentParser(description=script_description)
    
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    
    # Required parameters
    required.add_argument('-a', '--PATH_Scanpy_obj',      
                        help="Path to Scanpy object (.h5ad). Cell types must be stored in .obs['CellType'], with tumor cells labeled in the format 'Tumor *' (e.g., Tumor 1, Tumor 2)",
                        required=True,
                        )    
    
    required.add_argument('-s', '--PATH_str_tsv',  
                        help="PATH to STR list file generated using Krait", 
                        required=True,
                        )    
    
    required.add_argument('-r', 
                        "--PATH_reference_genome", 
                        help="PATH to reference genome (must be same as the one used for generating BAM)",
                        required=True,
                        )
    
    required.add_argument('-n', '--filename',      
                        help="Name of the output files", 
                        required=True,
                        )
    
    # Optional parameters
    optional.add_argument('-b', '--PATH_bam',      
                        help="PATH to input BAM file (must be sorted and indexed, required if sample sheet is not provided)", 
                        required=False,
                        )
    
    optional.add_argument('-ss', '--PATH_sample_sheet',     
                        help="Tab-delimited file with two columns: (1) sample name and (2) path to its BAM file. Don't include column name. Required when running scMnT on multiple samples",
                        required=False, 
                        type=str,
                        default=None,
                        )
    
    optional.add_argument('-f', '--flanking',     
                        help="Length of flanking sequences to collect for each read (default: 6)", 
                        required=False, 
                        type=int,
                        default=6,
                        )
    
    optional.add_argument('-m', '--minimum_loci',     
                        help="Minimum number of MS loci per cell required for MSI/MSS identification (default: 10). Cells with fewer loci are excluded.",
                        required=False, 
                        type=int,
                        default=10,
                        )
    
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
    
    optional.add_argument('--resume',
                        action='store_true',
                        help="Resume from the previous run (default: False)",
                        )
    
    # Multiple sample
    # optional.add_argument('-sk', '--sample_key',     
    #                     help="Column in adata.obs that identifies samples. Required when running scMnT on multiple samples (default: None).",
    #                     required=False, 
    #                     type=str,
    #                     default=None,
    #                     )
    
    optional.add_argument('-t', '--threads',      
                        help="Number of threads (default: 4)",
                        required=False, 
                        type=int, 
                        default=4,
                        )
    
    optional.add_argument('-out', '--DIR_out',
                        help='Directory to write output files (default: current directory)',
                        required=False,
                        type=str,
                        default=os.getcwd(),
                        )
    
    args = vars(parser.parse_args())

    PATH_Scanpy_obj         = args['PATH_Scanpy_obj']
    PATH_bam                = args['PATH_bam']
    PATH_str_tsv            = args['PATH_str_tsv']
    PATH_reference_genome   = args['PATH_reference_genome']
    filename                = args['filename']
    flanking_length         = args['flanking']
    minimum_loci            = args['minimum_loci']
    minimum_MS_length   = args['minimum_MS_length']
    maximum_MS_length   = args['maximum_MS_length']
    resume              = args['resume']
    # sample_key              = args['sample_key']  # Fix to SampleID
    PATH_sample_sheet       = args['PATH_sample_sheet']
    threads                 = args['threads']
    DIR_out                 = args["DIR_out"]
        
    # Load packages
    import scMnT.getAlleleTable as getAlleleTable
    import scMnT.scMnT_score as scMnT_score
    import scMnT.scMnT_find as scMnT_find
    import scMnT.scMnT_utility as scMnT_utility

    # Create log file
    scMnT_utility.checkAndCreate(DIR_out)
    PATH_log = f"{DIR_out}/scMnT.{filename}.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"[scMnT] Listing inputs:")
    for k, v in args.items():
        logging.info(f'  {k}\t:\t{v}')
        
    mapq_threshold              = 1
    realignment_flanking_length = 50
    
    import pandas as pd 
    import scanpy as sc
    
    if PATH_sample_sheet==None: # Single sample
        logging.info(f"[scMnT] Single sample")
        # (1) Get allele table
        PATH_Allele_Table = f'{DIR_out}/{filename}.AlleleTable.tsv'
        if os.path.exists(PATH_Allele_Table)==False:
            getAlleleTable.runGetAlleleTable( PATH_bam, PATH_str_tsv, PATH_reference_genome, mapq_threshold, threads, flanking_length, realignment_flanking_length, start_time, filename, DIR_out )
        else:
            if resume==True:
                logging.info(f"Allele table found in: {PATH_Allele_Table}. Resuming previous run")
            else:
                logging.error(f"Allele table found in: {PATH_Allele_Table}! Use --resume parameter to continue from previous run")
                raise ValueError
            
        # (2) Label MSI scores
        adata = sc.read_h5ad(PATH_Scanpy_obj)
        PATH_MSI_score_labeled_adata = f'{DIR_out}/{filename}.scMnT.h5ad'
        if os.path.exists(PATH_MSI_score_labeled_adata)==False:
            scMnT_score.runscMSI_score( adata, PATH_Allele_Table, minimum_MS_length, maximum_MS_length, filename, DIR_out )
        else:
            if resume==True:
                logging.info(f"MSI score-labeled AnnData found in: {PATH_MSI_score_labeled_adata}. Resuming previous run")
            else:
                logging.error(f"MSI score-labeled AnnData found in: {PATH_MSI_score_labeled_adata}! Use --resume parameter to continue from previous run")
                raise ValueError

        # (3) Find MSI cell type    
        adata = sc.read_h5ad(PATH_MSI_score_labeled_adata)
        if 'CellType' not in list(adata.obs.columns):
            logging.error(f"[scMnT] Cell type information must be available at .obs['CellType']!")
            raise ValueError
        scMnT_find.runscMSI_find( adata, minimum_loci, filename, DIR_out )
        
        logging.info(f"[scMnT] Finished scMnT\t(Total time taken: {scMnT_utility.getElapsedTime(start_time)} seconds)")
    
    else:
        logging.info(f"[scMnT] Multiple samples")
        
        PATH_Allele_Table_merged = f'{DIR_out}/{filename}.AlleleTable.tsv'
        # (1) Run getAlleleTable for each sample
        if os.path.exists(PATH_Allele_Table_merged)==False:
        
            Allele_Tables_merged = list()
            sample_sheet = pd.read_csv(PATH_sample_sheet, sep='\t', header=None)
            sample_sheet.columns = ['SampleID', 'PATH_BAM']
            
            for tup in sample_sheet.itertuples():
                DIR_out_e = f'{DIR_out}/{tup.SampleID}'
                scMnT_utility.checkAndCreate(DIR_out_e)
                
                PATH_Allele_Table_e = f'{DIR_out_e}/{tup.SampleID}.AlleleTable.tsv'
                if os.path.exists(PATH_Allele_Table_e)==False:
                    getAlleleTable.runGetAlleleTable( tup.PATH_BAM, PATH_str_tsv, PATH_reference_genome, mapq_threshold, threads, flanking_length, realignment_flanking_length, start_time, tup.SampleID, DIR_out_e )
                else:
                    if resume==True:
                        logging.info(f"Allele table for {tup.SampleID} found in: {PATH_Allele_Table_e}. Resuming previous run")
                    else:
                        logging.error(f"Allele table for {tup.SampleID} found in: {PATH_Allele_Table_e}! Use --resume parameter to continue from previous run")
                        raise ValueError            
                    
                # Label sample 
                Allele_Table_e = pd.read_csv(PATH_Allele_Table_e, sep='\t')
                Allele_Table_e['SampleID'] = tup.SampleID
                Allele_Tables_merged.append( Allele_Table_e )
            
            Allele_Tables_merged = pd.concat(Allele_Tables_merged)
            Allele_Tables_merged.reset_index(inplace=True, drop=True)
            Allele_Tables_merged.to_csv(PATH_Allele_Table_merged, sep='\t', index=False)
        
        else:
            if resume==True:
                logging.info(f"Merged Allele Table found in: {PATH_Allele_Table_merged}. Resuming previous run")
            else:
                logging.error(f"Merged Allele Table found in: {PATH_Allele_Table_merged}! Use --resume parameter to continue from previous run")
                raise ValueError
            
        # (2) Run scMnT score
        adata = sc.read_h5ad(PATH_Scanpy_obj)
        PATH_MSI_score_labeled_adata = f'{DIR_out}/{filename}.scMnT.h5ad'
        if os.path.exists(PATH_MSI_score_labeled_adata)==False:
            scMnT_score.runscMSI_score_ms( adata, PATH_Allele_Table_merged, minimum_MS_length, maximum_MS_length, filename, DIR_out )
        else:
            if resume==True:
                logging.info(f"MSI score-labeled AnnData found in: {PATH_MSI_score_labeled_adata}. Resuming previous run")
            else:
                logging.error(f"MSI score-labeled AnnData found in: {PATH_MSI_score_labeled_adata}! Use --resume parameter to continue from previous run")
                raise ValueError
        
        # (3) Find MSI cell type    
        adata = sc.read_h5ad(PATH_MSI_score_labeled_adata)
        if 'CellType' not in list(adata.obs.columns):
            logging.error(f"[scMnT] Cell type information must be available at .obs['CellType']!")
            raise ValueError
        scMnT_find.runscMSI_find( adata, minimum_loci, filename, DIR_out )
        
        logging.info(f"[scMnT] Finished scMnT\t(Total time taken: {scMnT_utility.getElapsedTime(start_time)} seconds)")

    return 

if __name__ == "__main__":
    main()
