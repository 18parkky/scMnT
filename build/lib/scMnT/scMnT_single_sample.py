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
                        help="Path to Scanpy object (.h5ad). Cell types must be stored in .obs['CellType'], with tumor cells labeled in the format 'Tumor *' (e.g., Tumor 1, Tumor 2).",
                        required=True,
                        )    
    
    required.add_argument('-b', '--PATH_bam',      
                        help="PATH to input BAM file (must be sorted and indexed)", 
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
    optional.add_argument('-k', '--sample_key',     
                        help="Column in adata.obs that identifies samples. Required when running scMnT on multiple samples (default: None).",
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
    sample_key              = args['sample_key']
    flanking_length         = args['flanking']
    minimum_loci            = args['minimum_loci']
    minimum_MS_length   = args['minimum_MS_length']
    maximum_MS_length   = args['maximum_MS_length']
    threads                 = args['threads']
    DIR_out                 = args["DIR_out"]
        
    # Load packages
    import scMnT.getAlleleTable as getAlleleTable
    import scMnT.scMnT_score as scMnT_score
    import scMnT.scMnT_find as scMnT_find
    import scMnT.scMnT_utility as scMnT_utility

    # Create log file
    scMnT_utility.checkAndCreate(DIR_out)
    PATH_log = f"{DIR_out}/scmnt.{filename}.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"[scMnT] Listing inputs:")
    for k, v in args.items():
        logging.info(f'{k}\t:\t{v}')
        
    mapq_threshold              = 60
    realignment_flanking_length = 50
    
    # (1) Get allele table
    getAlleleTable.runGetAlleleTable( PATH_bam, PATH_str_tsv, PATH_reference_genome, mapq_threshold, threads, flanking_length, realignment_flanking_length, start_time, PATH_log, DIR_out )
    PATH_Allele_Table = f'{DIR_out}/{filename}.AlleleTable.tsv'
    
    # (2) Label MSI scores
    import scanpy as sc
    adata = sc.read_h5ad(PATH_Scanpy_obj)
    if 'CellType' not in list(adata.obs.columns):
        logging.error(f"[scMnT] Cell type information must be available at .obs['CellType']!")
        raise ValueError
    scMnT_score.runscMSI_score( adata, PATH_Allele_Table, minimum_MS_length, maximum_MS_length, filename, DIR_out )

    # (3) Find MSI cell type    
    adata = sc.read_h5ad(f'{DIR_out}/{filename}.scMnT.h5ad')
    scMnT_find.runscMSI_find( adata, minimum_loci, filename, DIR_out )
    
    logging.info(f"[scMnT] Finished scMnT\t(Total time taken: {scMnT_utility.getElapsedTime(start_time)} seconds)")
    
    return 

if __name__ == "__main__":
    main()