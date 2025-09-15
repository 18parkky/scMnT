import os
import time
import math 
import glob
import pysam
import logging
import argparse
import subprocess
import numpy as np
import Levenshtein
import pandas as pd
import multiprocessing

import NanoMnT.nanomnt_utility as nanomnt_utility

# Correct the allele motif by motif
def correctAllele(STR_allele_table, err_threshold=None):
    """
    description:    Supplementary figure
    return:         (dictionary, list)
    """
    dict_allele_to_correctedAllele = dict()
    removed = list() 
        
    allele_set  = list( set(STR_allele_table.allele) )
    allele_set.sort(key=len)
    repeat_unit = str(STR_allele_table.iloc[0].repeat_unit)
    length_nt   = int(STR_allele_table.iloc[0].reference_STR_allele) * len(repeat_unit)

    allele_query_range = range( 0, int( len(allele_set[-1])/len(repeat_unit) ) + 5 )

    if err_threshold == None:
        err_threshold = int(length_nt / 2)

    for allele in allele_set:
        
        allele = allele.replace("*", "").upper()
        
        dict_potAllele_to_dist = dict()
        
        for potential_allele_size in allele_query_range:
            potential_allele = repeat_unit * potential_allele_size
            dist = Levenshtein.distance( potential_allele, allele )
            dict_potAllele_to_dist[potential_allele] = dist
            if dist == 0:
                break

        dict_potAllele_to_dist = dict(sorted(dict_potAllele_to_dist.items(), key=lambda x:x[1]))
        for potAllele, dist in dict_potAllele_to_dist.items():
            candidateAllele = potAllele
            minDist     = dist
            break 

        if minDist == 0:
            pass 
        else:
            num_candidate_allele = 0 
            dict_candidateAllele_to_lengthDiff = dict()
            for potAllele, dist in dict_potAllele_to_dist.items():
                if dist == minDist:
                    num_candidate_allele += 1
                    dict_candidateAllele_to_lengthDiff[potAllele] = abs( len( allele ) - len( potAllele ) )
                    
            # if there is only one candidate allele, choose that allele
            if num_candidate_allele == 1:
                pass 
            else:
                # if there is more than one, choose the allele that has the closest length (e.g., if Observed allele = AAGAA and candidate alleles are [Ax4, Ax5], choose Ax5)
                dict_candidateAllele_to_lengthDiff = dict(sorted(dict_candidateAllele_to_lengthDiff.items(), key=lambda x:x[1]))
                for potAllele, lengthDiff in dict_candidateAllele_to_lengthDiff.items():
                    candidateAllele = potAllele
                    break

        if minDist > err_threshold:
            removed.append(allele) 
        else:
            dict_allele_to_correctedAllele[allele] = ( candidateAllele, minDist )

    return (dict_allele_to_correctedAllele, removed)


def multiprocess( STR_allele_table_chunk, distance, process_num, PATH_multiprocess_out ):
    
    df_corrected_alleleTable = []
    for locus, edf in STR_allele_table_chunk.groupby('locus'):
        
        dict_allele_to_correctedAllele, removed = correctAllele(edf, distance)

        # Fill in a_to_ca with removed alleles
        for removed_allele in removed:
            dict_allele_to_correctedAllele[removed_allele] = ( "uncor", -1 )
        
        allele_list = [ a.replace("*", "").upper() for a in edf.allele ]
        edf['corrected_allele'] = [dict_allele_to_correctedAllele[a][0] for a in allele_list]
        edf['editing_distance'] = [dict_allele_to_correctedAllele[a][1] for a in allele_list]
        
        df_corrected_alleleTable.append(edf)
    
    df_corrected_alleleTable = pd.concat(df_corrected_alleleTable) 
    df_corrected_alleleTable["read_STR_length"] = [ tup.corrected_allele.replace("*", "").upper().count(tup.repeat_unit) for tup in df_corrected_alleleTable.itertuples() ]
    df_corrected_alleleTable.to_csv(f"{PATH_multiprocess_out}/corrected_allele_table.{process_num}.tsv", sep='\t', index=False)

    return 

def run(dir_STR_allele_table, distance, threads, start_time, dir_log, PATH_out):
    logging.basicConfig(filename=dir_log, level=logging.INFO)

    at_filename = os.path.splitext(os.path.basename(dir_STR_allele_table))[0]
    
    # (1) Create folder for storing temporary files
    PATH_multiprocess_out = f"{PATH_out}/{at_filename}.multiprocessing_temp"
    nanomnt_utility.checkAndCreate( PATH_out )
    nanomnt_utility.checkAndCreate( PATH_multiprocess_out )
        
    # (2) Load and chunk allele table for multiprocessing
    STR_allele_table = pd.read_csv( dir_STR_allele_table, sep='\t' )
        
    try:
        chunk_size = math.ceil( len(STR_allele_table) / threads )
        STR_allele_table_chunks = [STR_allele_table[i:i+chunk_size] for i in range(0,STR_allele_table.shape[0],chunk_size)]
    except ValueError: # Occurs when chunk size < threads
        logging.error(f"Number of reads ({len(STR_allele_table)}) is less than number of threads to use! Try reducing the number of threads.\t(elapsed time: {nanomnt_utility.getElapsedTime(start_time)} seconds)")
        raise ValueError
    
    # (4) Multiprocess using multiprocessing
    processes = list()
    logging.info(f"Starting multiprocessing with {threads}")    
    for idx, STR_allele_table_chunk in enumerate(STR_allele_table_chunks):
        
        p = multiprocessing.Process(target=multiprocess,
                                    args=[STR_allele_table_chunk, distance, idx+1, PATH_multiprocess_out] )
        p.start()
        processes.append(p)
    for process in processes:
        process.join()
        
    # (5) Collect multiprocessing results & save to disk   
    list_dir_allele_tables = glob.glob( f"{PATH_multiprocess_out}/corrected_allele_table.*.tsv" )
    df_concat_allele_table = pd.concat( [pd.read_csv(dir_allele_table, sep='\t') for dir_allele_table in list_dir_allele_tables] )
    df_concat_allele_table.to_csv(f"{PATH_out}/{at_filename}.corrected.tsv", sep='\t', index=False)
    
    for dir_allele_table in list_dir_allele_tables:
        os.remove( dir_allele_table ) 
    os.rmdir( PATH_multiprocess_out )
    
    

def main():
    start_time  = time.time()

    github_link = "https://github.com/18parkky/NanoMnT"
    script_description = f"Given an Allele table, perform STR error correction for each read. For more information, see GitHub: {github_link}"
    parser = argparse.ArgumentParser(description=script_description)
    
    # Required parameters
    parser.add_argument('-at', '--dir_STR_allele_table',  
                        help="Directory of STR allele table",
                        required=True
                        )
    parser.add_argument('-e', '--distance',      
                        help="Allowed editing (Levenshtein) distance threshold (default: half of STR length)",
                        required=False, 
                        type=int, 
                        default=None,
                        )
    parser.add_argument('-t', '--threads',      
                        help="Number of threads (default: 4)",
                        required=False, 
                        type=int, 
                        default=4,
                        )
    parser.add_argument('-PATH', '--PATH_out',  
                        help='PATH of the output files (default: current directory)', 
                        required=False, 
                        type=str, 
                        default=os.getcwd(),
                        )

    args = vars(parser.parse_args())

    dir_STR_allele_table        = args["dir_STR_allele_table"]
    distance            = args["distance"]
    threads             = args["threads"]
    PATH_out            = args["PATH_out"]

    at_filename = os.path.splitext(os.path.basename(dir_STR_allele_table))[0]
    
    # Create log file
    nanomnt_utility.checkAndCreate(PATH_out)
    dir_log = f"{PATH_out}/nanomnt.correctAllele.{at_filename}.log"
    logging.basicConfig(filename=dir_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'{k}\t:\t{v}')
        
    run(dir_STR_allele_table, distance, threads, start_time, dir_log, PATH_out)
    logging.info(f"Finished correctAllele.py\t(Total time taken: {nanomnt_utility.getElapsedTime(start_time)} seconds)")
        
if __name__ == "__main__":
    main()