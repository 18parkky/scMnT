import os
import time
import pysam
import logging
import argparse
import numpy as np
import pandas as pd

def calcSeqComplexity( sequence, k=2 ):
    try:
        kmer_freq = dict()
        for i in range(0, len(sequence)):
            kmer = sequence[i:i+k]
            if len(kmer) != k: continue

            try:
                kmer_freq[kmer] += 1
            except KeyError:
                kmer_freq[kmer] = 1
        
        seq_complexity = (len( sequence ) - 1)**2 / sum( [n**2 for n in kmer_freq.values()] )
        # return 1 - 1 / seq_complexity
        return seq_complexity
    except:
        return None

def checkAndCreate( PATH ):
    """
    description:    Check if a PATH exists and if not, create PATH
    return:         None
    """
    import os
    if not os.path.exists( PATH ): os.mkdir( PATH )
    
    
def getElapsedTime( start_time, rounding_point=2 ):
    """ 
    description:    Given a starting time, return elapsed time (seconds)
    return:         float
    """
    import time
    return round(time.time() - start_time, rounding_point)

def main():
    start_time  = time.time()

    github_link = "https://github.com/18parkky/NanoMnT"
    script_description = f"Given a Krait output file, filter out STR with low-complex flanking sequences"
    parser = argparse.ArgumentParser(description=script_description)
    
    # Required parameters
    parser.add_argument('-k', '--dir_krait',      
                        help="Directory of Krait output file", 
                        required=True,
                        )
    
    parser.add_argument('-r', '--dir_ref_genome',      
                        help="Directory of reference genome", 
                        required=True,
                        )
    
    parser.add_argument('-rn', '--reference_genome_name',      
                        help="Name of the reference genome", 
                        required=True,
                        )
    
    parser.add_argument('-PATH', '--PATH_out',  
                        help='PATH of the output files (default: current directory)', 
                        required=False, 
                        type=str, 
                        default=os.getcwd(),
                        )
    
    args = vars(parser.parse_args())

    dir_krait       = args["dir_krait"]
    dir_ref_genome  = args["dir_ref_genome"]
    ref_genome_name = args["reference_genome_name"]
    PATH_out        = args["PATH_out"]

    checkAndCreate(PATH_out)

    # Create log file    
    dir_log = f"{PATH_out}/nanomnt.filterKrait.{ref_genome_name}.log"
    logging.basicConfig(filename=dir_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'{k}\t:\t{v}')
        
    PATH_save_raw   = f"{PATH_out}/raw/{ref_genome_name}"
    PATH_save_filt  = f"{PATH_out}/filtered/{ref_genome_name}"

    checkAndCreate(f"{PATH_out}/raw")
    checkAndCreate(f"{PATH_out}/filtered")
    checkAndCreate(PATH_save_raw)
    checkAndCreate(PATH_save_filt)
        
    df_krait = pd.read_csv(dir_krait, sep='\t')
    
    # Label flanking sequences
    ref_genome = pysam.FastaFile( dir_ref_genome )

    flanking_length = 12

    list_lf_col, list_rf_col = list(), list()

    for tup in df_krait.itertuples():
        chromosome  = tup.sequence
        start_pos   = int( tup.start )
        end_pos     = int( tup.end )

        lf = str()
        try:
            for nt in ref_genome.fetch( str(chromosome), start_pos - flanking_length - 1, start_pos - 1 ):
                lf += nt.upper() 
            
            rf = str() 
            for nt in ref_genome.fetch( str(chromosome), end_pos, end_pos + flanking_length ):
                rf += nt.upper()

            list_lf_col.append( lf )
            list_rf_col.append( rf )
        except Exception as e:        
            list_lf_col.append( None )
            list_rf_col.append( None )

    df_krait["left_flanking_seq"] = list_lf_col
    df_krait["right_flanking_seq"] = list_rf_col
    
    # Calculate 2mer complexity
    df_krait["lf_2mer_complexity"] = [ calcSeqComplexity(lf, 2) for lf in df_krait["left_flanking_seq"] ]
    df_krait["rf_2mer_complexity"] = [ calcSeqComplexity(rf, 2) for rf in df_krait["right_flanking_seq"] ]
    
    # Save raw krait 
    for t, edf in df_krait.groupby("type"):
        edf.to_csv(f"{PATH_save_raw}/{ref_genome_name}-{t}bp_STR.tsv.gz", compression="gzip", sep='\t', index=False)
        
    # Save filtered krait 
    for t, edf in df_krait.groupby("type"):
        if t not in [1, 2, 3]: continue
        
        lf_avg, rf_avg = np.mean(edf["lf_2mer_complexity"]), np.mean(edf["rf_2mer_complexity"])
        lf_std, rf_std = np.std(edf["lf_2mer_complexity"]), np.std(edf["rf_2mer_complexity"])
        
        edf_oi = edf[(edf["lf_2mer_complexity"]>=lf_avg + 0.25*lf_std) & (edf["rf_2mer_complexity"]>=rf_avg + 0.25*rf_std)]
        edf_oi.to_csv(f"{PATH_save_filt}/{ref_genome_name}-{t}bp_STR.tsv.gz", compression="gzip", sep='\t', index=False)
    
    logging.info(f"Finished filterKrait.py\t(Total time taken: {getElapsedTime(start_time)} seconds)")

if __name__ == "__main__":
    main()