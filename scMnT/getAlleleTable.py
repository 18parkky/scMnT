import gc
import os
import gzip
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

import seaborn as sns 
import matplotlib.pyplot as plt

import NanoMnT.nanomnt_utility as nanomnt_utility

PILEUP_LIMIT = 5000000 # pysam's pileup's default value for max_depth is 8000

def correctAllele(STR_allele_table, dist_threshold=None):
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

    if dist_threshold == None:
        dist_threshold = int(length_nt / 2)

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

        if minDist > dist_threshold:
            removed.append(allele) 
        else:
            dict_allele_to_correctedAllele[allele] = ( candidateAllele, minDist )

    return (dict_allele_to_correctedAllele, removed)

def getFlankingSequence( pysam_fasta_obj, contig, start, end, flanking_length ):

    try:
        lf = str()
        for nt in pysam_fasta_obj.fetch( contig, start - flanking_length - 1, start - 1 ):
            lf += nt.upper() 
        
        rf = str() 
        for nt in pysam_fasta_obj.fetch( contig, end, end + flanking_length ):
            rf += nt.upper()
    except ValueError:
        lf, rf = None, None
        
    return (lf, rf)


# def load_SSR_table( PATH_str_tsv, required_columns, PATH_log=None ):
#     """ 
#     Determine whether the given Krait (or pytrf) file is appropriate.
#     """
#     STR_input_type = None
#     try:
#         ssr_tsv = open( PATH_str_tsv, "r" )
#     except:
#         ssr_tsv = gzip.open( PATH_str_tsv, "r" )
        
#     for line in ssr_tsv:
#         columns = line.strip().split("\t")
#         if len( set(columns).intersection( set(required_columns) ) ) == 0:                    
#             STR_input_type = "pytrf"
#             break
#         else:
#             STR_input_type = "krait"
#             break
                
#     if STR_input_type == None: # Raise error if STR_input_type wasn't determinedß
#         raise ValueError
#     if PATH_log != None:
#         logging.basicConfig(filename=PATH_log, level=logging.INFO)
#         logging.info( f"SSR table identified as:\t{STR_input_type}-generated" )

#     if STR_input_type == "pytrf":
#         STR_table = pd.read_csv(PATH_str_tsv, sep='\t', names=required_columns) 
#     elif STR_input_type == "krait":
#         STR_table = pd.read_csv(PATH_str_tsv, sep='\t')
    
#     ssr_tsv.close()
    
#     STR_table = STR_table.sample( frac=1, random_state=0 ).reset_index(drop=True) # shuffle STR table for even distribution of loci among chunks
#     return STR_table


def realign( bamfile, region, PATH_reference_genome, flanking_length, mapq_threshold, fasta_flanking_length, placeholder_length, sc, PATH_temp, ):

    # (1) Collect reads that align to the given STR region
    dict_reads = dict()
    for pileupcolumn in bamfile.pileup(region[0], region[3] - 1 - flanking_length, region[4] + flanking_length, truncate=True, min_base_quality = 0, min_mapping_quality = mapq_threshold, max_depth=PILEUP_LIMIT): #!#!#!#!#!#!#!#!#!#! min_base_quality 1 -> 0
        for pileupread in pileupcolumn.pileups:
            if (pileupread.alignment.is_secondary == False and pileupread.alignment.is_supplementary == False):
                if pileupread.alignment.query_name not in dict_reads.keys():
                    
                    if sc == True:
                        
                        try:
                            CB = pileupread.alignment.get_tag("CB")
                        except: 
                            CB = None 
                            
                        try:
                            UMI = pileupread.alignment.get_tag("UB")
                        except: 
                            UMI = None      
                            
                    else:
                        CB  = None 
                        UMI = None
                        
                    dict_reads[pileupread.alignment.query_name] = [pileupread.alignment.get_forward_sequence(), CB, UMI]
    
    # (2) Create temporary Fasta file that contains the given STR region sequence and its flankings
    STR_name = f"{region[0]}_{region[1]}x{region[2]}_{region[3]}"
    pysam_genome = pysam.FastaFile(PATH_reference_genome)
    with open(f"{PATH_temp}/{STR_name}.ref.fasta", "w") as fasta:
        lf, rf  = getFlankingSequence( pysam_genome, region[0], region[3], region[4], fasta_flanking_length )
        placeholder_sequence = region[1] * placeholder_length
        fasta.write(f">{region[0]}\n")
        fasta.write(f"{lf}{placeholder_sequence}{rf}\n")
    
    # (3) Create temporary Fasta file for reads that align to this STR region
    with open(f"{PATH_temp}/{STR_name}.reads.fasta", "w") as fasta:
        for readname, readinfo in dict_reads.items():
            sequence        = readinfo[0]
            fasta.write(f">{readname}\n")
            fasta.write(f"{sequence}\n")
            
    # (4) Align using minimap2
    STR_name = f"{region[0]}_{region[1]}x{region[2]}_{region[3]}"
    realigned_bam_out = f"{PATH_temp}/{STR_name}.sorted.bam"
    command = f"minimap2 -a -A 4 -B 10 -t 1 {PATH_temp}/{STR_name}.ref.fasta {PATH_temp}/{STR_name}.reads.fasta | samtools view -Sb -F0x900 | samtools sort > {realigned_bam_out}"
    subprocess.call(command, shell=True)
    pysam.index(realigned_bam_out)
 
    # (5) If sc==True, set CB and UMI tags to realigned BAM
    if sc == True:
        realigned_bamfile = pysam.AlignmentFile(realigned_bam_out, "rb")
        realigned_tagged_bam_out = f"{PATH_temp}/{STR_name}.sorted.tagged.bam"
        realigned_tagged_bamfile = pysam.AlignmentFile(realigned_tagged_bam_out, "wb", template=realigned_bamfile)
        for read in realigned_bamfile.fetch():
            read_name = read.query_name
            read.set_tag( "CB", dict_reads[read_name][1] )
            read.set_tag( "UB", dict_reads[read_name][2] )
            
            realigned_tagged_bamfile.write(read)
            
        realigned_tagged_bamfile.close()

        # Delete realigned BAM
        os.remove( realigned_bam_out )
        os.rename( realigned_tagged_bam_out, realigned_bam_out)
        pysam.index(realigned_bam_out)

    # Delete temporary FASTA files
    tempFiles = glob.glob(f"{PATH_temp}/{STR_name}*.fasta")
    for file in tempFiles:
        os.remove( file )
        
    return  

def errorCorrection_multiprocess( STR_allele_table_chunk, distance, process_num, PATH_multiprocess_out ):
    
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
    df_corrected_alleleTable["read_STR_allele"] = [ tup.corrected_allele.replace("*", "").upper().count(tup.repeat_unit) for tup in df_corrected_alleleTable.itertuples() ]
    df_corrected_alleleTable.to_csv(f"{PATH_multiprocess_out}/corrected_allele_table.{process_num}.tsv", sep='\t', index=False)

    return 

def extractSTR_multiprocess( PATH_bamfile, STR_table, PATH_reference_genome, flanking_length, realignment_f, mapq_threshold, sc, PATH_temp ):
    
    for tup in STR_table.itertuples():
        
        # (1) Realign reads using minimap2 via subprocess
        STR_region = [ str(tup.sequence), tup.motif, tup.repeat, tup.start, tup.end ]
        fasta_flanking_length = realignment_f
        bamfile = pysam.AlignmentFile(PATH_bamfile, "rb")
        placeholder_length = 50
        realign(bamfile, STR_region, PATH_reference_genome, flanking_length, mapq_threshold, fasta_flanking_length, placeholder_length, sc, PATH_temp )
        
        # (2) Get the realigned STR information of each read
        STR_name = f"{STR_region[0]}_{STR_region[1]}x{STR_region[2]}_{STR_region[3]}"
        realigned_bamfile = pysam.AlignmentFile(f"{PATH_temp}/{STR_name}.sorted.bam", "rb")
        
        dict_readname_to_STRinfo = dict() # key: read_name, value: [STR_seq, base quality, int_left_padding, int_right_padding, flag]
        bool_firstIteration = False
        
        pileup_start_pos    = fasta_flanking_length - flanking_length
        pileup_end_pos      = fasta_flanking_length + flanking_length + (placeholder_length * len(tup.motif))
    
        for pileupcolumn in realigned_bamfile.pileup(STR_region[0], pileup_start_pos, pileup_end_pos, truncate=True, min_base_quality = 0, min_mapping_quality = mapq_threshold, max_depth=PILEUP_LIMIT): 

            for pileupread in pileupcolumn.pileups:

                read_name = pileupread.alignment.query_name

                # First iteration 
                if bool_firstIteration == False and pileupcolumn.pos == fasta_flanking_length - flanking_length:
            
                    if sc == True:
                        try:
                            CB = pileupread.alignment.get_tag("CB")
                        except: 
                            CB = None 
                            
                        try:
                            UMI = pileupread.alignment.get_tag("UB")
                        except: 
                            UMI = None 
                    
                    else:
                        CB = None 
                        UMI = None 
                        
                    dict_readname_to_STRinfo[read_name] = ["", flanking_length, flanking_length, pileupread.alignment.flag, CB, UMI,]
                
                if read_name in dict_readname_to_STRinfo.keys():
                    if not pileupread.is_refskip:

                        if pileupread.is_del:
                            dict_readname_to_STRinfo[read_name][0] += "*"
                            
                        elif pileupread.is_tail: # Reads that simply finished or have clippings within the region are filtered in this line
                            del dict_readname_to_STRinfo[read_name] 
                            
                        else: 
                            if pileupread.indel == 0:
                                dict_readname_to_STRinfo[read_name][0] += pileupread.alignment.query_sequence[pileupread.query_position]
                            elif pileupread.indel > 0: # insertion(s) ahead
                                
                                if pileupcolumn.pos == fasta_flanking_length - 1:
                                    dict_readname_to_STRinfo[read_name][0] += ( pileupread.alignment.query_sequence[pileupread.query_position] + pileupread.alignment.query_sequence[pileupread.query_position+1:pileupread.query_position+pileupread.indel+1])

                                elif pileupcolumn.pos < fasta_flanking_length - 1: # When insertion is before STR seq
                                    dict_readname_to_STRinfo[read_name][0] += ( pileupread.alignment.query_sequence[pileupread.query_position] + pileupread.alignment.query_sequence[pileupread.query_position+1:pileupread.query_position+pileupread.indel+1].lower())
                                    dict_readname_to_STRinfo[read_name][1] += pileupread.indel

                                elif pileupcolumn.pos > fasta_flanking_length + flanking_length - 1: # When insertion is after STR seq
                                    dict_readname_to_STRinfo[read_name][0] += ( pileupread.alignment.query_sequence[pileupread.query_position] + pileupread.alignment.query_sequence[pileupread.query_position+1:pileupread.query_position+pileupread.indel+1].lower())
                                    dict_readname_to_STRinfo[read_name][2] += pileupread.indel

                            elif pileupread.indel < 0: # Deletions
                                dict_readname_to_STRinfo[read_name][0] += pileupread.alignment.query_sequence[pileupread.query_position]  
      
    
        # (3) Store read information to pickle and write to temporary PATH
        alleleTable_entries = list()                    
        for readname, strInfo in dict_readname_to_STRinfo.items():
            
            left_flanking_seq   = strInfo[0][ :strInfo[1] ]
            str_seq             = strInfo[0][ strInfo[1] : -strInfo[2] ].replace("*", "")
            right_flanking_seq  = strInfo[0][ -strInfo[2]: ]
            reference_STR_allele = int( (STR_region[4] - STR_region[3] + 1) / len(STR_region[1]) )
            flag = strInfo[3]
            CB, UMI = strInfo[4], strInfo[5] 
            
            alleleTable_entries.append( [readname, f"{STR_region[0]}:{STR_region[3]}-{STR_region[4]}", STR_region[1], str_seq, reference_STR_allele, left_flanking_seq, right_flanking_seq, flag, CB, UMI] )
        
        if len(alleleTable_entries) != 0:
            nanomnt_utility.saveWithPickle(alleleTable_entries, PATH_temp, STR_name)

        # Delete realigned BAM file
        os.remove(f"{PATH_temp}/{STR_name}.sorted.bam")
        os.remove(f"{PATH_temp}/{STR_name}.sorted.bam.bai")

    return  

def runGetAlleleTable( PATH_bam, PATH_str_tsv, PATH_reference_genome, mapq_threshold, threads, flanking_length, realignment_f, sc, start_time, PATH_log, DIR_out ):
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    
    bam_filename = os.path.splitext(os.path.basename(PATH_bam))[0]
    
    # (1) Create folder for storing temporary files
    PATH_multiprocess_out = f"{DIR_out}/{bam_filename}.multiprocessing_temp"
    nanomnt_utility.checkAndCreate( DIR_out )
    nanomnt_utility.checkAndCreate( PATH_multiprocess_out )
        
    # (2) Prepare and chunk Krait SSR table    
    required_columns = ["sequence", "start", "end", "motif", "type", "repeat", "length"]
    # STR_table = load_SSR_table( PATH_str_tsv, required_columns, PATH_log )
    STR_table = pd.read_csv(PATH_str_tsv, sep='\t')
    STR_table = STR_table.sample(frac=1) # shuffle STR_table
        
    try:
        chunk_size = math.ceil( len(STR_table) / threads )
        STR_table_chunks = [STR_table[i:i+chunk_size] for i in range(0,STR_table.shape[0],chunk_size)]
    except ValueError: # Occurs when chunk size < threads
        logging.error(f"Number of input STR loci ({len(STR_table)}) is less than number of threads to use! Try reducing the number of threads to use.\t(elapsed time: {nanomnt_utility.getElapsedTime(start_time)} seconds)")
        raise ValueError
    
    # (4) Multiprocess using multiprocessing
    processes = list()
    logging.info(f"Starting multiprocessing with {threads}")
    
    # bamfile = pysam.AlignmentFile(PATH_bam, "rb")
    
    for idx, STR_chunk in enumerate(STR_table_chunks):
        PATH_temp = f"{PATH_multiprocess_out}/thread_{idx+1}"
        nanomnt_utility.checkAndCreate(PATH_temp)
        
        p = multiprocessing.Process(target=extractSTR_multiprocess,
                                    args=[PATH_bam, STR_chunk, PATH_reference_genome, flanking_length, realignment_f, mapq_threshold, sc, PATH_temp] )
        p.start()
        processes.append(p)
    for process in processes:
        process.join()
    
    # (5) Collect multiprocessing results  
    alleleTable_entries = list()
    for i in range(len( STR_table_chunks )):
        idx = i + 1
        
        list_PATH_pickles = glob.glob(f"{PATH_multiprocess_out}/thread_{idx}/*.pickle")
        
        for PATH_pickle in list_PATH_pickles:
            cur_alleleTable_entries = nanomnt_utility.loadFromPickle(PATH_pickle)
            for entry in cur_alleleTable_entries:
                alleleTable_entries.append( entry )
            
            # Delete pickle files
            os.remove(PATH_pickle)
        
        os.rmdir( f"{PATH_multiprocess_out}/thread_{idx}" )

    df_alleleTable = pd.DataFrame(alleleTable_entries, columns=["read_name", "locus", "repeat_unit", "allele", "reference_STR_allele", "left_flanking_seq", "right_flanking_seq", "flag", "CB", "UMI"])
    logging.info(f"Total of {len(set(df_alleleTable.locus))} microsatellite loci detected.\t(elapsed time: {nanomnt_utility.getElapsedTime(start_time)} seconds)")

    # (6) Prepare error correction 
    try:
        chunk_size = math.ceil( len(df_alleleTable) / threads )
        STR_allele_table_chunks = [df_alleleTable[i:i+chunk_size] for i in range(0,df_alleleTable.shape[0],chunk_size)]
    except ValueError: # Occurs when chunk size < threads
        logging.error(f"Number of reads ({len(df_alleleTable)}) is less than number of threads to use! Multiprocessing won't be used for error-correction.")
    
    # (7) Start error correction
    processes = list()
    logging.info(f"Starting multiprocessing with {threads}")    
    for idx, STR_allele_table_chunk in enumerate(STR_allele_table_chunks):
        
        p = multiprocessing.Process(target=errorCorrection_multiprocess,
                                    args=[STR_allele_table_chunk, None, idx+1, PATH_multiprocess_out] )
        p.start()
        processes.append(p)
    for process in processes:
        process.join()
    
    # (8) Collect multiprocessing results & save results to disk   
    list_PATH_allele_tables = glob.glob( f"{PATH_multiprocess_out}/corrected_allele_table.*.tsv" )
    df_concat_allele_table = pd.concat( [pd.read_csv(PATH_allele_table, sep='\t') for PATH_allele_table in list_PATH_allele_tables] )
    num_failed_reads    = df_concat_allele_table[(df_concat_allele_table["editing_distance"]==-1)].shape[0]
    correction_rate     = 100 * ( 1 - ( num_failed_reads / len(df_concat_allele_table) ) )
    logging.info(f"Finished error-correction. (Correction rate: {round(correction_rate, 2)}%) (elapsed time: {nanomnt_utility.getElapsedTime(start_time)} seconds)")

    bam_filename = os.path.splitext(os.path.basename(PATH_bam))[0]
    logging.info(f"Writing STR allele table to disk")
    df_concat_allele_table.to_csv(f"{DIR_out}/{bam_filename}.STR_allele_table.tsv", sep='\t', index=False)
    

    for PATH_allele_table in list_PATH_allele_tables:
        os.remove( PATH_allele_table )
    
    # list_temp_PATH = glob.glob(f"{PATH_multiprocess_out}/*")
    # deleteAllPATHs = False
    # for temp_PATH in list_temp_PATH:
    #     if len(os.listdir(temp_PATH)) == 0:
    #         os.rmdir(temp_PATH)
    #         deleteAllPATHs = True
    # if deleteAllPATHs == True:
    #     os.rmdir(PATH_multiprocess_out)   
    os.rmdir(PATH_multiprocess_out)
    
    return


def main():
    start_time  = time.time()

    github_link = "https://github.com/18parkky/NanoMnT"
    script_description = f"Given a BAM file, generate STR allele table. For more information, see GitHub: {github_link}"
    parser = argparse.ArgumentParser(description=script_description)
    
    # script_PATH     = os.path.dirname(os.path.dirname(__file__))

    # reference_PATH  = os.path.join( script_PATH, "ref/genome" )
    # krait_PATH      = os.path.join( script_PATH, "ref/krait" )

    # PATH_default_krait           = f"{krait_PATH}/T2TCHM13_SSR.tsv.gz"

    # Read paramters
    # Required parameters
    parser.add_argument('-b', '--PATH_bam',      
                        help="PATH to input BAM file (must be sorted and indexed)", 
                        required=True,
                        )
    
    parser.add_argument('-s', '--PATH_str_tsv',  
                        help="PATH to STR list file generated using either Krait or Pytrf (.tsv)", 
                        required=True,
                        )    
    
    parser.add_argument('-r', 
                        "--PATH_ref_genome", help="PATH to reference genome (must be same as the one used for generating BAM)",
                        required=True,
                        )
    
    # Optional parameters
    parser.add_argument('-m', '--mapping_quality', 
                        help="Mapping quality threshold when processing reads (default: 60)", 
                        required=False, 
                        type=int, 
                        default=60,
                        )
    
    parser.add_argument('-f', '--flanking',     
                        help="Length of flanking sequences to collect for each read (default: 12)", 
                        required=False, 
                        type=int,
                        default=12,
                        )
    
    parser.add_argument('-rf', '--realignment_flanking',     
                        help="Length of flanking sequences of STR to use as reference, during re-alignment process (default: 20000)", 
                        required=False, 
                        type=int,
                        default=20000,
                        )
    
    parser.add_argument('-t', '--threads',      
                        help="Number of threads (default: 4)",
                        required=False, 
                        type=int, 
                        default=4,
                        )
    
    parser.add_argument('--sc',                 
                        help="Input is single-cell data. Each read in BAM must contain CB and UB tag. (default: False)", 
                        action='store_true',
                        )    
    
    parser.add_argument('-out', '--DIR_out',  
                        help='Directory to write output files (default: current directory)', 
                        required=False, 
                        type=str, 
                        default=os.getcwd(),
                        )

    args = vars(parser.parse_args())

    PATH_bam        = args["PATH_bam"]
    PATH_str_tsv    = args["PATH_str_tsv"]
    PATH_ref_genome = args["PATH_ref_genome"]
    mapq_threshold  = args["mapping_quality"]
    num_threads     = args["threads"]
    flanking_length = args["flanking"]
    realignment_f   = args["realignment_flanking"]
    sc              = args["sc"]
    DIR_out         = args["DIR_out"]
    
    bam_filename = os.path.splitext(os.path.basename(PATH_bam))[0]
    
    # Create log file
    nanomnt_utility.checkAndCreate(DIR_out)
    PATH_log = f"{DIR_out}/nanomnt.getAlleleTable.{bam_filename}.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'{k}\t:\t{v}')
        
    runGetAlleleTable(PATH_bam, PATH_str_tsv, PATH_ref_genome, mapq_threshold, num_threads, flanking_length, realignment_f, sc, start_time, PATH_log, DIR_out)
    logging.info(f"Finished getAlleleTable.py\t(Total time taken: {nanomnt_utility.getElapsedTime(start_time)} seconds)")

        
if __name__ == "__main__":
    main()