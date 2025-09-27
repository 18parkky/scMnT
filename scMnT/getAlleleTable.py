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
import Levenshtein
import pandas as pd
import multiprocessing

import scmnt_utility

PILEUP_LIMIT = 5000000 # pysam's pileup's default value for max_depth is 8000

def create_STR_loci_batches(STR_table_chunk, adjacency_threshold):

    batches = [ list() ]
    for sequence, edf in STR_table_chunk.groupby('sequence', observed=True):
        for tup in edf.sort_values('start').itertuples():
            
            # Assign locus to batch
            for batch_idx in range(len(batches)):
                if len(batches[batch_idx])==0:
                    batches[batch_idx].append( tup.locus ) 
                    break
                else:
                    locusA_end      = int(batches[batch_idx][-1].split('_')[3])
                    locusB_start    = int(tup.start)
                    distance = locusB_start - locusA_end
                    if distance < adjacency_threshold: # If loci are adjacent, current locus needs to be assigned to another batch
                        if batch_idx == len(batches)-1:
                            batches.append( [tup.locus] )   # Create new batch if out of batch
                            break
                        else:
                            continue
                    else:
                        batches[batch_idx].append( tup.locus )
                        break
                    
    STR_table_chunk_batches = list()
    for batch in batches:
        STR_table_chunk_batches.append( STR_table_chunk[(STR_table_chunk['locus'].isin( batch ))] )
    
    return STR_table_chunk_batches

def identify_poly_loci_reads(PATH_bam, STR_table, flanking_length, mapq_threshold):
    bamfile = pysam.AlignmentFile(PATH_bam, 'rb')
    dict_readname_to_loci = dict()
    for tup in STR_table.itertuples():
        chrom = tup.sequence 
        start, end = tup.start, tup.end
        for read in bamfile.fetch(chrom, start - 1 - flanking_length, end + flanking_length):
            if read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue
            if read.mapping_quality < mapq_threshold:
                continue
            try:
                dict_readname_to_loci[read.query_name].append( tup.locus )
            except KeyError:
                dict_readname_to_loci[read.query_name] = [ tup.locus ]
    
    poly_loci_reads = set()
    dict_locus_to_polyLociReads = dict()
    for readname, loci in dict_readname_to_loci.items():
        if len(loci)==1: continue 
        
        poly_loci_reads.add(readname)
        for locus in loci:
            try: dict_locus_to_polyLociReads[locus].append(readname)
            except KeyError: dict_locus_to_polyLociReads[locus] = [ readname ]
        
    return poly_loci_reads, dict_locus_to_polyLociReads

def chunk_genes_by_counts(STR_table, dict_locus_to_readcount, threads):
    # Initialize partitions
    loci_chunks = [dict() for i in range(threads)]
    chunk_sums  = [0 for i in range(threads)]
    
    # Sort genes by read count
    dict_locus_to_readcount_sorted = sorted(dict_locus_to_readcount.items(), key=lambda x: x[1], reverse=True)
    
    # Greedy assignment: put each gene in the partition with the smallest sum
    for locus, readcount in dict_locus_to_readcount_sorted:
        idx = chunk_sums.index(min(chunk_sums)) # Find index of partition with the lowest sum so far
        loci_chunks[idx][locus] = readcount
        chunk_sums[idx] += readcount
        
    chunks = list()
    for loci_chunk in loci_chunks:
        chunks.append( STR_table[(STR_table['locus'].isin( loci_chunk ))] )
    
    return chunks

def correctAllele(STR_allele_table, dist_threshold=None):
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

def extract_STR_reads(STR_table, poly_loci_reads, PATH_temp, PATH_bamfile, flanking_length, mapq_threshold):
    bamfile = pysam.AlignmentFile(PATH_bamfile, "rb")
    PATH_MS_reads_fasta = f'{PATH_temp}/reads.fasta'
    
    dict_readnames_to_CB_UMI = dict()
    
    with open(PATH_MS_reads_fasta, 'w') as fasta:
    
        for tup in STR_table.itertuples():
            chrom, start, end = str(tup.sequence), tup.start, tup.end
            
            for read in bamfile.fetch(chrom, start - 1 - flanking_length, end + flanking_length):
                if read.is_secondary or read.is_supplementary or read.is_duplicate:
                    continue
                if read.mapping_quality < mapq_threshold:
                    continue

                qname = read.query_name
                CB = read.get_tag("CB") if read.has_tag("CB") else None
                UMI = read.get_tag("UB") if read.has_tag("UB") else None

                dict_readnames_to_CB_UMI[qname] = [CB, UMI]
                
                if qname in poly_loci_reads: continue 
                fasta.write(f'>{qname}\n{read.get_forward_sequence()}\n')

    return dict_readnames_to_CB_UMI

def create_STR_reference(STR_table, realignment_flanking_length, PATH_reference_genome, placeholder_length, PATH_temp):
    
    PATH_fasta = f"{PATH_temp}/ref.fasta"
    pysam_genome = pysam.FastaFile(PATH_reference_genome)
    
    with open(PATH_fasta, "w") as fasta: 
        for tup in STR_table.itertuples():
            chrom = tup.sequence 
            motif = tup.motif 
            repeat = tup.repeat 
            start, end = tup.start, tup.end
            
            STR_name = f'{chrom}_{motif}x{repeat}_{start}_{end}'
            lf, rf  = getFlankingSequence( pysam_genome, chrom, start, end, realignment_flanking_length )
            placeholder_sequence = motif * placeholder_length
            fasta.write(f">{STR_name}\n")
            fasta.write(f"{lf}{placeholder_sequence}{rf}\n")
    return

def tag_BAM(dict_readnames_to_CB_UMI, PATH_realigned_bam_out, PATH_realigned_tagged_bam_out,):
    realigned_bamfile = pysam.AlignmentFile(PATH_realigned_bam_out, "rb")
    
    realigned_tagged_bamfile = pysam.AlignmentFile(PATH_realigned_tagged_bam_out, "wb", template=realigned_bamfile)
    for read in realigned_bamfile.fetch():
        read_name = read.query_name
        read.set_tag( "CB", dict_readnames_to_CB_UMI[read_name][0] )
        read.set_tag( "UB", dict_readnames_to_CB_UMI[read_name][1] )
        realigned_tagged_bamfile.write(read)
        
    realigned_tagged_bamfile.close()

    # Delete realigned BAM
    os.remove( PATH_realigned_bam_out )
    os.rename( PATH_realigned_tagged_bam_out, PATH_realigned_bam_out)
    pysam.index(PATH_realigned_bam_out)
    return

def collect_realignment_info(PATH_realigned_bam_out, placeholder_length, realignment_flanking_length, flanking_length, mapq_threshold):
    dict_readname_to_STRinfo = dict()
    realigned_bamfile = pysam.AlignmentFile(PATH_realigned_bam_out, 'rb')

    for STR_locus in realigned_bamfile.references:
        
        chrom, motif_repeat, start, end = STR_locus.split('_')
        motif, repeat = motif_repeat.split('x')
        
        bool_firstIteration = False
        
        pileup_start_pos    = realignment_flanking_length - flanking_length
        pileup_end_pos      = realignment_flanking_length + flanking_length + (placeholder_length * len(motif))

        for pileupcolumn in realigned_bamfile.pileup(STR_locus, pileup_start_pos, pileup_end_pos, truncate=True, min_base_quality = 0, min_mapping_quality = mapq_threshold, max_depth=PILEUP_LIMIT): 
            for pileupread in pileupcolumn.pileups:
                
                read_name = pileupread.alignment.query_name
                
                # First iteration 
                if bool_firstIteration == False and pileupcolumn.pos == realignment_flanking_length - flanking_length:
            
                    try:
                        CB = pileupread.alignment.get_tag("CB")
                    except: 
                        CB = None 
                        
                    try:
                        UMI = pileupread.alignment.get_tag("UB")
                    except: 
                        UMI = None 

                        
                    dict_readname_to_STRinfo[read_name] = ["", flanking_length, flanking_length, pileupread.alignment.flag, CB, UMI, STR_locus]
                
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
                                
                                if pileupcolumn.pos == realignment_flanking_length - 1:
                                    dict_readname_to_STRinfo[read_name][0] += ( pileupread.alignment.query_sequence[pileupread.query_position] + pileupread.alignment.query_sequence[pileupread.query_position+1:pileupread.query_position+pileupread.indel+1])

                                elif pileupcolumn.pos < realignment_flanking_length - 1: # When insertion is before STR seq
                                    dict_readname_to_STRinfo[read_name][0] += ( pileupread.alignment.query_sequence[pileupread.query_position] + pileupread.alignment.query_sequence[pileupread.query_position+1:pileupread.query_position+pileupread.indel+1].lower())
                                    dict_readname_to_STRinfo[read_name][1] += pileupread.indel

                                elif pileupcolumn.pos > realignment_flanking_length + flanking_length - 1: # When insertion is after STR seq
                                    dict_readname_to_STRinfo[read_name][0] += ( pileupread.alignment.query_sequence[pileupread.query_position] + pileupread.alignment.query_sequence[pileupread.query_position+1:pileupread.query_position+pileupread.indel+1].lower())
                                    dict_readname_to_STRinfo[read_name][2] += pileupread.indel

                            elif pileupread.indel < 0: # Deletions
                                dict_readname_to_STRinfo[read_name][0] += pileupread.alignment.query_sequence[pileupread.query_position]
    return dict_readname_to_STRinfo

def summarize_STR_info(dict_readname_to_STRinfo, batch_idx, PATH_temp):
    AlleleTable_i = list()
    for readname, STRinfo in dict_readname_to_STRinfo.items():
        left_flanking_sequence  = STRinfo[0][ :STRinfo[1] ]
        right_flanking_sequence = STRinfo[0][ -STRinfo[2]: ]
        STR_sequence            = STRinfo[0][ STRinfo[1] : -STRinfo[2] ].replace("*", "")
        
        locus_parts = STRinfo[6].split('_')
        chrom, start, end = str(locus_parts[0]), int(locus_parts[2]), int(locus_parts[3])
        motif = locus_parts[1].split('x')[0]
        
        flag                    = STRinfo[3]
        CB, UMI                 = STRinfo[4], STRinfo[5] 
        reference_STR_allele = int(locus_parts[1].split('x')[1])
        if CB==None or UMI==None: continue
        
        AlleleTable_i.append( [readname, f"{chrom}:{start}-{end}", motif, STR_sequence, reference_STR_allele, left_flanking_sequence, right_flanking_sequence, flag, CB, UMI] )
        
    if len(AlleleTable_i)!=0:
        scmnt_utility.saveWithPickle(AlleleTable_i, PATH_temp, f'AlleleTable_temp_batch{batch_idx}')
        
    # Delete realigned BAM file
    os.remove(f"{PATH_temp}/realigned.sorted.bam")
    os.remove(f"{PATH_temp}/realigned.sorted.bam.bai")
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

def extractSTR_multiprocess( PATH_bamfile, STR_table_chunk, PATH_reference_genome, flanking_length, realignment_flanking_length, mapq_threshold, PATH_temp ):
    
    STR_table_chunk_batches = create_STR_loci_batches(STR_table_chunk, realignment_flanking_length * 2)
    
    for batch_idx, STR_table in enumerate(STR_table_chunk_batches):
        
        # (0) Identify reads that map to a single or multiple STR loci
        poly_loci_reads, dict_locus_to_polyLociReads = identify_poly_loci_reads(PATH_bamfile, STR_table, flanking_length, mapq_threshold)
        
        # (1) Extract reads that map to STR loci
        dict_readnames_to_CB_UMI = extract_STR_reads(STR_table, poly_loci_reads, PATH_temp, PATH_bamfile, flanking_length, mapq_threshold)

        # (2) Create custom reference FASTA for STR loci
        placeholder_length = 50
        create_STR_reference(STR_table, realignment_flanking_length, PATH_reference_genome, placeholder_length, PATH_temp)
        
        # (3) Minimap2
        PATH_realigned_bam_out = f"{PATH_temp}/realigned.sorted.bam"
        command = f"minimap2 -a -A 4 -B 10 -t 1 {PATH_temp}/ref.fasta {PATH_temp}/reads.fasta 2>/dev/null | samtools view -Sb -F0x900 | samtools sort > {PATH_realigned_bam_out}"   # Filter out supplementary & secondary
        subprocess.call(command, shell=True)
        pysam.index(PATH_realigned_bam_out)
        
        # (4) Set CB and UMI tags to realigned BAM
        tag_BAM(dict_readnames_to_CB_UMI, PATH_realigned_bam_out, f'{PATH_temp}/realigned.sorted.tagged.bam',)
            
        # (5) Delete FASTA files since we now have the realigned BAM
        tempFiles = glob.glob(f"{PATH_temp}/*.fasta")
        for file in tempFiles:
            os.remove( file )
            
        # (6) Collect STR information for each re-aligned read
        # key: read_name, value: [STR_seq, base quality, int_left_padding, int_right_padding, flag]
        dict_readname_to_STRinfo = collect_realignment_info(PATH_realigned_bam_out, placeholder_length, realignment_flanking_length, flanking_length, mapq_threshold,)
        
        # (7) Process poly locus reads 
        bamfile = pysam.AlignmentFile(PATH_bamfile, 'rb')
        pysam_genome = pysam.FastaFile(PATH_reference_genome)

        for tup in STR_table[(STR_table['locus'].isin( list(dict_locus_to_polyLociReads.keys()) ))].itertuples():
            chrom, start, end, motif = tup.sequence, tup.start, tup.end, tup.motif

            read_fasta = open(f'{PATH_temp}/{tup.locus}_reads.fasta', 'w')
            for read in bamfile.fetch(chrom, start - 1 - flanking_length, end + flanking_length):
                if read.is_secondary or read.is_supplementary or read.is_duplicate:
                    continue
                if read.mapping_quality < mapq_threshold:
                    continue
                if read.query_name not in poly_loci_reads: continue 
                read_fasta.write(f'>{read.query_name}\n{read.get_forward_sequence()}\n')
            read_fasta.close()
            
            ref_fasta = open(f'{PATH_temp}/{tup.locus}_ref.fasta', 'w')
            lf, rf  = getFlankingSequence( pysam_genome, chrom, start, end, realignment_flanking_length )
            placeholder_sequence = motif * placeholder_length
            ref_fasta.write(f">{tup.locus}\n")
            ref_fasta.write(f"{lf}{placeholder_sequence}{rf}\n")
            ref_fasta.close()
            
            PATH_s_realigned_bam_out = f'{PATH_temp}/{tup.locus}_realigned.sorted.bam'
            command = f"minimap2 -a -A 4 -B 10 -t 1 {PATH_temp}/{tup.locus}_ref.fasta {PATH_temp}/{tup.locus}_reads.fasta 2>/dev/null | samtools view -Sb -F0x900 | samtools sort > {PATH_s_realigned_bam_out}"   # Filter out supplementary & secondary
            subprocess.call(command, shell=True)
            pysam.index(PATH_s_realigned_bam_out)
            
            PATH_s_realigned_tagged_bam_out = f'{PATH_temp}/{tup.locus}.sorted.tagged.bam'
            try:
                tag_BAM(dict_readnames_to_CB_UMI, PATH_s_realigned_bam_out, PATH_s_realigned_tagged_bam_out)
            except ValueError: 
                continue 
            
            dict_s_readname_to_STRinfo = collect_realignment_info(PATH_s_realigned_bam_out, placeholder_length, realignment_flanking_length, flanking_length, mapq_threshold)
                
            for k,v in dict_s_readname_to_STRinfo.items(): 
                dict_readname_to_STRinfo[k]=v 
                
            os.remove(f'{PATH_temp}/{tup.locus}_ref.fasta')
            os.remove(f'{PATH_temp}/{tup.locus}_reads.fasta')
            os.remove(PATH_s_realigned_bam_out)
            os.remove(f'{PATH_s_realigned_bam_out}.bai')

        # (8) Summarize realignment results
        summarize_STR_info(dict_readname_to_STRinfo, batch_idx, PATH_temp)

def cleanup_multiprocessing_results(PATH_multiprocess_out, threads, start_time):
    # (4) Collect multiprocessing results  
    AlleleTable = list()
    for i in range(threads):
        idx = i + 1
        list_PATH_pickles = glob.glob(f"{PATH_multiprocess_out}/thread_{idx}/*.pickle")
        
        for PATH_pickle in list_PATH_pickles:
            cur_alleleTable_entries = scmnt_utility.loadFromPickle(PATH_pickle)
            for entry in cur_alleleTable_entries:
                AlleleTable.append( entry )
            
            # Delete pickle files
            os.remove(PATH_pickle)
        
        os.rmdir( f"{PATH_multiprocess_out}/thread_{idx}" )
        
    AlleleTable = pd.DataFrame(AlleleTable, columns=["read_name", "locus", "repeat_unit", "allele", "reference_STR_allele", "left_flanking_seq", "right_flanking_seq", "flag", "CB", "UMI"])
    logging.info(f"Total of {len(set(AlleleTable['locus']))} microsatellite loci detected.\t(elapsed time: {scmnt_utility.getElapsedTime(start_time)} seconds)")
    
    return AlleleTable 

def prepare_STR_chunks( PATH_str_tsv, PATH_bamfile, mapq_threshold, threads, start_time ):
    STR_table = pd.read_csv(PATH_str_tsv, sep='\t')
    STR_table['locus'] = [ f'{tup.sequence}_{tup.motif}x{tup.repeat}_{tup.start}_{tup.end}' for tup in STR_table.itertuples() ]
    # STR_table = STR_table.sample(frac=1) # shuffle STR_table
    loci_count = len(STR_table)
        
    # Chunk based on gene expressions
    logging.info(f"Estimating read counts for each MS locus. (elapsed time: {scmnt_utility.getElapsedTime(start_time)} seconds)")
    bamfile = pysam.AlignmentFile(PATH_bamfile, 'rb')
    dict_locus_to_readcount = dict()
    for tup in STR_table.itertuples():
        chrom = tup.sequence 
        start, end = tup.start, tup.end
        dict_locus_to_readcount[tup.locus] = 0
        for read in bamfile.fetch(chrom, start - 1, end):
            if read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue
            if read.mapping_quality < mapq_threshold:
                continue 
            dict_locus_to_readcount[tup.locus] += 1
            
    dict_locus_to_nonzero_readcount = { k:v for k,v in dict_locus_to_readcount.items() if v>0 }
    
    STR_table = STR_table[(STR_table['locus'].isin( dict_locus_to_nonzero_readcount.keys() ))].copy()
    
    logging.info(f"Processing {len(dict_locus_to_nonzero_readcount)} loci instead of {loci_count} loci (elapsed time: {scmnt_utility.getElapsedTime(start_time)} seconds)")

    STR_table_chunks = chunk_genes_by_counts(STR_table, dict_locus_to_nonzero_readcount, threads)
    return STR_table_chunks

def runGetAlleleTable( PATH_bam, PATH_str_tsv, PATH_reference_genome, mapq_threshold, threads, flanking_length, realignment_flanking_length, start_time, filename, DIR_out ):
        
    # (1) Create folder for storing temporary files
    PATH_multiprocess_out = f"{DIR_out}/{filename}.multiprocessing_temp"
    scmnt_utility.checkAndCreate( DIR_out )
    scmnt_utility.checkAndCreate( PATH_multiprocess_out )
        
    # (2) Prepare and chunk Krait SSR table
    STR_table_chunks = prepare_STR_chunks( PATH_str_tsv, PATH_bam, mapq_threshold, threads, start_time )
    
    # (3) Multiprocessing
    processes = list()
    logging.info(f"Starting multiprocessing with {threads}")    
    for idx, STR_chunk in enumerate(STR_table_chunks):
        PATH_temp = f"{PATH_multiprocess_out}/thread_{idx+1}"
        scmnt_utility.checkAndCreate(PATH_temp)
        
        p = multiprocessing.Process(target=extractSTR_multiprocess,
                                    args=[PATH_bam, STR_chunk, PATH_reference_genome, flanking_length, realignment_flanking_length, mapq_threshold, PATH_temp] )
        p.start()
        processes.append(p)
    for process in processes:
        process.join()
    
    # (4) Collect multiprocessing results  
    AlleleTable = cleanup_multiprocessing_results(PATH_multiprocess_out, threads, start_time)
    
    # (6) Prepare error correction 
    try:
        chunk_size = math.ceil( len(AlleleTable) / threads )
        STR_allele_table_chunks = [AlleleTable[i:i+chunk_size] for i in range(0,AlleleTable.shape[0],chunk_size)]
    except ValueError: # Occurs when chunk size < threads
        logging.error(f"Number of reads ({len(AlleleTable)}) is less than number of threads to use! Multiprocessing won't be used for error-correction.")
    
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
    logging.info(f"Finished error-correction. (Correction rate: {round(correction_rate, 2)}%) (elapsed time: {scmnt_utility.getElapsedTime(start_time)} seconds)")

    logging.info(f"Writing STR allele table to disk")
    df_concat_allele_table.to_csv(f"{DIR_out}/{filename}.AlleleTable.tsv", sep='\t', index=False)

    for PATH_allele_table in list_PATH_allele_tables:
        os.remove( PATH_allele_table )
        
    os.rmdir(PATH_multiprocess_out)
    return

def main():
    start_time  = time.time()

    github_link = "https://github.com/18parkky/scMnT"
    script_description = f"[scMnT] Given a BAM file, generate STR allele table. For more information, see GitHub: {github_link}"
    parser = argparse.ArgumentParser(description=script_description)
    
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    # Required parameters
    required.add_argument('-b', '--PATH_bam',      
                        help="PATH to input BAM file (must be sorted and indexed)", 
                        required=True,
                        )
    
    required.add_argument('-s', '--PATH_str_tsv',  
                        help="PATH to STR list file generated using either Krait or Pytrf (.tsv)", 
                        required=True,
                        )    
    
    required.add_argument('-r', '--PATH_reference_genome', 
                        help="PATH to reference genome (must be same as the one used for generating BAM)",
                        required=True,
                        )
    
    required.add_argument('-n', '--filename',      
                        help="Name of the output files", 
                        required=True,
                        )
    
    # Optional parameters
    optional.add_argument('-f', '--flanking',     
                        help="Length of flanking sequences to collect for each read (default: 6)", 
                        required=False, 
                        type=int,
                        default=6,
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

    PATH_bam                = args['PATH_bam']
    PATH_str_tsv            = args['PATH_str_tsv']
    PATH_reference_genome   = args['PATH_reference_genome']
    filename                = args['filename']
    threads                 = args['threads']
    flanking_length         = args['flanking']
    DIR_out                 = args['DIR_out']
        
    # Create log file
    scmnt_utility.checkAndCreate(DIR_out)
    PATH_log = f"{DIR_out}/scmnt.getAlleleTable.{filename}.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'{k}\t:\t{v}')
        
    mapq_threshold              = 60
    realignment_flanking_length = 50
    
    runGetAlleleTable( PATH_bam, PATH_str_tsv, PATH_reference_genome, mapq_threshold, threads, flanking_length, realignment_flanking_length, start_time, filename, DIR_out )
    logging.info(f"Finished getAlleleTable.py\t(Total time taken: {scmnt_utility.getElapsedTime(start_time)} seconds)")

if __name__ == "__main__":
    main()