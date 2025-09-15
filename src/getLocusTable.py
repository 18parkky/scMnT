import os 
import gc
import time
import glob
import math
import logging
import argparse
import numpy as np
import pandas as pd
import multiprocessing

from scipy.signal import find_peaks
from scipy.optimize import fsolve

import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("darkgrid", {"grid.color": ".6", "grid.linestyle": ":"})

import NanoMnT.nanomnt_utility as nanomnt_utility

dict_repeatUnit_to_prominence = {
    "A" : (0.2, 0.4),
    "T" : (0.2, 0.4),
}

min_numEA, max_numEA = 4.96442686645189, 8.763010850981194

def getFunc_numEA_to_prominence(  numEA, numEA_max, numEA_min, P_max, P_min ):
    slope       = ( P_max - P_min ) / ( numEA_min - numEA_max )
    y_intercept = P_max - slope * numEA_min
    return slope * numEA + y_intercept

def generateHistogram( x_peak, y_peak, histogram_length, prominence ):    

    # Define the function based on the equation
    def getR(x, n, a, S):
        return 2*a*x**n - S*x + (S - 2*a)

    def geometricProgression( a, r, n ):
        return a * r**(n-1)
    
    if ( 0 < prominence <= 1 ) is False:
        raise ValueError

    H = y_peak                      # Peak y-value
    S = 1 - H                       # Sum of all y-values excluding peak y-value
    a = H * (1-prominence)          # peak neighboring dot yvalue
    n = histogram_length            # histogram length

    # Use fsolve to find the root
    r = fsolve(getR, x0=0, args=(n, a, S))[0]
    # print(f"r:\t{r}")
    
    dict_allele_histogram = dict()

    for i in range(1, n + 1):
        an = geometricProgression( a , r, i )
        dict_allele_histogram[i + x_peak] = an * 100
        dict_allele_histogram[-i + x_peak] = an * 100

    dict_allele_histogram[0 + x_peak] = H * 100
    
    return dict_allele_histogram

def findBestPeak( find_peaks_out, ref_length, filling_range ):
    c = ref_length - filling_range 
    peaks = find_peaks_out[0]
    prominences = find_peaks_out[1]["prominences"]
    
    dict_peak_to_prominence = { peak : prominences[idx] for idx, peak in enumerate(peaks) }
    dict_peak_to_prominence = dict(sorted(dict_peak_to_prominence.items(), key=lambda x:x[1], reverse=True))
    
    for k, v in dict_peak_to_prominence.items():
        return k + c, v
    
def process_STR_locus(STR_allele_table, batch_id, read_selection, 
                    min_coverage, filename, PATH_log, get_allele_size_info,
                    DIR_multiprocess_out, DIR_alleleSizeDist_out):
    
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
        
    # Iterate through each STR locus and generate locus table
    list_of_lociInfos = [ ]
    for locus, edf_at in STR_allele_table.groupby('locus'):
        
        chromosome  = str( locus.split(":")[0] )
        start, end  = int( locus.split(":")[1].split("-")[0] ), int( locus.split(":")[1].split("-")[1] )
        
        repeat_unit = str( edf_at.iloc[0].repeat_unit )
        reference_STR_allele = int( ( int( locus.split(':')[1].split('-')[1] ) - int( locus.split(':')[1].split('-')[0] ) + 1 ) / len(repeat_unit) )
        
        if repeat_unit == "A" or repeat_unit == "T":
            # (2-1) Select edf based on read selection approach
        
            if read_selection == "fw_rv_reads":
                df_at_of_choice = edf_at[(edf_at["flag"]==nanomnt_utility.dict_repeatUnit_to_favorableFlag[repeat_unit])]
            if read_selection == "all_reads":
                df_at_of_choice = edf_at 
        else:
    
            df_at_of_choice = edf_at
            
        df_at_of_choice = df_at_of_choice[~(df_at_of_choice["corrected_allele"].isna()) & (df_at_of_choice["editing_distance"]!=-1)]
        # Skip locus if coverage is insufficient
        if len(df_at_of_choice) <= min_coverage:
            continue

        # (2-2) Calculate percentage of reference allele
        allele_count = pd.Series( dict( (e, c) for e, c in zip( * np.unique( [ ca.upper() for ca in df_at_of_choice.corrected_allele ], return_counts = True ) ) ) )
        try: perc_ref_allele = allele_count[repeat_unit * reference_STR_allele] / allele_count.sum( )
        except KeyError: perc_ref_allele = None

        # (2-4) Genotype locus and create allele size distribution information
        observed_allele_histogram = nanomnt_utility.getAlleleHistogram( df_at_of_choice, nanomnt_utility.allele_histogram_width )
        observed_allele_histogram = dict(sorted(observed_allele_histogram.items(), key=lambda x:x[0]))
        
            # (2-4-1) Main-approach: Create synthetic histograms and find best match
        allele_freq = allele_count / allele_count.sum( )    
        numEA = 1 / np.sum( allele_freq ** 2 )
        
        
        # list_allele_counts = list( observed_allele_histogram.values() )
            # Generate synthetic histograms
        dict_length_to_alleleHistogram = dict()
        allele_search_range = range( 1, 100 + 1 )   # STR whose size exceeds 100 may not be processed properly
        for i in allele_search_range:
            try:
                min_prominence = dict_repeatUnit_to_prominence[repeat_unit][0]
                max_prominence = dict_repeatUnit_to_prominence[repeat_unit][1]
            except:
                min_prominence = 0.4
                max_prominence = 0.6
            
            allele_histogram = generateHistogram( i, max(observed_allele_histogram.values())/100, 10, getFunc_numEA_to_prominence(numEA, max_numEA, min_numEA, max_prominence, min_prominence) ) # (0.4 ,0.2) works very well for A-repeats
            dict_length_to_alleleHistogram[i] = dict(sorted(allele_histogram.items(), key=lambda x : x[0]))
            # Calculate distance for each synthetic histogram
        dict_allele_to_histogramDist = dict()
        for putative_allele in allele_search_range:
            putative_allele_histogram = dict_length_to_alleleHistogram[putative_allele]
            dist = nanomnt_utility.calcHistogramDistance( observed_allele_histogram, putative_allele_histogram, "cosine" )
            dict_allele_to_histogramDist[putative_allele] = dist
            # Get the histogram with smallest distance (closest histogram)
        dict_allele_to_histogramDist = dict(sorted(dict_allele_to_histogramDist.items(), key=lambda x : x[1]))
        for allele in dict_allele_to_histogramDist.keys():
            genotyped_allele = allele
            break 
        # genotyped_allele = 0

            # (2-4-2) Sub-approach: Find most prominent peak using scipy's find_peak()
        peaks = find_peaks( list(observed_allele_histogram.values()), prominence=1 )
        try:
            majorPeak, prominence = findBestPeak(peaks, reference_STR_allele, nanomnt_utility.allele_histogram_width)
        except TypeError: # No peaks have been found
            majorPeak, prominence = None, None

        if get_allele_size_info == True:

            DIR_alleleSizeDist_subdirectory = f"{DIR_alleleSizeDist_out}/{repeat_unit}/{chromosome}"

            # (File 1) Allele size distribution histogram
            coverage = len( df_at_of_choice[(df_at_of_choice["corrected_allele"] != "uncor")] )

            f = sns.lineplot( observed_allele_histogram )
            f = sns.scatterplot( observed_allele_histogram )

            plt.axvline( x=int(genotyped_allele), color="dodgerblue", linestyle=':', linewidth=2, alpha=0.7 )
            plt.axvline( x=int(reference_STR_allele), color="red", linestyle=':', linewidth=2, alpha=0.7 )
            # plt.axvline( x=int(majorPeak), color="gold", linestyle=':', linewidth=2, alpha=0.8 )
            
            f.set(title=f"{repeat_unit}x{genotyped_allele}_cov:{coverage}")
            f.set_xlabel("STR allele")
            f.set_ylabel("Percentage (%)")

            fig = f.get_figure()
            fig.savefig( f"{DIR_alleleSizeDist_subdirectory}/{repeat_unit}x{reference_STR_allele}_{locus}.png" )
            plt.clf()

                # (File 2) Allele length distribution text file
            with open(f"{DIR_alleleSizeDist_subdirectory}/{repeat_unit}x{reference_STR_allele}_{locus}.txt", "w") as log:
                log.write( f"allele\tfrequency\n" )
                for allele, freq in observed_allele_histogram.items():
                    log.write(f"{allele}\t{freq}\n")
        
        # (2-5) Calculate relative allele size 
        int_relative_allele_size    = genotyped_allele - reference_STR_allele
        histogram_area_relative_to_reference  = 0
        for allele, observed_freqeuncy in observed_allele_histogram.items():
            histogram_area_relative_to_reference += (allele - reference_STR_allele) * observed_freqeuncy
        histogram_area_relative_to_reference /= 100

        list_of_lociInfos.append( [ 
            chromosome, start, end,             # locus position
            locus,                              # locus position (string; e.g., chr1:54123-55123)
            repeat_unit,                        # repeat unit (motif)
            reference_STR_allele,               # reference allele
            perc_ref_allele,                    # percentage_of_reference_allele
            allele_count.sum() / len( df_at_of_choice ), # correction rate
            allele_count.sum( ),                # corrected read count
            genotyped_allele,                   # STR genotype result
            majorPeak,                          # most prominent peak
            prominence,                         # peak prominence
            len(peaks[0]),                      # number of peaks
            int_relative_allele_size,           # relative allele size (allele - reference allele)
            histogram_area_relative_to_reference# relative allele size (histogram area relative to reference allele)
        ])
    
    STR_locus_table = pd.DataFrame( list_of_lociInfos, 
                                  columns = [ 'chromosome', 'start', 'end', 'locus', 'repeat_unit', 'reference_STR_allele', 'percentage_of_reference_allele', 'correction_rate', 
                                             'corrected_read_count', "allele", "peak_pos", "peak_prominence", "num_peaks", "int_relative_allele_size", "histogram_area_relative_to_reference"] )
    
    STR_locus_table.to_csv(f"{DIR_multiprocess_out}/{filename}.getLocusTable.temp_{batch_id+1}.tsv", sep='\t', index=False)
    return

def chunk_STR_allele_table(STR_allele_table, num_loci, threads):
    loci_per_thread = math.ceil( num_loci / threads )
    
    list_chunk  = list()
    chunk_temp  = list()
    loci_count  = 0
    
    for locus, edf in STR_allele_table.groupby("locus"):
        chunk_temp.append( edf )
        loci_count += 1

        if loci_count % loci_per_thread == 0:
            loci_count = 0 
            list_chunk.append( pd.concat( chunk_temp ) )
            chunk_temp = list()
    try:
        list_chunk.append( pd.concat( chunk_temp ) )
    except ValueError:
        pass
    
    return list_chunk

def runGetLocusTable( STR_allele_table, filename, threads,
        get_allele_size_info, read_selection, 
        min_coverage, start_time, PATH_log, DIR_out ):
    
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    
    # (1) Create a new folder for (a) multiprocessing and (b) saving STR loci allele size distribution files
    DIR_multiprocess_out = f"{DIR_out}/{filename}.multiprocessing_temp"
    nanomnt_utility.checkAndCreate( DIR_multiprocess_out )    

    if get_allele_size_info == True:
        DIR_alleleSizeDist_out = f"{DIR_out}/{filename}.allele_size_distribution"
        logging.info(f"Creating directories(s) for saving allele size distribution files: {DIR_alleleSizeDist_out}/*")
        logging.info(f"Do NOT delete or create new folder in this directory!")
        nanomnt_utility.checkAndCreate( DIR_alleleSizeDist_out )
        
        # Create subdirectory for allele length dist info files: e.g., DIR_out/repeat_unit/chromosome/
        STR_allele_table["chromosome"] = [ locus.split(":")[0] for locus in STR_allele_table.locus ]
        for repeat_unit, edf in STR_allele_table.groupby("repeat_unit"):
            for chromosome, edf2 in edf.groupby("chromosome"):
                DIR_alleleSizeDist_subdirectory_A = f"{DIR_alleleSizeDist_out}/{repeat_unit}"
                DIR_alleleSizeDist_subdirectory_B = f"{DIR_alleleSizeDist_out}/{repeat_unit}/{chromosome}"
                nanomnt_utility.checkAndCreate(DIR_alleleSizeDist_subdirectory_A)
                nanomnt_utility.checkAndCreate(DIR_alleleSizeDist_subdirectory_B)
    else:
        DIR_alleleSizeDist_out = None

    # (2) Chunk STR allele table
    num_loci    = len(set(STR_allele_table.locus))
    list_STR_allele_table_chunks = chunk_STR_allele_table( STR_allele_table, num_loci, threads )

    # (3) Use multiprocessing to iterate thorugh STR_allele_table and build STR_locus_table for each locus
    # Start multiprocessing
    processes = list()
    logging.info(f"Starting multiprocessing with {threads}\t(elapsed time: {nanomnt_utility.getElapsedTime(start_time)} seconds)")
    for idx, STR_allele_table_chunk in enumerate(list_STR_allele_table_chunks):
        p = multiprocessing.Process(target=process_STR_locus, 
                                    args=[STR_allele_table_chunk, idx, read_selection, min_coverage, filename, PATH_log, 
                                          get_allele_size_info, DIR_multiprocess_out, DIR_alleleSizeDist_out] )
        p.start()
        processes.append(p)
    for process in processes:
        process.join()

    del list_STR_allele_table_chunks 
    gc.collect()

    # (4) Merge temporary STR locus table from multiprocessing to create final STR locus table 
    logging.info(f"Finished multiprocessing. Merging temporary files.\t(elapsed time: {nanomnt_utility.getElapsedTime(start_time)} seconds)")
    list_PATH_tsvs   = glob.glob(f"{DIR_multiprocess_out}/{filename}.getLocusTable.temp_*.tsv")
    STR_locus_table  = list()
    for PATH_tsv in list_PATH_tsvs:
        STR_locus_table.append( pd.read_csv( PATH_tsv, sep='\t' ) )
        os.remove(PATH_tsv)
    os.rmdir(DIR_multiprocess_out)
    
    STR_locus_table = pd.concat( STR_locus_table )
    STR_locus_table["read_selection"] = read_selection

    # logging.info(f"Creating KDE plot of peak prominences")
    # for repeat_unit, edf in STR_locus_table.groupby("repeat_unit"):
    #     f = sns.kdeplot(x=edf.reset_index(drop=True)["peak_prominence"], bw_adjust=0.2, fill=True)
            
    #     f.set(title=f"{repeat_unit}-repeat STR (n={len(edf)})")
    #     f.set_xlabel("Peak promience")
    #     f.set_ylabel("Relative density")
    #     fig = f.get_figure()
    #     fig.savefig( f"{DIR_out}/{repeat_unit}_repeat_prominence_KDE.png" )
    #     plt.clf()

    logging.info(f"Writing STR locus table to disk")
    STR_locus_table.to_csv(f"{DIR_out}/{filename}.STR_locus_table.tsv", sep='\t', index=False)
    
    # (5) Erase any empty subdirectories
    if DIR_alleleSizeDist_out != None:
        logging.info(f"Deleting empty directories in: {DIR_alleleSizeDist_out}/*")
        if get_allele_size_info == True:
            try:
                empty_directories = nanomnt_utility.find_empty_directories( DIR_alleleSizeDist_out )
                for directory in empty_directories:
                    os.rmdir( directory )
            except FileNotFoundError:
                logging.warn( f"Could not delete empty directories in: {DIR_alleleSizeDist_out}/*. This is usually not a problem" )
        
    return 

def main():
    start_time  = time.time()

    github_documentation_link = ""

    script_description = f"Using STR allele table, generate STR locus table and allele length distribution information. For more information, see GitHub documentation: {github_documentation_link}"
    parser = argparse.ArgumentParser(description=script_description)
    
    PATH_script     = os.path.dirname(os.path.dirname(__file__))

    # Required parameters
    parser.add_argument('-at', '--PATH_STR_allele_table',  
                        help="PATH to STR allele table",
                        required=True
                        )
    parser.add_argument('-read_selection',
                        help="Read selection approach for STR genotyping: (1) all_reads (2) fw_rv_reads )", 
                        required=True
                        )    
    # Optional parameters
    parser.add_argument('--get_allele_size_info',  
                        help="Generate allele size distribution information of each locus (default: False)", 
                        action='store_true'
                        )    
    parser.add_argument('-cov', '--min_coverage',       
                        help="Minimum coverage required for analyzing STR locus. Loci with coverage lower than this value will not be processed (default: 20)", 
                        required=False, type=int, default=20)
    parser.add_argument('-t', '--threads',              
                        help="Number of threads (default: 4)", 
                        required=False, 
                        type=int, 
                        default=4
                        )
    parser.add_argument('-out', '--DIR_out',          
                        help='Directory to write the output files (default: current directory)', 
                        required=False, 
                        type=str, 
                        default=os.getcwd() 
                        )

    args = vars(parser.parse_args())
    
    PATH_STR_allele_table       = args["PATH_STR_allele_table"]
    read_selection              = args["read_selection"]
    get_allele_size_info        = args["get_allele_size_info"]
    min_coverage        = args["min_coverage"]
    threads             = args["threads"]
    DIR_out            = args["DIR_out"]
    
    filename    = os.path.splitext( os.path.basename( PATH_STR_allele_table ))[0]
    nanomnt_utility.checkAndCreate( DIR_out )
    
    PATH_log = f"{DIR_out}/nanomnt.getLocusTable.{filename}.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing parameters:")
    for k, v in args.items():
        logging.info(f'{k}\t:\t{v}')
    
    if read_selection not in ["fw_rv_reads", "all_reads"]:
        logging.error(f"Invalid read_selection input: {read_selection} must be either (1) fw_rv_reads, (2) all_reads")
        raise ValueError

    STR_allele_table = pd.read_csv( PATH_STR_allele_table, sep='\t' )

    runGetLocusTable( STR_allele_table, filename, threads,
        get_allele_size_info, read_selection, 
        min_coverage, start_time, PATH_log, DIR_out )
    
    logging.info(f"Finished getLocusTable.py\t(Total time taken: {nanomnt_utility.getElapsedTime(start_time)} seconds)")

if __name__ == "__main__":
    main()