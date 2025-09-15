import math
import pysam
import pickle

allele_histogram_width = 20
dict_repeatUnit_to_favorableFlag = {
    "A" : 0,
    "T" : 16,
}

def calcHistogramDistance( allele_histogram_1, allele_histogram_2, method="euclidean"):
    
    key_1 = set(allele_histogram_1.keys())
    key_2 = set(allele_histogram_2.keys())
    if key_1 != key_2:
        for allele in key_1.union(key_2):
            if allele not in key_1:
                allele_histogram_1[allele] = 0 
            elif allele not in key_2:
                allele_histogram_2[allele] = 0  
    
    dist = 0 
    
    if method == "euclidean":
        for allele in allele_histogram_1.keys():
            freq_1 = allele_histogram_1[allele]
            freq_2 = allele_histogram_2[allele]
            dist += (freq_1 - freq_2)**2
        return math.sqrt( dist )
    
    elif method == "pearson":
        from scipy.stats import pearsonr  
        allele_histogram_1 = dict(sorted(allele_histogram_1.items(), key=lambda x:x[0]))
        allele_histogram_2 = dict(sorted(allele_histogram_2.items(), key=lambda x:x[0]))
        return 1 - pearsonr( list(allele_histogram_1.values()), list(allele_histogram_2.values()) )[0]
    
    elif method == "cosine":
        from scipy.spatial import distance
        allele_histogram_1 = dict(sorted(allele_histogram_1.items(), key=lambda x:x[0]))
        allele_histogram_2 = dict(sorted(allele_histogram_2.items(), key=lambda x:x[0]))
        return distance.cosine(list(allele_histogram_1.values()), list(allele_histogram_2.values()))
    
    raise ValueError

def getAlleleHistogram(df_STR_table_oi, filling_range = 0):
    repeat_count = list()
    repeat_unit = df_STR_table_oi.iloc[0].repeat_unit 
    for ca in df_STR_table_oi[df_STR_table_oi["corrected_allele"] != -1].corrected_allele:
        repeat_count.append( ca.upper().count(repeat_unit) )

    total_read_count = len(repeat_count)
    dict_repeat_count = {e : 100 * repeat_count.count(e) / total_read_count for e in repeat_count}

    dict_repeat_count_min_key = min(dict_repeat_count.keys())
    dict_repeat_count_max_key = max(dict_repeat_count.keys())
    for i in range(dict_repeat_count_min_key, dict_repeat_count_max_key):
        if i not in dict_repeat_count.keys():
            dict_repeat_count[i] = float(0)

    if filling_range == 0:
        return dict(sorted(dict_repeat_count.items()))
    else:
        dict_repeat_count = dict(sorted(dict_repeat_count.items()))
        reference_STR_allele = int( df_STR_table_oi.iloc[0].reference_STR_allele )   #! Change reference_MS_length to reference_STR_allele
    
        dict_repeat_count_n = dict()
        for filling_allele in range( reference_STR_allele - filling_range, reference_STR_allele + filling_range + 1 ):
            # if filling_allele < 0:
            #     continue
            
            if filling_allele in dict_repeat_count.keys():
                dict_repeat_count_n[filling_allele] = dict_repeat_count[filling_allele]
            else:
                dict_repeat_count_n[filling_allele] = float(0)
        
        return dict_repeat_count_n

def saveWithPickle(obj, PATH_out, filename="saveWithPickle"):
    import pickle
    with open(f'{PATH_out}/{filename}.pickle', 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

def loadFromPickle(dir_pickle):
    import pickle
    with open(dir_pickle, 'rb') as handle:
        unserialized_pickle = pickle.load(handle)
    return unserialized_pickle

def getCurDate(short=True):
    """
    description:    Returns the current time
    return:         current time (type:string)
    """
    import datetime
    date            = datetime.datetime.now()
    logging_time    = f"{date.year}-{date.month}-{date.day}_{date.hour}-{date.minute}-{date.second}-{date.microsecond}"
    if short == True:
        return "-".join( logging_time.split("-")[:-1] )
    else:
        return logging_time

def hamming_distance(chaine1, chaine2): # https://stackoverflow.com/questions/54172831/hamming-distance-between-two-strings-in-python 
    """
    description:    Given two strings, calculate the hamming 
                    distance of the two strings and return it
    return:         int
    """
    # return len(list(filter(lambda x : ord(x[0])^ord(x[1]), zip(chaine1, chaine2))))
    return len(list(filter(lambda x : ord(x[0])^ord(x[1]), zip(chaine1, chaine2)))) + abs(len(chaine1) - len(chaine2)) # The last abs(len(chaine1) - len(chaine2)) was added

def SAM_to_BAM( dir_sam, filename=None ):
    """
    description:    Given a directory of a SAM file, convert it to BAM file using samtools 
                    via the subprocess module, and delete the original SAM file.
    return:         None
    """
    import os
    import subprocess

    if filename == None:
        filename = os.path.splitext(os.path.basename(dir_sam))[0]

    subprocess.call(f"samtools view -Sb {dir_sam} -o {filename}.bam", shell=True)
    os.remove( f"{dir_sam}")
    return

def sort_and_index_BAM( dir_bam, filename=None ):
    """
    description:    Given a directory of a BAM file, sort and index it using samtools 
                    via the subprocess module, and delete the original BAM file.
    return:         None
    """
    import os
    import subprocess

    if filename == None:
        filename = os.path.splitext(os.path.basename(dir_bam))[0]

    subprocess.call(f"samtools sort {dir_bam} -o {filename}.sorted.bam", shell=True)
    subprocess.call(f"samtools index {filename}.sorted.bam", shell=True)
    os.remove( f"{dir_bam}")
    return
  
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

def check_empty_PATH( dir_PATH ):
    """
    description:    Given a PATH, check if it's empty
    return:         bool or None
    """
    import os
    try:
        return not next(os.scandir(dir_PATH), None)
    except FileNotFoundError:
        return None

def find_empty_directories(directory):
    """ 
    description:    Given a PATH, recursively check if it and its child PATHs 
                    are empty, and return PATHs that aren't empty.
    return:         list
    """
    import os
    empty_directories = []

    # Check if directory is empty
    if not any(os.scandir(directory)):
        empty_directories.append(directory)

    # Recursively check child directories
    for item in os.scandir(directory):
        if item.is_dir():
            empty_directories.extend(find_empty_directories(item.path))

    return empty_directories