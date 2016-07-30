'''
Created on 2016/02/26

@author: shu
'''

import sys
from optparse import OptionParser

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from src.lib.astral_utility import get_classification
from src.lib.blast_result_reader import BlastResultReader
import math
import json

def get_options():
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input_filename",
                  help="the output of BLAST", metavar="FILE")
    
    parser.add_option("-a", "--astral_filename", dest="astral_filename",
                  help="ASTRAL filename", metavar="FILE")
    
    parser.add_option("-n", "--number_of_false_positives", dest="number_false_positives",
                  help="the number of false positives",  action="store", type="int", default=5000)
    
    parser.add_option("-r", "--output_raw_data", dest="output_raw_data",
                  help="output raw data for ROC",  action="store_true")
        
    return parser.parse_args()


def is_true_positive(sequence1_classification, sequence2_classification):
    ret = sequence1_classification["class"] == sequence2_classification["class"]
    ret = ret and sequence1_classification["fold"] == sequence2_classification["fold"] 
    ret = ret and sequence2_classification["fold"] and sequence1_classification["superfamily"] == sequence2_classification["superfamily"]
    return ret

def is_false_positive(sequence1_classification, sequence2_classification):
    ret = sequence1_classification["class"] != sequence2_classification["class"]
    ret = ret or sequence1_classification["fold"] != sequence2_classification["fold"]
    return ret

        
def build_classification_map(astral_file_handle):
    file_format = "fasta"
    classification_map = {}
    for record in SeqIO.parse(astral_file_handle, file_format):
        classification = get_classification(record.description)
        classification_map[record.id] = classification

    return classification_map

def get_key(classification):
    return (classification["class"], classification["fold"] , classification["superfamily"])


def count_true_positive(classification_map):
    counts= {}
    for classification in classification_map.values() :
        key = get_key(classification)
        #print record
        if counts.has_key(key) :
            counts[key] += 1
        else :
            counts[key] = 1
    
    T = 0
    for count in counts.values() :
        T += count*(count -1)
    return T

def compute_roc_mean(n, T, true_positive_counts):
    true_positive_count_sum = 0
    for i in range(n):
        true_positive_count_sum += true_positive_counts[i]
        
    return float(true_positive_count_sum)/float(n*T)

def compute_roc_sigma(n, T, true_positive_counts):
    print n, len(true_positive_counts)
    assert((n + 1) == len(true_positive_counts))
    true_positive_count_n_1 = true_positive_counts[n]
    sum_difference_pow2 = 0
    for i in range(n):
        sum_difference_pow2 += (true_positive_count_n_1 - true_positive_counts[i])**2
        
    return math.sqrt(float(sum_difference_pow2)/(float(n*T))**2)

def compute_roc_score(blast_records, classification_map, n, T):    
    sorted_blast_records = sorted(blast_records, key=lambda x:x["e_value"])
    
    true_positive_counts = []
    current_true_positive_count = 0
    false_positive_count = 0
    for record in sorted_blast_records: 
        #print record["e_value"]
        if false_positive_count > n:
            break
        
        query_name = record["query_name"]
        subject_name = record["subject_name"]
        if query_name == subject_name :
            continue
        
        query_classification = classification_map[query_name]
        subject_classification = classification_map[subject_name]
        #print query_name, "(", query_classification, ")",  subject_name, "(", subject_classification, ")", record["e_value"]

        
        if is_true_positive(query_classification, subject_classification) :
            current_true_positive_count += 1
            #print True
            
        elif is_false_positive(query_classification, subject_classification) :
            false_positive_count += 1
            true_positive_counts.append(current_true_positive_count)
            #print False
            
        else :
            pass
            #print None

    if false_positive_count < (n + 1) :
        for i in range(false_positive_count, n + 1) :
            true_positive_counts.append(current_true_positive_count)

    mean_roc = compute_roc_mean(n, T, true_positive_counts)
    sigma_roc = compute_roc_sigma(n, T, true_positive_counts)
    return mean_roc, sigma_roc


def compute_roc_curve(blast_records, classification_map, max_n, T):    
    sorted_blast_records = sorted(blast_records, key=lambda x:x["e_value"])
    
    true_positive_counts = []
    current_true_positive_count = 0
    false_positive_count = 0
    hit_checking = {}
    duplicate_count = 0
    for record in sorted_blast_records: 
        #print record["e_value"]
        if false_positive_count >= max_n:
            break
        
        query_name = record["query_name"]
        subject_name = record["subject_name"]
        if query_name == subject_name :
            continue
        
        if hit_checking.has_key((query_name, subject_name)) :
            duplicate_count += 1
            #print "hit duplicate :" , query_name,  subject_name
            continue
            #pass
        hit_checking[(query_name, subject_name)] = 1
        query_classification = classification_map[query_name]
        subject_classification = classification_map[subject_name]
        #print query_name, "(", query_classification, ")",  subject_name, "(", subject_classification, ")", record["e_value"]

        
        if is_true_positive(query_classification, subject_classification) :
            current_true_positive_count += 1
            #print True
            
        elif is_false_positive(query_classification, subject_classification) :
            false_positive_count += 1
            true_positive_counts.append(current_true_positive_count)
            #print False
            
        else :
            pass
            #print None

    if false_positive_count < max_n :
        for i in range(false_positive_count, max_n) :
            true_positive_counts.append(current_true_positive_count)
    

    print "duplicate_count", duplicate_count
    assert(len(true_positive_counts) == max_n)
    return true_positive_counts
            
        
    


if __name__ == '__main__':
    (options, args) = get_options()
    astral_file_handle = open(options.astral_filename)
    classification_map = build_classification_map(astral_file_handle)
    astral_file_handle.close()
    n = options.number_false_positives
    T = count_true_positive(classification_map)
    input_handle = open(options.input_filename)

    blast_reader = BlastResultReader(input_handle)
    blast_records = []
    for recode in blast_reader :
        blast_records.append(recode)
    input_handle.close()
    
    #print "n=", n 
    #print "T=", T 
    mean_roc, sigma_roc = compute_roc_score(blast_records, classification_map, n, T)
    ret = {}
    ret["mean_roc"] = mean_roc
    ret["sigma_roc"] = sigma_roc
    if options.output_raw_data :
        ret["true_positive_counts"] = compute_roc_curve(blast_records, classification_map, n, T)
    print json.dumps(ret, indent=2)
    
