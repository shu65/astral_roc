'''
Created on 2016/02/25

@author: shu
''' 
import sys
from optparse import OptionParser
from operator import itemgetter, attrgetter

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from lib.astral_utility import get_classification


def get_options():
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input_filename",
                  help="input ASTRAL fasta", metavar="FILE")
    
    parser.add_option("-o", "--output_filename", dest="output_filename",
                  help="Output fasta filename", metavar="FILE")
    
        
    return parser.parse_args()


def get_key(record):
    classification = get_classification(record.description)
    return (classification["class"], classification["fold"] , classification["superfamily"])

def generate_queries(input_handle, output_handle):
    file_format = "fasta"
    counts = {}
    for record in SeqIO.parse(input_handle, file_format):
        key = get_key(record)
        #print record
        if counts.has_key(key) :
            counts[key] += 1
        else :
            counts[key] = 1
        
    
    input_handle.seek(0)
    for record in SeqIO.parse(input_handle, file_format):
        key = get_key(record)
        if counts[key] > 1 :
            SeqIO.write(record, output_handle, file_format)
        

if __name__ == '__main__':
    (options, args) = get_options()
    
    input_handle = open(options.input_filename)
    output_handle = open(options.output_filename, "w")
    generate_queries(input_handle, output_handle)

                
    input_handle.close()
    output_handle.close()
    