'''
Created on 2012/11/05

@author: shu
'''

import csv

class BlastResultReader:

    def __init__(self, file):
        self.__prev_query_name = ""
        self.__prev_ranking = -1
        self.__file = csv.reader(file, delimiter='\t')
        
    def __iter__(self):
        return self

    def next(self):
        row = self.__file.next()
        while (((not (len(row) == 0)) and (row[0][0] == '#')) or (len(row) < 11)) :
            row = self.__file.next()
            
        if len(row) == 0:
            raise StopIteration
        if self.__prev_query_name != row[0] :
            self.__prev_query_name = row[0]
            self.__prev_ranking = 0
        
        self.__prev_ranking = self.__prev_ranking + 1
        recode = {}
        recode["query_name"] = row[0].split(" ")[0];
        recode["subject_name"] = row[1].split(" ")[0];
        recode["ranking"] = self.__prev_ranking;
        recode["query_start"] = row[6];
        recode["query_end"] = row[7];
        recode["subject_start"] = row[8];
        recode["subject_end"] = row[9];
        recode["e_value"] = float(row[10]);
        recode["identity"] = float(row[2]);
        recode["bit_score"] = float(row[11]);
        return recode
    
if __name__ == "__main__":
    blast_reader = BlastResultReader(open("/mnt/fs/suzuki/tmp/result"))
    for recode in blast_reader :
        print recode
        

        
