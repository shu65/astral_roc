'''
Created on 2016/02/26

@author: shu
'''

def get_classification(header):
    classification_string = header.split(" ")[1].strip()
    classifications =  classification_string.split(".")
    return {"class": classifications[0], "fold":classifications[1] , "superfamily": classifications[2], "family":classifications[3]}