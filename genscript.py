#!/usr/bin/env python
import argparse

'''read data list and generate pipline script'''
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",help="file to be copied")
parser.add_argument("-p1","--pathofdatabase",help="paht of database")
parser.add_argument("-p2","--pathofdata",help="path of archieved data")
args = parser.parse_args()

filename = args.input
path_database = args.pathofdatabase
path_data = args.pathofdata
path_pipeline = "/home/hxmt/tuoyl/HXMT_process/hxmtanalysis/he_craftdata_anal.py"

with open(filename,'r')as f:
    file = f.read()
file = file.split('\n')[0:-1]
with open(filename+'_script.sh','w')as f:
    for i in file:
        write_str = 'python %s -i %s -o %s --hxbary --clobber --fig\n'%(path_pipeline, path_database+i, path_data+i)
        f.write(write_str)
