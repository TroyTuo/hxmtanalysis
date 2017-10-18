#!/usr/bin/env python
import argparse

'''read data list and generate pipline script'''
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",help="data archived path")
parser.add_argument("-n","--name",help="source name")
args = parser.parse_args()
filename = args.input
sourcename = args.name
with open(filename,'r')as f:
    file = f.read()
file = file.split('\n')[0:-1]
with open(filename+'_script.sh','w')as f:
    for i in file:
        write_str = 'python   
