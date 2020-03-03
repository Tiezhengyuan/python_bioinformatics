#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 14:20:18 2019
Pipeline for statistical analysis
@author: Tiezheng Yuan
"""


import argparse
import os
import sys
import subprocess

#
import altriaseq.utils.basic as ab
###############################################################################
def pass_args():
    #pass arguments
    parser=argparse.ArgumentParser(description="A pipeline for datastatistical analysis. v.1.0, 2018-2019")
    parser.add_argument('-o', dest='dir_results', action='store', help='Path of results')
    parser.add_argument('-m', dest='mode', action='store', type=str, required=True,
        help='-m deseq: deseq')
    parser.add_argument('-i', dest='file_rc', type=str, required=True,
        help='Path of the read counts file')
    parser.add_argument('-p', dest='file_pair', action='store', required=True,
        help='Path of the file of sample names')
    args=parser.parse_args()
    #
    args.dir_src=os.path.dirname(os.path.realpath(__file__))+'/'
    args.dir_altriaseq=os.path.abspath(os.path.join(args.dir_src, os.pardir))+'/'
    args.dir_R=args.dir_altriaseq+'R/'
    
    #
    if args.file_rc is None:
        print("Error! the parameter -o is missing")
        sys.exit()
    elif not os.path.isfile(args.file_rc):
        print("Error! the input of the read counts file is wrong")
        sys.exit()
    #
    if args.file_pair is None:
        print("Error! the parameter -p is missing")
        sys.exit()
    elif not os.path.isfile(args.file_pair):
        print("Error! the input of the sample names file is wrong")
        sys.exit()
        
    #update namespace
    if args.dir_results is None:
        args.dir_results = os.path.dirname(args.file_rc)
    args.dir_results = ab.basic().format_dir(args.dir_results)
    
     
    #
    par={}
    for key,value in vars(args).items():
        par[key]= value
        #print(key, ':\t', value)
    
    #
    return args
    


###############################################################################
def main(args):
    #detect novel transcripts
    if args.mode == 'deseq':
        rscript=args.dir_R+'DESeq.R'
        command="Rscript --vanilla {} {} {} {}".format(rscript, \
            args.file_rc, args.file_pair, args.dir_results)
        print("@@@@", command)
        subprocess.Popen(command, stdout=subprocess.PIPE, shell=True).stdout.read() 
            

    #
##############################


##################################
###parameters
args=pass_args()

##################################
main(args)

print('Great! Done!')
#end

