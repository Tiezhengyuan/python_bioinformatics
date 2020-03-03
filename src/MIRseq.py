#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
#last revision  2018/09/07
Created on Mon Jul 23 12:31:01 2018
@author: Tiezheng Yuan
version v 1.0
"""

import argparse
import os
import sys
import subprocess
#
import altriaseq.utils.basic as ab
import altriaseq.utils.genome as ag
import altriaseq.utils.outer as ao

###############################################################################
def pass_par():
    #pass arguments
    parser=argparse.ArgumentParser(description="A pipeline for miRNA-seq analysis. 05.2019")
    parser.add_argument('-o', dest='dir_results', action='store', required=True, 
                        help='Path of results')
    parser.add_argument('-q', dest='dir_rawdata', action='store', required=True, 
                        help='miRNA-seq data in fastq format')
    #parser.add_argument('-s', dest='paired', action='store_true', help='Single-end sequencing')
    parser.add_argument('-t', dest='threads_num', type=int, default=24, 
                        help='Number of multi-threads')
    parser.add_argument('--adapter-3end', dest='adapter_3end', action='store', 
                        default='TGGAATTC', help='Sequences of the 3-end adapter')
    parser.add_argument('-f', dest='genome_fa_file', action='store', 
                        help='Genome fasta file')    
    parser.add_argument('--no-S1', dest='S1_alignment', action='store_false', help='Skip step 1')
    parser.add_argument('--no-S2', dest='S2_merge', action='store_false', help='Skip step 2')
    args=parser.parse_args()
    #for key,value in vars(args).items():
    #    par[key]= value
    #print(par)
    
    #genome index
    args.dir_src=os.path.dirname(os.path.realpath(__file__))+'/'
    args.dir_altriaseq=os.path.abspath(os.path.join(args.dir_src, os.pardir))+'/'
    if args.genome_fa_file is not None and os.path.isfile(args.genome_fa_file):
        args.dir_genome=os.path.dirname(args.genome_fa_file)+'/'
        args.genome_index=os.path.splitext(args.genome_fa_file)[0]
    else:
        args.dir_genome=args.dir_altriaseq+'NT3.1/'
        args.genome_index=args.dir_genome+'NT3.1-PSC'
        args.genome_fa_file=args.genome_index+'.fa'
    
    #
    args.dir_results = ab.basic().format_dir(args.dir_results)
    args.file_samples=args.dir_results+'sample_info.csv'
    args.dir_rawdata = ab.basic().format_dir(args.dir_rawdata)
    args.paired=False

    #
    if not os.path.isdir(args.dir_rawdata):
        print('Error Input! No FASTQ directory.\n')
        sys.exit()  
    if not os.path.isfile(args.file_samples):
        print('Error Input! No sample_info.csv.\n')
        sys.exit() 
    return args

##############################################################################
        
def main(args):
    #uncompress *fz file if there are
    ao.tool(args).run_pigz(args.dir_rawdata+"*.gz", False)
             
    # get sample ino
    args.sample_R1, args.sample_R2=ag.genome(args).sample_info()
    args.sample_names=args.sample_R1.keys()
    ab.basic().print_dict(args.sample_R1)


    #1:alignment
    if args.S1_alignment:
        ab.basic(args).pp_map_threads(ao.tool(args).run_mirdeep2, args.sample_names)
        #for sample_name in args.sample_names:
        #    ao.tool(args).run_mirdeep2(sample_name)  
            
    #2: merge
    if args.S2_merge:
        ag.genome(args).merge_mirdeep2()
#####################################




##############################
###parameters
args=pass_par()

##################################
main(args)

print('Great! Done!')
#end

