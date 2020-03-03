#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
#last revision  2018/09/07
Created on Mon Jul 23 12:31:01 2018
Last revision June 2019

@author: Tiezheng Yuan
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
    parser=argparse.ArgumentParser(description="A pipeline for mRNA-seq using hisat2 and stringtie. 07.2018-06.2019")
    parser.add_argument('-o',dest='dir_results', action='store', required=True, help='Path of results')
    parser.add_argument('-q',dest='dir_rawdata', action='store', required=True, help='mRNA-seq data in fastq format')
    parser.add_argument('-s',dest='single_end', action='store_true', help='single-end sequencing')
    parser.add_argument('-t',dest='threads_num', type=int, default=24, help='Number of multi-threads')
    parser.add_argument('--no-S1', dest='S1_alignment', action='store_false', help='Skip step 1')
    parser.add_argument('--no-S2', dest='S2_assembly', action='store_false', help='Skip step 2')
    parser.add_argument('--no-S3', dest='S3_ballgown', action='store_false', help='Skip step 3')
    args=parser.parse_args()
    #for key,value in vars(args).items():
    #    par[key]= value
    #print(par)
    
    #genome index
    args.dir_src=os.path.dirname(os.path.realpath(__file__))+'/'
    args.dir_altriaseq=os.path.abspath(os.path.join(args.dir_src, os.pardir))+'/'
    args.genome_index=args.dir_altriaseq+'NT3.1/NT3.1-PSC'
    #args.genome_index='/mnt/rdedata12/NT3.1/NT3.1-PSC'
    args.genome_gtf_file=args.genome_index+'.gtf'
    args.genome_fa_file=args.genome_index+'.fasta'
  

    #
    args.dir_results = ab.basic().format_dir(args.dir_results)
    args.file_samples=args.dir_results+'sample_info.csv'
    if not os.path.isfile(args.file_samples):
        print('Error Input! No sample_info.csv in ', args.dir_results)
        sys.exit()  
    


    #
    args.paired= True if args.single_end is False else False
    
    #
    args.dir_rawdata = ab.basic().format_dir(args.dir_rawdata)
    if not os.path.isdir(args.dir_rawdata):
        print('Error Input! No FASTQ directory.\n')
        sys.exit()  
    
        #
    par={}
    for key,value in vars(args).items():
        par[key]= value
        print(key, ':\t', value)
    return args

##############################################################################
        
def main(args):
    #uncompress *fz file if there are
    ao.tool(args).run_pigz(args.dir_rawdata+"*.gz", False)
             
    # get sample ino
    args.sample_R1, args.sample_R2=ag.genome(args).sample_info()
    args.sample_names=args.sample_R1.keys()
    ab.basic().print_dict(args.sample_R1)
    ab.basic().print_dict(args.sample_R2)


    #1:alignment
    if args.S1_alignment:
        ab.basic(args).pp_map_threads(ao.tool(args).run_hisat2, args.sample_names)
        #for sample_name in args.sample_names:
        #    ao.tool(args).run_hisat2(sample_name)
 

    #2: 
    if args.S2_assembly:
        #3:samtootls
        for sample_name in args.sample_names:
            ao.tool(args).run_samtools(sample_name)  
            
        #4:stringtie
        ab.basic(args).pp_map_threads(ao.tool(args).run_stringtie, args.sample_names)
        #for sample_name in args.sample_names:
        #    ao.tool(args).run_stringtie(sample_name)
    
    #4:counts file
    if args.S3_ballgown:
        #5:merge
        ao.tool(args).run_merge()
        
        #multi-processes
        ab.basic(args).pp_map_threads(ao.tool(args).run_ballgown, args.sample_names)
        
        #export FPKM table
        command="Altria_ballgown {}".format(args.dir_results)
        print("@@@@", command)
        subprocess.Popen(command, stdout=subprocess.PIPE, shell=True).stdout.read() 
    #
#####################################




##############################
###parameters
args=pass_par()

##################################
main(args)

print('Great! Done!')
#end

