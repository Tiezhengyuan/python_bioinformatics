#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 14:20:18 2019

Pipeline for genome DNA assembly

@author: Tiezheng Yuan
"""


import argparse
import os
import sys

#
#import altriaseq
#import altriaseq.utils
import altriaseq.utils.basic as ab
import altriaseq.utils.circos as ac
import altriaseq.utils.genome as ag
import altriaseq.utils.outer as ao

###############################################################################
def pass_par():
    ###parameters
    #pass arguments
    parser=argparse.ArgumentParser(description="A pipeline for genome assembly and annotation. v.1.0, 2018-2019")
    parser.add_argument('-o',dest='dir_results', action='store', required=True, 
                        help='Path of results')
    parser.add_argument('-q',dest='file_fq', action='store', required=True, 
                        help='genome sequences in fastq format')
    parser.add_argument('-n',dest='genome_name', action='store', default='A', help='Genome name')
    parser.add_argument('--no-S1', dest='S1', action='store_false', help='Skip step 1 for debugging')
    parser.add_argument('--no-S2', dest='S2', action='store_false', help='Skip step 2 for debugging')
    parser.add_argument('--no-S3', dest='S3', action='store_false', help='Skip step 3 for debugging')

    args=parser.parse_args()
    #for key,value in vars(args).items():
    #    par[key]= value
    #print(args)
            
    #
    if not os.path.isfile(args.file_fq):
        print('Error Input! No FASTQ files.\n')
        sys.exit()   
    #default directories
    args.dir_results = ab.basic().format_dir(args.dir_results)


    return args


###############################################################################
def main(args):
    args.dir_genome = ab.basic().format_dir(args.dir_results+'unicycler')
    args.file_fa = args.dir_genome+'assembly.fasta'
    args.dir_annot = ab.basic().format_dir(args.dir_results+'prokka')
    args.dir_circos = ab.basic(args).format_dir(args.dir_results+'circos')
    args.dir_blast = ab.basic(args).format_dir(args.dir_results+'blast')
    
    #1:
    if args.S1:
        ao.tool(args).run_nanoplot()

    #2:
    if args.S2:
        #2.1 DNA chain assembly
        ao.tool(args).run_unicycler()
        #2.2: genome annotation
        ao.tool(args).run_prokka()
        #2.3 genome assessment
        ao.tool(args).run_quast()

        
    #3circos
    if args.S3:
        #3.1: blastn
        fa_name=os.path.splitext(os.path.basename(args.file_fa))[0]
        args.blast_db = args.dir_blast +'db_' + fa_name
        args.blast_out=args.dir_blast+'tblastn_'+fa_name+'.tsv'
        ao.tool(args).run_blastn()
        
        ##3.2:circos plotting
        args.chr_name=args.genome_name
        ac.circos(args).one_genome()
    #
##############################


##################################
###parameters
args=pass_par()

##################################
main(args)

print('Great! Done!')
#end

