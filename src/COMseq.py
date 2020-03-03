#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 17:26:08 2019

@author: Tiezheng Yuan

The tool is a pipeline for Genome Comparison
"""

import argparse
import copy
import os
import sys

import altriaseq.utils.basic as ab
import altriaseq.utils.circos as ac
import altriaseq.utils.genome as ag
import altriaseq.utils.biofile as bf
import altriaseq.utils.outer as ao

###############################################################################
def pass_par():
    #pass arguments
    parser = argparse.ArgumentParser(description="A pipeline for genome comparison. v1.0, 2019")
    parser.add_argument('-n',dest='genome_names', action='store', required=True, 
                        help='Genome names')
    parser.add_argument('-o',dest='dir_results', action='store', required=True, 
                        help='Path of results')
    parser.add_argument('-f',dest='fa_files', action='store', required=True, 
                        help='genome sequences in fasta format')
    parser.add_argument('-g',dest='gff_files', action='store', 
                        help='genome sequences in GFF format')
    parser.add_argument('-a',dest='faa_files', action='store', 
                        help='genome amino acid sequences in FASTA format')
    parser.add_argument('-s',dest='min_size', default=10, type=int, action='store', 
                        help='The minimum of genes required for a collinear block')
    parser.add_argument('--gff-tag-name', dest='gff_tag_name', default='ID', action='store', 
        help='The tag name for eliciting geneID in column 9 from GFF file')
    parser.add_argument('--no-S1', dest='S1_annotation', action='store_false', 
                        help='Skip step 1 for debugging')
    parser.add_argument('--no-S2', dest='S2_compare', action='store_false', 
                        help='Skip step 2 for debugging')
    parser.add_argument('--no-S3', dest='S3_collinear', action='store_false', 
                        help='Skip step 3 for debugging')
    args=parser.parse_args()
    #print(args)

    #
    args.dir_src=os.path.dirname(os.path.realpath(__file__))+'/'
    args.dir_altriaseq=os.path.abspath(os.path.join(args.dir_src, os.pardir))+'/'
    #default directories
    args.dir_results = ab.basic().format_dir(args.dir_results)
    if not os.path.isdir(args.dir_results):
        print('Error! No input path for storing results: {}\n'.format(args.dir_results))
        sys.exit()
    #genome pairs
    args.genome_names=args.genome_names.split(',')
    args.fa_files=args.fa_files.split(',')
    if len(args.genome_names)!=len(args.fa_files):
        print('Error Input! Wrong input of genome pairs!\n')
        sys.exit()
    #
    if args.gff_files is None:
        args.gff_files=[]
    else:
        args.gff_files=args.gff_files.split(',')
    if args.faa_files is None:
        args.faa_files=[]
    else:
        args.faa_files=args.faa_files.split(',')
 
    #set dir for each genome:
    args_dict={}
    for i, name in enumerate(args.genome_names):
        sub_args=copy.deepcopy(args)
        if i <len(args.fa_files):
            sub_args.file_fa=args.fa_files[i]
        if i <len(args.gff_files):
            sub_args.file_gff=args.gff_files[i]
        if i <len(args.faa_files):
            sub_args.file_faa=args.faa_files[i]
        sub_args.chr_name = name
        sub_args.genome_name = name
        sub_args.dir_annot = ab.basic(args).format_dir(args.dir_results+name)
        sub_args.dir_circos = sub_args.dir_annot
        args_dict[name]=sub_args
    args.sub_args=args_dict
        
    
    #
    #for key,value in sorted(vars(args).items()):
    #    print(key,':\t', value)
    
    return args

    
    
###############################################################################
def main(args):
                
    #1: genome annotation
    if args.S1_annotation:
        print('\t###Annotate genome and visualization for each genome.')
        for i, name in enumerate(args.genome_names):
            args.file_fa=args.fa_files[i]
            #1.1: run prokka if not gff files
            if  len(args.gff_files)==0:
                args.dir_annot = ab.basic(args).format_dir(args.dir_results+name)
                ao.tool(args).run_prokka()
            #1.2: ciros plotting
            args.dir_circos=args.dir_annot
            ac.circos(args.sub_args[name]).one_genome()

    
    ##2:genome comparison
    if len(args.genome_names)==1:
        print('Only one genome. No genome comparison but jut draw genome using circos.\n')
        sys.exit()
    if args.S2_compare:
        #2.1: run mauve for genome DNA alignments
        args.file_backbone=args.dir_results+'mauve/backbone.txt'
        args.dir_mauve = ab.basic(args).format_dir(args.dir_results+'mauve')
        ao.tool(args).run_mauve()  

        #2.2: run circos for multiple based on mauve results
        if os.path.isfile(args.file_backbone):
            args.in_circos=[ args.dir_results+x+'/' for x in args.genome_names]
            args.dir_circos = ab.basic(args).format_dir(args.dir_results+'circos/')
            ac.circos(args).multiple_genome()
        
    #3:  collinear
    if args.S3_collinear:
        #update annotation
        args.sub_args=ag.genome(args).elicit_annotation()
        
        #3.1: run MCScanX
        args.dir_collinear = ab.basic(args).format_dir(args.dir_results+'synteny/')
        args.collinear_prefix=args.dir_collinear+'_'.join(args.genome_names)
        ao.tool(args).run_mcscanx()
        
        #3.2: build control files and draw dot plot
        ao.tool(args).run_mcscanx_plotter()
        
        #3.3: draw circos
        ac.circos(args).collinear()

#        

################################

##############################
###parameters
args=pass_par()

##################################
##run main script
main(args)

print('Great! Done!')
#end