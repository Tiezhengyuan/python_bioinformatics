#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 14:20:18 2019
Pipeline for Nanopore sequencing data analysis
@author: Tiezheng Yuan
"""


import argparse
import os
import sys

#
#import altriaseq
#import altriaseq.utils
import altriaseq.utils.basic as ab
#import altriaseq.utils.biofile as bf
import altriaseq.utils.fast5 as af
import altriaseq.utils.genome as ag
import altriaseq.utils.outer as ao

###############################################################################
def pass_args():
    #pass arguments
    parser=argparse.ArgumentParser(description="A pipeline for Nanopore sequencing data analysis. v.1.0, 2018-2019")
    parser.add_argument('-o', dest='dir_results', action='store', required=True, 
        help='Path of results')
    parser.add_argument('-m', dest='nano_mode', action='store', type=str, required=True,
        help='-m cdna: detect novel isoforms, -m quant: quantify isoforms of mutiple samples')
    #cDNA mode
    parser.add_argument('--fast5', dest='dir_fast5', type=str, 
        help='Path of FAST5 reads sequences (That will be searched recursively)')
    parser.add_argument('-q', dest='file_fq', action='store', 
        help='sequencing data in fastq format')
    parser.add_argument('-f', dest='genome_fa_file', action='store', 
        help='Genome DNAs files in fasta format')
    parser.add_argument('-g', dest='genome_gtf_file', action='store', 
        help='Genome annotation file in gtf/gff format')
    #differential mode
    parser.add_argument('--fastq', dest='dir_fastq', type=str, 
         help='Path of FASTQ reads sequences (That will be searched recursively)')
    parser.add_argument('-c', dest='transcripts_fa_file', action='store', 
        help='Transcripts DNAs (cDNAs) in fasta format')
    '''
    #Never use! External validation
    #pre-processing
    #af.fast5(args).decompose_fast5()
    parser.add_argument('-t', dest='nano_type', type=str, default='best', 
        help='Nanopore sequence entires: best, both strands(=2D), template strand(=fwd), complement strand(=rev)')
    parser.add_argument('--min_length', dest='min_length', type=int, default=1000,
        help='Exclude reads shorter than this length (in bp)')
    parser.add_argument('--min_mean_qual', dest='min_mean_qual', type=float, default=9, 
        help='Exclude reads with a mean qscore less than this value')
    parser.add_argument('--min_qual_window', dest='min_qual_window', type=float, default=0, 
        help='Exclude reads where their mean qscore in a sliding window drops below this value')
    parser.add_argument('--window_size', dest='window_size', type=int, default=50, 
        help='The size of the sliding window used for --min_qual_window')
    #default target base>500Mbp determined by the sequencing depth
    parser.add_argument('--target_bases', dest='target_bases', type=int, 
        help='Discard the worst reads based on minimum qscore in a sliding window') 
    '''
    parser.add_argument('--no-S1', dest='S1_pre', action='store_false', help='Skip step 1')
    parser.add_argument('--no-S2', dest='S2_process', action='store_false', help='Skip step 2')
    parser.add_argument('--no-S3', dest='S3_post', action='store_false', help='Skip step 3')
    args=parser.parse_args()
    
    #update namespace
    if args.dir_results is None:
        print("Error! the parameter -o is missing")
        sys.exit()
    else:
        vars(args)['dir_results'] = ab.basic().format_dir(args.dir_results)
    
    #genome index
    args.dir_src=os.path.dirname(os.path.realpath(__file__))+'/'
    args.dir_altriaseq=os.path.abspath(os.path.join(args.dir_src, os.pardir))+'/'
    args.dir_exe=args.dir_altriaseq+'exe/'
    if args.nano_mode == 'cdna':
        if args.genome_fa_file is not None and os.path.isfile(args.genome_fa_file):
            args.genome_index= os.path.splitext(args.genome_fa_file)[0]
        else:
            args.genome_index=args.dir_altriaseq+'NT3.1/NT3.1-PSC'
            args.genome_fa_file=args.genome_index+'.fasta'
        if args.genome_gtf_file is None:
            args.genome_gtf_file=args.genome_index+'.gtf'
        #default references
        args.ref_file=args.genome_fa_file
        #fast5 or fastq
        if args.dir_fast5 is None:
            if not os.path.isfile(args.file_fq):
                print("Error! Need to specify -q (FASTQ file).")
                sys.exit()
        else:
            if os.path.isdir(args.dir_fast5):
                args.dir_fast5 = os.path.abspath(args.dir_fast5)
                fast5_files = af.fast5(args).find_fast5()
                fast5_num=len(fast5_files)
                if  fast5_num>0 :
                    print("Then number of FAST5 files is {}".format(fast5_num),file=sys.stderr)
                    args.fast5_files=fast5_files
                    if args.file_fq is None:
                        args.file_fq=args.dir_results+'fastq/pass.fq'
                else:
                    print("Error! No FAST5 files are found.")
                    sys.exit()
            else:
                print("Error! Wrong path of FAST5 files.")
                sys.exit()
    #quant mode
    elif args.nano_mode == 'quant':
        #reference
        if args.transcripts_fa_file is not None and os.path.isfile(args.transcripts_fa_file):
            args.transcripts_index= os.path.splitext(args.transcripts_fa_file)[0]
        else:
            args.transcripts_index=args.dir_altriaseq+'NT3.1/NT3.1-PSC_transcripts'
            args.transcripts_fa_file=args.transcripts_index+'.fa'
        #fastq directory
        if not os.path.isdir(args.dir_fastq):
            print("Error! Need to specify --fastq.")
            sys.exit()  
        #sample information
        args.file_samples=args.dir_results+'sample_info.csv'
        if not os.path.isfile(args.file_samples):
            print('Error Input! No sample_info.csv in ', args.dir_results)
            sys.exit() 
            
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
    if args.nano_mode == 'cdna':
        if args.S1_pre:
            #convert fast5 to fasta
            if args.dir_fast5 is not None:
                ao.tool(args).run_albacore()
            #QC
            ao.tool(args).run_nanoplot()
            
        #2:discover genes and annotate their isoforms
        if args.S2_process:
            args.dir_aln=ab.basic().format_dir(args.dir_results+'alignments')
            args.dir_pinfish = ab.basic().format_dir(args.dir_results+'pinfish')
            #1-2:alignment
            args.query_file=args.dir_results+'fastq/pass.fq'
            args.bam_file=args.dir_aln+'reads_aln_sorted.bam'
            ao.tool(args).run_minimap2_cdna()
            #3-6:gff
            ao.tool(args).run_pinfish()
            #7-:polished minimap2
            args.query_file=args.dir_pinfish+"polished_transcripts.fas"
            args.bam_file = args.dir_aln+"polished_reads_aln_sorted.bam"
            ao.tool(args).run_minimap2_cdna()        
            #8-11:gff
            #input bam file : polished_reads_aln_sorted.bam
            ao.tool(args).run_pinfish_polished()

        #3:
        if args.S3_post:
            ao.tool(args).run_gffcompare()
    
    #quantitative isoforms analysis
    elif args.nano_mode == 'quant':
        #1:#convert fast5 to fastq
        if args.S1_pre:
            args.file_fq=args.dir_results+'fastq/classified.fq'
            #1.1: demultiplexing
            ao.tool(args).run_guppy_barcoder()
            #1.2: QC
            ao.tool(args).run_nanoplot()
            
        #2:discover genes and annotate theri isoforms
        if args.S2_process:
            #retrieve fastq files
            args.dir_rawdata=args.dir_results+'fastq'
            sample_R1, sample_R2=ag.genome(args).sample_info()
            ab.basic().print_dict(sample_R1)
            args.sample_names=sample_R1.keys()
            #
            args.ref_file=args.transcripts_fa_file
            #sequence alignment
            for sample in args.sample_names:
                args.dir_sample=ab.basic().format_dir(args.dir_results+sample)
                #2.1:build index, alignment, bam conversion
                args.query_file=sample_R1[sample]
                args.bam_file=args.dir_sample+'reads_aln_sorted.bam'
                ao.tool(args).run_minimap2_quant()
                #2.2: counts reads
                ao.tool(args).run_salmon()
            #2.3: merge data
            ag.genome(args).merge_salmon()
            

    #
##############################


##################################
###parameters
args=pass_args()

##################################
main(args)

print('Great! Done!')
#end

