#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 18:07:41 2019

@author: Tiezheng Yuan
"""

import os
import sys


import altriaseq.utils.basic as ab

##############################################################################
#main
class par_mode:
    def __init__(self, par=None):
        self.par=par
        
    def select_app(self):
        mode=['exit','mRNA-Seq','Genome DNA assembly', 'Genome DNA sequence comparison']
        print('\nSelect what you want:')
        for i in range(len(mode)):
            print("\t{}:  {}".format(i, mode[i]))
        tag=0
        while(1):
            tag=input('Enter the number:')
            try:
                tag=int(tag)
            except ValueError:
                pass
            else:
                if 0<=tag<4:
                    break
        #print(tag)
        return tag
        
    def select_dir(self, instr):
        print('\nSelect the directory for {}:'.format(instr))
        while(1):
            indir=input('Enter the path:')
            if os.path.isdir(indir):
                indir=os.path.abspath(indir)+'/'
                break
            else:
                print('Wrong input. Try again!')
        return indir

    def set_dir(self, instr):
        print('\nSet the directory for {}:'.format(instr))
        while(1):
            indir=input('Enter the path:')
            if os.path.isdir(indir):
                indir=os.path.abspath(indir)
                break
            else:
                try: 
                    indir=ab.basic().format_dir(indir)
                except ValueError:
                    print('Wrong input. Try again!')
                else:
                    break
        return indir
    
    def set_name(self, instr):
        print('\n{},'.format(instr))
        while(1):
            name=input('Enter:')
            if len(name)>0:
                break
            else:
                print('Wrong input. Try again!')
        return name
        
    def select_file(self, instr):
        print('\nSet the file {}:'.format(instr))
        while(1):
            infile=input('Enter the path:')
            if os.path.isfile(infile):
                infile=os.path.abspath(infile)
                break
            else:
                print('Wrong input. Try again!')
        return infile

#
    def judge_COMseq(self):
        a=par['genome_names'].split(',')
        b=par['fa_files'].split(',')    
        tag=0
        #
        if len(a)!=len(b):
            print('Error! The genome name and its fasta file should be one one one.\n')
            tag=1
        #
        if len(a)<2 or len(b)<2:
            print('Error! At least two genomes should be set.\n')
            tag=1
        #
        for f in b:
            if not os.path.isfile(f):
                print('Error! The file {} doesnot exist.\n'.format(f))
                tag=1
        if tag==1:
            sys.exit()
        #
        
#
                
##############################################################################
par={'app':3, 'infile':None, 'indir':None, 'outdir':None, 'name':None,
     'src_dir':os.path.dirname(os.path.abspath(__file__))+'/'}

#selct running mode
par['app']=par_mode(par).select_app()

#mode
if par['app']==0:
    sys.exit()
#mRNA-seq
elif par['app']==1:
    par['indir']=par_mode(par).select_dir('storing fastq files')#fastq directory
    par['outdir']=par_mode(par).select_dir('storing results')#results directory
    if not os.path.isfile(par['outdir']+'/sample_info.csv'):
        print('Error: No file known as sample_info.csv is in', par['outdir'])
        sys.exit()
    command="yuan_RNAseq -q {} -o {}".format(par['indir'], par['outdir'])

#genome assembly
elif par['app']==2:
    par['infile']=par_mode(par).select_file('storing sequence reads')#fastq file
    par['genome_name']=par_mode(par).set_name('Set a genome name')#results
    par['outdir']=par_mode(par).select_dir('storing results')#results directory
    par['outdir']=ab.basic().format_dir(par['outdir']+par['genome_name'])
    command="yuan_GENseq -q {} -o {}".format(par['infile'], par['outdir'])
#genome comparison
elif par['app']==3:
    par['outdir']=par_mode(par).select_dir('storing results')#results directory
    par['genome_names']=par_mode(par).set_name('Enter at least two genome names seperated by comma')
    par['fa_files']=par_mode(par).set_name('Enter genome sequence in FASTA format seperated by comma')
    #judge input
    par_mode(par).judge_COMseq()
    #launch pipeline
    command="yuan_COMseq -o {} -n {} -f {}".format(par['outdir'],par['genome_names'],par['fa_files'])
else:
    print('Error!. Wrong input.\n')
    sys.exit()

###
#launch pipeline
#print(par)   
print("@@@@@@@@@@@@", command)
os.system(command)

#
print('\nDone. Great!')
#end