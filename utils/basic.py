#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 16:11:20 2019

@author: Tiezheng Yuan
"""

import gzip
import os
import re
import multiprocessing as mp   # multi process
import multiprocessing.dummy as mpd #multi-threading

###############################################################################
class basic:
    def __init__(self, args=None):
        self.args=args

#format directory
    def format_dir(self, indir):
        #judge if absolute dir
        if indir.find('/')==0:
            pass
        elif indir.find('~')==0:
            indir = re.sub('~', os.getenv('HOME'), indir)
        elif indir is None:
            indir = os.getcwd()
        else: # ./test or test
            indir = os.path.abspath(indir)
        #alway followed by '/'
        if not indir[-1]=='/':
            indir += '/'
            
        #create directory if not exists
        if not os.path.isdir(indir):
            os.mkdir(indir, 0o777)
        return indir
        
    def print_dict(self, indict):
        n = 1
        for key in sorted(indict.keys()):
            print('{:5}: {:10}\t{}'.format(n, key, indict[key]))
            n += 1

#ouput is read-only file handle
    def readonly_handle(self, infile):
        if infile.endswith('.gz'):
            in_obj = gzip.open(infile, 'rt')
        else:
            in_obj = open(infile, 'rt')
        return in_obj
        
#read certain column in a file
    def to_list(self, infile):
        outlist=[]
        try: 
            in_obj = self.readonly_handle(infile)
            for line in in_obj:
                line = line.strip()
                outlist.append(line)
            in_obj.close()
        except FileNotFoundError:
            print(infile, 'didnot exit!')
            pass
        return outlist

#export list to a text file seperated by return
    def list_to_file(self, inlist, out_file):
        out_obj = open(out_file, 'wt')
        for key in inlist:
            out_obj.write(str(key)+'\n')
        out_obj.close()
        print('write a list to ', out_file)

#
    
#list all files with a given directory and sub directories
#get all files
    def recrusive_files(self, indir): 
        all_files=[]
        for root, dirs, files in os.walk(indir):
            #print('########', root, dirs, len(all_files) )
            for filename in list(files):
                out_file = os.path.join(root, filename)
                if os.path.isfile(out_file) and out_file.find('/.') == -1:
                    all_files.append(out_file)
                    #print(len(all_files), out_file)
        return all_files

#multiple threads by threads pool
    def pp_map_threads(self, func, args_list):
        t = self.args.threads_num
        print("Multiple threads = ", t)
        #t=1
        if t < 2:
            for par in args_list:
                func(par)
        else:
            #threads number
            pool = mpd.Pool(t)
            #pass one argument at a time 
            pool.map(func, args_list) 
            pool.close()
            pool.join()  
            
#multiple process by process pool
#Note: pass a function - and not a method - to the worker processes
    def pp_map_process(self, func, args_list):
        t = self.args.threads_num
        print("Multiple processes: ", t)
        #t=1
        if t < 2:
            for par in args_list:
                func(par)
        else:
            #process number
            pool = mp.Pool(processes=t)
            #pass one argument at a time 
            pool.map(func, args_list)
            pool.close()
            pool.join()  
            
    def seek_fq(self,dir_raw_data):
        print('Retrieve all *.fastq files under', dir_raw_data)
        raw_files = []
        #get all files
        all_files = self.recrusive_files(dir_raw_data) 
        #print(all_files)
        #find file with .fastq or .fq
        for af in all_files:
            m = re.search(r'fastq$|fq$|fastq.gz$|fq.gz$', af)
            if m:
                #print('raw data:',af)
                raw_files.append(af)
        #print(raw_files)
        R2_files=[x for x in raw_files if '_R2' in x]
        R1_files=list(set(raw_files)-set(R2_files))
        #print(R2_files)
        return R1_files, R2_files, raw_files

#read sample_info.csv
    def read_sample_file(self, infile):
        sample_info={}
        #
        try: 
            in_obj = self.readonly_handle(infile)
            for line in in_obj:
                line=line.strip()
                sample_name, name_part, on = line.split(',')
                if on=='yes':
                    sample_info[sample_name]=name_part
            in_obj.close()
        except FileNotFoundError:
            print(infile, 'didnot exit!')
            pass
        #print(sample_info)
        return(sample_info)
