#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 16:30:13 2019

@author: superuser
"""

import os
import sys
import h5py
import statistics
#
import altriaseq.utils.basic as ab


class fast5:
    def __init__(self, args):
        self.args=args
        self.filtered_fast5_files=[]
        self.total_bases = 0
        self.filtered_bases = 0


#
    def decompose_fast5(self):
        #print(self.args.fast5_files)
        
        #filter out by length, qual, and window
        for fast5_file in self.args.fast5_files:
            passes,length=self.check_filters(fast5_file)
            self.total_bases += length
            #print(passes,length)
            if passes:
                self.filtered_fast5_files.append(fast5_file)
                self.filtered_bases += length
        #
        print('Total reads={}, {}bp.'.format(len(self.args.fast5_files), self.total_bases))
        print('Reserved reads={}, {}bp.'.format(len(self.filtered_fast5_files), self.filtered_bases))
        if not self.filtered_fast5_files:
            print("WARNING! No reads are left after filtering processing!\n\n")
            sys.exit()       


        #filter out fast5 based on target_bases
        good_fast5_files=set()
        good_bases = 0
        if self.args.target_bases and self.filtered_bases > self.args.target_bases:
            print('Aim at target bases {} bp'.format(self.args.target_bases))
            quals_lengths =[ self.min_window_qual_and_length(f) for f in self.filtered_fast5_files]
            quals_lengths.sort(reverse=True)
            min_window_qual_threshold =0.0
            for min_window_qual, length, fast5_file in quals_lengths:
                good_bases += length
                min_window_qual_threshold = min_window_qual
                good_fast5_files.add(fast5_file)
                print(fast5_file, self.good_bases, min_window_qual)
                if good_bases > self.args.target_bases:
                    break
            #
            self.filtered_fast5_files=good_fast5_files
            self.filtered_bases=good_bases
            print('Min window quality threshold = ' + '%.2f' % min_window_qual_threshold)
            print('{} reads remain ({}bp)\n'.format(len(good_fast5_files),good_bases) ) 
        else:
            print('WARNING! Skip long reads screening due to the lack of --target_bases or low sequencing depth\n')
            
        #export to fastq
        self.export_fq_fa(self.filtered_fast5_files, self.args.dir_results+'filtered')
        failed_fast5_files=list(set(self.args.fast5_files)-set(self.filtered_fast5_files))
        self.export_fq_fa(failed_fast5_files, self.args.dir_results+'failed')

    def export_fq_fa(self, fast5_files, file_prefix):
        if len(fast5_files)>0:
            print('Export fast5 into the files:{}.fq and {}.fa'.format(file_prefix,file_prefix))
            fq_obj=open(file_prefix+'.fq', 'w')
            fa_obj=open(file_prefix+'.fa', 'w')
            for fast5_file in fast5_files:
                try:
                    hdf5_file = h5py.File(fast5_file, 'r')
                    basecall_location = self.get_fastq_hdf5_location(hdf5_file)
                    if basecall_location:
                        items=hdf5_file[basecall_location].value.decode()
                        #print(hdf5_file[basecall_location].value.decode(), end='', flush=True)
                        fq_obj.write(items)
                        L=items.split("\n")
                        #print(L[4])
                        fa_obj.write(">{}\n{}\n".format(L[0].split(' ')[1], L[1]))
                except IOError:
                    pass
            fq_obj.close()
            fa_obj.close()

    #
    def hdf5_names(self, hdf5_file):
        names = []
        hdf5_file.visit(names.append)
        return names

#
    def check_filters(self, fast5_file):
        tag=(False,0)
        try:
            hdf5_file = h5py.File(fast5_file, 'r')# file handle
            #print(list(hdf5_file.keys()))
            basecall_location = self.get_fastq_hdf5_location(hdf5_file)
            if basecall_location:
                fastq_str = hdf5_file[basecall_location].value
                try:
                    parts = fastq_str.split(b'\n')
                    seq, quals = parts[1], parts[3]
                except IndexError:
                    fastq_str, seq, quals = '', '', ''
                #
                mean_qscore=self.mean_qscore(quals) > self.args.min_mean_qual
                min_length=len(seq) > self.args.min_length
                window_qscore=self.min_window_qscore(quals) > self.args.min_qual_window
                if fastq_str and seq and mean_qscore and min_length and window_qscore :
                    tag=(True, len(seq)) #filter out 
        except (IOError, RuntimeError):
            pass
        return tag


#This function returns the path in the FAST5 file to the best FASTQ. If there are multiple
#basecall locations, it returns the last one (hopefully from the most recent basecalling).
    def get_fastq_hdf5_location(self, hdf5_file):
        names = self.hdf5_names(hdf5_file)
        #print(names)
        basecall_locations = sorted([x for x in names if x.upper().endswith('FASTQ')])
        #print(basecall_locations)
        #two strands
        two_d_locations = [x for x in basecall_locations if 'BASECALLED_2D' in x.upper()]
        #template strand
        template_locations = [x for x in basecall_locations if 'TEMPLATE' in x.upper()]
        #reverse complement strand
        complement_locations = [x for x in basecall_locations if 'COMPLEMENT' in x.upper()]
        #
        if self.args.nano_type == '2D' and two_d_locations:
            return two_d_locations[-1]
        elif self.args.nano_type == 'fwd' and template_locations:
            return template_locations[-1]
        elif self.args.nano_type == 'rev' and complement_locations:
            return complement_locations[-1]
        
        # first choose  2D basecalling
        if two_d_locations:
            return two_d_locations[-1]

        # choose the best based on mean qscore.
        elif template_locations and complement_locations:
            template_location = template_locations[-1]
            complement_location = complement_locations[-1]
            mean_template_qscore = self.mean_score(hdf5_file, template_location)
            mean_complement_qscore = self.mean_score(hdf5_file, complement_location)
            if mean_template_qscore >= mean_complement_qscore:
                return template_location
            else:
                return complement_location

        # If the read has only template basecalling (normal for 1D) or only complement, then that's what we use.
        elif template_locations:
            return template_locations[-1]
        elif complement_locations:
            return complement_locations[-1]

    # If the read has none of the above, but still has a fastq value in its hdf5, that's weird, but
    # we'll consider it a 1d read and use it.
        elif basecall_locations:
            return basecall_locations[-1]

        return None

#Returns the mean qscore over the entire length of the qscore string.
    def mean_qscore(self,quals):
        try:
            return sum([q - 33 for q in quals]) / len(quals)
        except ZeroDivisionError:
            return 0.0

#Returns the minimum mean qscore over a sliding window.
    def min_window_qscore(self,quals):
        quals = [q - 33 for q in quals]  # covert to numbers
        current_window_qscore = statistics.mean(quals[:self.args.window_size])
        shift_count = len(quals) - self.args.window_size
        if shift_count < 1:
            return current_window_qscore
        min_window_qscore = current_window_qscore
        for i in range(shift_count):
            leaving_window = quals[i]
            entering_window = quals[i + self.args.window_size]
            current_window_qscore += (entering_window - leaving_window) / self.args.window_size
            if current_window_qscore < min_window_qscore:
                min_window_qscore = current_window_qscore
        return min_window_qscore

    def min_window_qual_and_length(self, fast5_file):
        try:
            hdf5_file = h5py.File(fast5_file, 'r')
            basecall_location = self.get_fastq_hdf5_location(hdf5_file)
            if basecall_location:
                fastq_str = hdf5_file[basecall_location].value
                try:
                    parts = fastq_str.split(b'\n')
                    seq, quals = parts[1], parts[3]
                    return self.min_window_qscore(quals), len(seq), fast5_file
                except IndexError:
                    pass
        except (IOError, RuntimeError):
            pass
        return 0.0, 0, fast5_file


    #
#
    def find_fast5(self):
        fast5_files=[]
        for dir_name, _, filenames in os.walk(self.args.dir_fast5):
            for filename in filenames:
                if filename.endswith('.fast5'):
                    fast5_files.append(os.path.join(dir_name, filename))
        return fast5_files










#end
