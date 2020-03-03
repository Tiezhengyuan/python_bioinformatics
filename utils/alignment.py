#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 21:53:24 2019

@author: superuser
"""

import itertools
import os
import re
#
import altriaseq.utils.basic as ab
###############################################################################    
class alignment:
    def __init__(self, par):
        self.par=par

#
    def parse_blat(self):
        alignment_file = self.par['dir_results'] + '/Splint_to_read_alignments.psl'
        fasta_file =  self.par['dir_results'] + '/R2C2_temp_for_BLAT.fasta'
        adapter_dict = {}
        length = 0
        for line in open(fasta_file):
            length += 1
        iterator = 0
        infile = open(fasta_file, 'r')
        while iterator < length:
            line = infile.readline()
            sequence = infile.readline()
            name = line[1:].strip()
            adapter_dict[name] = {}
            adapter_dict[name]['+'] = []
            adapter_dict[name]['-'] = []
            adapter_dict[name]['+'].append(('-', 1, 0))
            adapter_dict[name]['-'].append(('-', 1, len(sequence)))
            iterator += 2

        burn_dict = {}
        for line in open(alignment_file):
            a = line.strip().split('\t')
            read_name, adapter, strand = a[9], a[13], a[8]
            gaps, score = float(a[5]), float(a[0])
            sequence_length = int(a[10])
            if gaps < 50 and score > 50: # Looks for unspliced quality alignment
                if strand == '+':
                    start = int(a[11]) - int(a[15])
                    end = int(a[12]) + int(a[14]) - int(a[16])
                if strand == '-':
                    start = int(a[11]) - (int(a[14]) - int(a[16]))
                    end = int(a[12]) + int(a[15])
                position = min(max(0, int(start+((end-start)/2))), sequence_length-1)
                adapter_dict[read_name][strand].append((adapter, float(a[0]), position))
        return adapter_dict

#