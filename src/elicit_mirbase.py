#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 15:12:46 2019

@author: superuser
"""
import argparse
import os
import sys

#
import altriaseq.utils.biofile as bf
#import altriaseq.utils.basic as ab


####
def split_mirbase(fa_dict, specie, out_prefix):
    OUT1=open(out_prefix+specie+'.fa', 'w')
    OUT2=open(out_prefix+'others.fa', 'w')
    m=0
    n=0
    for mir_id, mir_seq in fa_dict.items():
        seq=''.join(mir_seq).replace('U', 'T')
        outline=">{}\n{}\n".format(mir_id, seq)
        if mir_id.startswith(specie):
            #print(mir_id)
            OUT1.write(outline)
            m += 1
        else:
            OUT2.write(outline)
            n += 1
    OUT1.close()
    OUT2.close()
    print(m,n)
#####
fa_mature='/mnt/rdedata12/yuan/mirbase.org/pub/mirbase/CURRENT/mature.fa'
fa_hairpin='/mnt/rdedata12/yuan/mirbase.org/pub/mirbase/CURRENT/hairpin.fa'
dir_out='/mnt/rdedata12/yuan/AltriaSeq/NT3.1/'

#mature
fa_dict, ref_list=bf.biofile().fa_dict(fa_mature)
split_mirbase(fa_dict, 'nta', dir_out+'mirbase_mature_')

#hairpin
fa_dict, ref_list=bf.biofile().fa_dict(fa_hairpin)
split_mirbase(fa_dict, 'nta', dir_out+'mirbase_hairpin_')


#end