#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 17:31:17 2019

@author: superuser
"""
import re

IN=open('/mnt/rdedata2/Synteny/LargeComparison/N_langsdorfii.gff','r')
OUT=open('/mnt/rdedata2/Synteny/LargeComparison/N_langsdorfii.gtf','w')
for L in IN:
    L=L.rstrip()
    items=L.split('\t')
    annot=items[8]
    if items[2] in ('CDS'):
        #   print(items[2],' ; '.join(new))
        #export CDS
        pattern=re.compile(r'ID (cds.evm.model.NL1.0_scaffold\d+.\d+)')
        transcript_id=pattern.findall(annot)[0]
        items[8]="transcript_id \"{}\"".format(transcript_id)
        #print(annot, items[8])
        OUT.write('\t'.join(items)+"\n")
        
        #exon row
        if items[6]=='-':
            exon=re.compile(r'(\d+)..(\d+)\): (\d+)')
            strand='-'
        else:
            exon=re.compile(r'(\d+)..(\d+): (\d+)') 
            strand='+'
        #print(annot,exon.findall(annot), '\n')
        for exon_pos in exon.findall(annot):
            items[2]='exon'
            items[3], items[4], items[7] = exon_pos
            #print(items)
            #new=[items[0],items[1],'exon',start,end,'.',strand,copy,items[8]]
            OUT.write('\t'.join(items)+"\n")
                    
IN.close()
OUT.close()    
    
    
#

