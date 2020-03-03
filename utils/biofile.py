#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 18:29:24 2019

@author: superuser
"""

#import itertools
import os
import re
#
import altriaseq.utils.basic as ab

###############################################################################    
class biofile:
    def __init__(self, args=None):
        self.args=args
        self.record_num=0


#read fasta used for karyotype.txt of circos       
    def read_fa(self, file_fa=None):
        infile = self.args.file_fa if file_fa is None else file_fa
        print('elicit dna information from ', infile)
        #store accession: sequences
        ref_dict = {'chr':''}
        ref_list = ['chr']
        in_obj=ab.basic().readonly_handle(infile)
        #read fa file
        for line in in_obj:
            line = line.rstrip()
            if line.find(">") == 0:
                #fa_id = self.acc_number(line)
                fa_id = re.sub(r"^>", "", line.split(" ")[0])
                ref_list.append(fa_id)
                ref_dict[fa_id] = []
                self.record_num += 1
            else:
                ref_dict[fa_id].append(line)
        in_obj.close()
        
        #arrage ref
        total_len=0
        pos_start=0
        pos_end=-1
        ref_start={}
        ref_end={}
        for fa_id in ref_list[1:]:
            #print fa_id, ref_dict[fa_id][:3], "\n"
            ref_dict[fa_id] = ''.join(ref_dict[fa_id])
            ref_dict['chr'] += ref_dict[fa_id]
            #
            pos_start = pos_end+1
            ref_start[fa_id]=pos_start
            pos_end = total_len + len(ref_dict[fa_id])
            ref_end[fa_id]=pos_end
            #print(ref_start[fa_id], ref_end[fa_id], len(ref_dict[fa_id]))
            
            #update for the next
            total_len += len(ref_dict[fa_id])
            
            #print fa_id, len(re.findall('^N*', ref_dict[fa_id])[0]), len(ref_dict[fa_id]),'\n'
        ref_start['chr']=0
        ref_end['chr']=total_len
        #print(ref_end['chr'])
        return ref_dict, ref_list, ref_start,ref_end
        
#elicit seq and acc from fasta file
    def fa_dict(self, file_fa):
        n=0
        #store accession: sequences
        ref_dict = {}
        ref_list = []
        in_obj=ab.basic().readonly_handle(file_fa)
        #read fa file
        for line in in_obj:
            line = line.rstrip()
            if line.find(">") == 0:
                #acc = self.acc_number(line)
                acc = re.sub(r"^>", "", line.split(" ")[0])
                ref_list.append(acc)
                ref_dict[acc] = []
                n+=1
            else:
                ref_dict[acc].append(line)
        in_obj.close()
        print('Number of sequences:', n)
        return ref_dict, ref_list
    
#outfile is combine them into one fasta file 
#used by mcscanx
    def read_fa_files(self, outfile):
        fa_dict={}
        fa_obj=open(outfile, 'w')
        for name in self.args.genome_names:
            sub_args=self.args.sub_args[name]
            #read AA
            ref_dict, ref_list=self.fa_dict(sub_args.file_faa)
            #print('\n\n\n%%%', ref_list)
            #export
            for ID in ref_list:
                key=name+'_'+ID
                #if key in aa_gene:
                #    print(aa_gene[key]['genome'], key, name)
                fa_dict[key]={'geneID':ID, 'aa_seq':''.join(ref_dict[ID]), 'gff':None}
                fa_obj.write(">{}\n{}\n".format(key, ''.join(ref_dict[ID])))
        fa_obj.close()
        return fa_dict

#outfile is combine them into one gff file
#elicit from fasta and gff
#used for mcscanx
    def read_gff_files(self, outfile, fa_dict):
        out_obj=open(outfile, 'w')  
        for name in self.args.genome_names:
            #read gene position
            in_obj=open(self.args.sub_args[name].file_gff, 'r')
            for line in in_obj:
                line = line.rstrip()
                items=line.split('\t')
                if len(items)==9:#remove comments line
                    annot=re.split(';| ; ', items[8]) #column #9
                    for one in annot:
                        #tage could be geneID or Accession
                        start=re.search('=| ', one).start()
                        tag_name,tag=one[:start], one[(start+1):]
                        ID=name+'_'+tag
                        #print('##{}##{}##{}##'.format(one, tag_name, tag))
                        if tag_name == self.args.gff_tag_name and ID in fa_dict.keys():
                            if fa_dict[ID]['gff'] is None:#unique line
                                out=[name+'_'+items[0], ID, items[3],items[4]]
                                fa_dict[ID]['gff']=out
                                out_obj.write("{}\n".format('\t'.join(out)))
                                break
            in_obj.close()
                    

        out_obj.close()
        return fa_dict        
        
        
# Takes a FASTQ file and returns dictionary of lists
#    readDict {'name_root':['full_header', seq, quality]...}
    def read_fq(self):
        readDict = {}
        lineNum, lastPlus, lastHead, skip = 0, False, '', False
        for line in open(self.par['file_fq']):
            line = line.rstrip()
            if not line:
                continue

            if lineNum % 4 == 0 and line[0] == '@':
                name = line[1:].split()[0]
                readDict[name], lastHead = [], name

            if lineNum % 4 == 1:
                readDict[lastHead].append(line)

            if lineNum % 4 == 2:
                lastPlus = True

            if lineNum % 4 == 3 and lastPlus:
                avgQ = sum([ord(x)-33 for x in line])/len(line)
                sLen = len(readDict[lastHead][-1])
                if avgQ >= self.par['quality_cutoff'] and sLen >= self.par['read_length_cutoff']:
                    readDict[lastHead].append(line)
                    readDict[lastHead] = tuple(readDict[lastHead])
                else:
                    del readDict[lastHead]
                lastPlus, lastHead = False, ''

            lineNum += 1
        return readDict

    def write_fq(self, adapter_dict, reads):
        success = 0
        os.system('mkdir ' + self.par['dir_results']  + '/splint_reads')
        for read in reads:
            name, sequence, quality = read, reads[read][0], reads[read][1]
            adapter_plus = sorted(adapter_dict[name]['+'],key=lambda x: x[1], reverse=True)
            adapter_minus=sorted(adapter_dict[name]['-'], key=lambda x: x[1], reverse=True)
            plus_list_name, plus_list_position = [], []
            minus_list_name, minus_list_position = [], []

            for adapter in adapter_plus:
                if adapter[0] != '-':
                    plus_list_name.append(adapter[0])
                    plus_list_position.append(adapter[2])
            for adapter in adapter_minus:
                if adapter[0] != '-':
                    minus_list_name.append(adapter[0])
                    minus_list_position.append(adapter[2])

            if len(plus_list_name) > 0 or len(minus_list_name) > 0:
                success += 1
                splint_file=self.par['dir_results'] + 'splint_reads/' + str(int(success/4000)) + '/R2C2_raw_reads.fastq'
                try:
                    out_fastq = open(splint_file, 'a')
                except:
                    os.system('mkdir ' + self.par['dir_results']  + '/splint_reads/' + str(int(success/4000)))
                    out_fastq = open(splint_file, 'w')
                    list_pos=  str(plus_list_position[0]) if len(plus_list_name) > 0 else str(minus_list_position[0])
                    out_fastq.write('@' + name + '_' + list_pos + '\n' + sequence + '\n+\n' + quality + '\n')
            else:
                no_splint_file=self.par['dir_results'] + 'splint_reads/No_splint_reads.fastq'
                try:
                    out_fastq = open(no_splint_file, 'a')
                except:
                    out_fastq = open(no_splint_file, 'w')
                    out_fastq.write('>' + name + '\n' + sequence + '\n+\n' + quality + '\n')

#
    def file_to_dict(self, infile, pattern=","):
        outdict={}
        IN=open(infile, 'r')
        for line in IN:
            line = line.rstrip()
            if line.startswith('#'):
                continue
            else:
                items=line.split(pattern)
                outdirct[items[0]]=items[1]
        IN.close()
        return outdirct

#demultiplexing
    def demultiplexing(self, indir):
        pass
        
        

        
#end