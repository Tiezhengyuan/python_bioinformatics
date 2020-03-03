#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 17:14:55 2019

@author: Tiezheng Yuan
"""
import itertools
import os
import re
import pandas as pd
#
import altriaseq.utils.basic as ab
import altriaseq.utils.biofile as bf

###############################################################################    
class genome:
    def __init__(self, args):
        self.args=args
        self.record_num = 0

        
# elicit annotations from DNA fasta files, AA fasta files, gff files
#type in (DNA, GFF, AA): args.file_fa, args.file_gff, args.file_faa
#args.dir_annot
    def elicit_annotation(self):
        for name in self.args.sub_args.keys():
            print (name)
            sub_args=self.args.sub_args[name] #namespace
            #read fasta
            if sub_args.file_fa:
                ref_dict,ref_list, ref_start, ref_end=bf.biofile(sub_args).read_fa()
                self.args.sub_args[name].ref_dict=ref_dict
                self.args.sub_args[name].ref_list=ref_list
                self.args.sub_args[name].ref_start=ref_start
                self.args.sub_args[name].ref_end=ref_end
                print(self.args.sub_args[name].ref_list)
                
            #read gff
            #print(vars(sub_args).keys())
            if 'file_gff' in vars(sub_args).keys():
                pass
            else:
                files=[ x for x in os.listdir(sub_args.dir_annot) if x.endswith('.gff')]
                if len(files)>0:
                    self.args.sub_args[name].file_gff=sub_args.dir_annot+files[0]
                    #print('%%%', self.args.sub_args[name].file_gff)

            #read gff
            #print(vars(sub_args).keys())
            if 'file_faa' in vars().keys():
                pass
            else:
                files=[ x for x in os.listdir(sub_args.dir_annot) if x.endswith('.faa')]
                if len(files)>0:
                    self.args.sub_args[name].file_faa=sub_args.dir_annot+files[0]
                    #print('%%%', self.args.sub_args[name].file_gff)
            #
            if os.path.isfile(sub_args.dir_annot+'karyotype.txt'):
                self.args.sub_args[name].file_karyotype=sub_args.dir_annot+'karyotype.txt'
                
        return self.args.sub_args
            

#
    def read_gff(self):
        gff_dict={}
        file_collinear=self.args.collinear_prefix+'.gff'
        in_obj=open(file_collinear, 'r')
        for line in in_obj:
            line=line.rstrip()
            contig_name, geneID, start,end=line.split('\t')
            if contig_name not in gff_dict.keys():
                gff_dict[contig_name]={}
            gff_dict[contig_name][geneID]=(int(start),int(end))
        in_obj.close()
        #print([(x, len(gff_dict[x])) for x in gff_dict.keys()])
        return gff_dict
#
    def read_collinearity(self):
        tag=0
        align, contig_A,contig_B= 0, '', ''
        collinear={}
        file_collinear=self.args.collinear_prefix+'.collinearity'
        in_obj=open(file_collinear, 'r')
        for line in in_obj:
            line=line.rstrip()
            if 'Alignment' in line:
                #if align in collinear.keys():
                #    print(align, collinear[align]['A'], collinear[align]['B'],\
                #          len(collinear[align]['pairs']))
                #print(line)
                n=0
                titles=line.split(' ')
                contig_A,contig_B=titles[6].split('&')
                align=titles[2]
                #initiate
                collinear[align]={'A':contig_A, 'B':contig_B, 'strand':titles[7], 'pairs':[]}
                tag=1
            elif tag==1:
                items=line.split('\t')
                collinear[align]['pairs'].append((items[1],items[2]))
                n+=1
                #print(items)
                
        in_obj.close()
        return collinear
        
#sample names
    def sample_info(self):
        #fastq files
        R1_files, R2_files, raw_files=ab.basic().seek_fq(self.args.dir_rawdata)

        #read sample names
        sample_info=ab.basic().read_sample_file(self.args.file_samples)
        sample_names=sample_info.keys()
        #ab.basic().print_dict(sample_info)

        #sample_info                            
        sample_R1={}
        sample_R2={}
        #connect raw file to sample name
        for sample_name in sample_names:
            name_part=sample_info[sample_name]
            R1_str=','.join([ x for x in R1_files if name_part in x ])
            #print(name_part, R1_str)
            sample_R1[sample_name]=R1_str
            if hasattr(self.args, "paired") and self.args.paired is True:
                sample_R2[sample_name]=R1_str.replace('_R1','_R2')
        #
        #ab.basic().print_dict(sample_R1)
        #ab.basic().print_dict(sample_R2)
        return sample_R1, sample_R2
        

#CDS
    def CDS_file(self):
        #export
        file_sense=self.args.dir_circos+'highlight_senseCDS.txt'
        f1_obj = open(file_sense, 'wt')
        file_antisense=self.args.dir_circos+'highlight_antisenseCDS.txt'
        f2_obj = open(file_antisense, 'wt')
        print('Save into ', file_sense, file_antisense)
        #
        in_obj=ab.basic().readonly_handle(self.args.file_gff)
        #read fa file
        for line in in_obj:
            if '\tCDS\t' in line:
                items=line.split('\t')
                contig_id=items[0]
                #print(contig_id)
                start=self.args.ref_start[contig_id]+int(items[3])-1
                end=self.args.ref_start[contig_id]+int(items[4])-1
                cds='{} {} {} black'.format(self.args.chr_name, start,end)
                if items[6]=='+':
                    f1_obj.write(cds+'\n')
                    #print(cds)
                elif items[6]=='-':
                    f2_obj.write(cds+'\n')
        f1_obj.close()
        f2_obj.close()    
        in_obj.close()

#GC_content
#require DNA fasta
    def GC_file(self):
        seq=self.args.ref_dict['chr']
        seq_len=len(self.args.ref_dict['chr'])
        #print(seq_len)
        GC_total=len([x for x in seq if x=='G' or x=='C'])/seq_len
        #export
        file_f1=self.args.dir_circos+'plot_GCcontent.txt'
        f1_obj = open(file_f1, 'wt')
        file_f2=self.args.dir_circos+'plot_GCskew.txt'
        f2_obj = open(file_f2, 'wt')
        start=0
        bin_size=int(seq_len/10000)#default 10000 bins
        while start<seq_len:
            end=start+bin_size if start+bin_size<seq_len else seq_len
            sub_seq=seq[start:end]
            GC_perc=len([x for x in sub_seq if x=='G' or x=='C'])/bin_size
            GC="{} {} {} {}\n".format(self.args.chr_name, start,end,round(GC_perc,4))
            f1_obj.write(GC)
            skew="{} {} {} {}\n".format(self.args.chr_name, start,end,round(GC_perc-GC_total,4))
            f2_obj.write(skew)
            start=end
            #print(GC)
        f1_obj.close()    
        f2_obj.close()
#
    def read_backbone_file(self):
        #return combinations
        comb_index=[x for x in itertools.combinations(range(len(self.args.genome_names)), 2)]
        print(comb_index)
        #
        print('@@@Read ',self.args.file_backbone)
        in_obj=ab.basic().readonly_handle(self.args.file_backbone)
        first_line=in_obj.readline()
        #read fa file
        F1 = open(self.args.dir_circos+'link_forward.txt', 'wt')
        F2 = open(self.args.dir_circos+'link_reverse.txt', 'wt')
        for line in in_obj:
            line=line.strip()
            pos=[int(x) for x in line.split('\t')]
            #
            for a_index,b_index in comb_index:
                a_name=self.args.genome_names[int(a_index)]
                b_name=self.args.genome_names[int(b_index)]
                a_start=pos[int(a_index)*2]
                a_end=pos[int(a_index)*2+1]
                b_start=pos[int(b_index)*2]
                b_end=pos[int(b_index)*2+1]
                #
                if a_end>0 and b_end>0:
                    F1.write("{} {} {} {} {} {}\n".format(\
                    a_name,a_start,a_end,b_name,b_start,b_end))
                elif a_end<0 and b_end<0:
                    F1.write("{} {} {} {} {} {}\n".format(\
                    a_name,abs(a_start),abs(a_end),b_name,abs(b_start),abs(b_end)))
                elif a_end*b_end<0:
                    pos=[abs(x) for x in pos]
                    F2.write("{} {} {} {} {} {}\n".format(\
                    a_name,abs(a_start),abs(a_end),b_name,abs(b_start),abs(b_end)))
        in_obj.close()
        F1.close()
        F2.close()

#merge mirdeep2 files into one RC table
    def merge_mirdeep2(self):
        mir_dict={} #mature miRNA
        pre_dict={} #precursor miRNA
        #
        for sample_name in self.args.sample_names:
            sample_dir=self.args.dir_results+sample_name
            all_files=ab.basic(self.args).recrusive_files(sample_dir)
            mir_file=[ x for x in all_files if 'miRNAs_expressed_all' in x][0]
            if os.path.isfile(mir_file):
                print('%%%%%%%', mir_file)
                #read 
                indf=pd.read_csv(mir_file, sep='\t')
                indf.columns=['mature','RC','precursor', 'total','seq','seq_norm']
                #mature
                mature_df=indf.groupby('mature').agg('mean')
                mir_dict[sample_name]=dict(mature_df['RC'])
                #print(sample_name, len(mature_df),':', len(mir_dict[sample_name]))
                #precursor
                precursor_df=indf.groupby('precursor').agg('mean')
                pre_dict[sample_name]=dict(precursor_df['RC'])
                #print(sample_name, len(precursor_df),':', len(pre_dict[sample_name]))

            else:
                print('Warning! No mirRNA detected in ', sample_dir)
        #export
        outdir=ab.basic().format_dir(self.args.dir_results+'stat')
        #
        mir_df=pd.DataFrame.from_dict(mir_dict, orient='index')
        mir_df=mir_df.transpose().sort_index()
        #print(mir_df)
        mir_df.to_csv(outdir+'matured_miRNA_RC.txt', sep='\t',index=True, \
                      header=True, index_label='matured_miRNA')
        #
        pre_df=pd.DataFrame.from_dict(pre_dict, orient='index')
        pre_df=pre_df.transpose().sort_index()
        #print(pre_df)
        pre_df.to_csv(outdir+'precursor_miRNA_RC.txt', sep='\t',index=True,\
                      header=True, index_label='precursor_miRNA')        

#merge quant.sf files into one RC table determined by the tool salmon of nanopore
    def merge_salmon(self):
        TPM_dict={} #scaling RC (RC/total-RC)
        RC_dict={} #read counts
        #
        for sample_name in self.args.sample_names:
            sample_dir=self.args.dir_results+sample_name
            infile=sample_dir+'/quant.sf'
            if os.path.isfile(infile):
                print('%%%%%%%', infile)
                #read 
                indf=pd.read_csv(infile, sep='\t', index_col=0)
                indf.columns=['length','effective_length', 'TPM','RC']
                #ids=list(indf['transcripts_id'])
                #print(indf)
                #
                TPM_dict[sample_name]=dict(indf['TPM'])
                RC_dict[sample_name]=dict(indf['RC'])
 
            else:
                print('Warning! No mirRNA detected in ', sample_dir)
        #export
        outdir=ab.basic().format_dir(self.args.dir_results+'stat')
        #
        #print(RC_dict)
        RC_df=pd.DataFrame.from_dict(RC_dict)
        RC_df=RC_df.fillna(0)
        #print(RC_df)
        #print(RC_df.shape)
        #RC_df.index=ids
        RC_df.to_csv(outdir+'transcripts_RC.txt', sep='\t',index=True, \
                      header=True, index_label='transcripts_id')
        #
        TPM_df=pd.DataFrame.from_dict(TPM_dict)
        TPM_df=TPM_df.fillna(0)
        #print(TPM_df)
        TPM_df.to_csv(outdir+'transcripts_TPM.txt', sep='\t',index=True, \
                      header=True, index_label='transcripts_id')
       

def read_novel(self):
    un=pd.read_csv('/mnt/rdedata4/miRNA/results2/L42-1/result_11_05_2019_t_18_29_24.csv', \
               sep='\t',header=None, names=list(range(17)))

    a=list(un[un[0].str.contains('novel miRNAs')].index)[0]
    b=list(un[un[0].str.contains('mature miRBase')].index)[0]
    novel=un.iloc[(a+2):(b-1),:]
    novel.columns=list(un.iloc[(a+1),:])
    novel
    novel.dtypes
    novel['miRDeep2 score']=novel['miRDeep2 score'].astype(float)
    novel['mature read count']=novel['mature read count'].astype(int)
    novel['significant randfold p-value']=novel['significant randfold p-value'].astype('str')
    
    novel=novel[novel['miRDeep2 score']<=10]
    novel=novel[novel['miRDeep2 score']>=-10]
    novel=novel[novel['significant randfold p-value']=='yes']
    novel
    novel.columns
    novel['provisional id']
    novel['mature read count']
    novel['example miRBase miRNA with the same seed']
    novel['precursor coordinate']

#
#end