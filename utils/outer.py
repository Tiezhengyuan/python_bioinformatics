#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 16:12:54 2019

@author: Tiezheng Yuan
"""
import os
import subprocess
#
import altriaseq.utils.basic as ab
#import altriaseq.utils.genome as ag
import altriaseq.utils.biofile as bf

###############################################################################    
class tool:
    def __init__(self, args):
        self.args=args

#
    def run_tool(self,command):
        print("@@@@@@@@@@@@", command)
        output="NA"
        output=subprocess.Popen(command, stdout=subprocess.PIPE, shell=True).stdout.read()  
        #        
        return output

#
    def run_mcscanx(self):
        #combine AA and gfff files in fasta formaat
        self.args.file_faa=self.args.collinear_prefix+'.faa'
        fa_dict=bf.biofile(self.args).read_fa_files(self.args.file_faa)
        #                
        self.args.file_gff=self.args.collinear_prefix+'.gff'
        bf.biofile(self.args).read_gff_files(self.args.file_gff,fa_dict)
        
        #run blastp
        self.run_blastp(self.args.file_faa)
        #run mcscanx
        exe = self.args.dir_altriaseq+'exe/MCScanX'
        command = "{} -s {} {}".format(exe, self.args.min_size, self.args.collinear_prefix)
        self.run_tool(command)

#build control files for mcscanx
    def run_mcscanx_plotter(self):
        file_dot_control = self.args.dir_collinear+'dot_plotter.ctr'
        print(file_dot_control)
        F1 = open(file_dot_control, 'w')
        F1.write('{}\n{}\n'.format(800,800))
        F2 = open(self.args.dir_collinear+'circle_plotter.ctr', 'w')
        F2.write('{}\n'.format(800))
        all_chrs = []
        for name in self.args.genome_names:
            sub_args = self.args.sub_args[name]
            #print(vars(sub_args).keys())
            chrs = [name+'_'+ref for ref in sub_args.ref_list[1:]]
            F1.write('{}\n'.format(','.join(chrs)))
            all_chrs += chrs
        F2.write('{}\n'.format(','.join(all_chrs)))
        F1.close()
        F2.close()
        
        #dot plot
        #Usage: java dot_plotter -g gff_file -s collinearity_file -c control_file -o output_PNG_file
        exe = self.args.dir_altriaseq+"exe/mcscanx_dot_plotter"
        file_gff = self.args.collinear_prefix+'.gff'
        file_collinear = self.args.collinear_prefix+'.collinearity'
        dot_png = self.args.dir_collinear+'synteny_dotplot.png'
        command="{} -g {} -s {} -c {} -o {}".format(exe, file_gff, \
                 file_collinear, file_dot_control, dot_png)
        self.run_tool(command)        
        
        #
    def run_blastp(self, file_ref):
        blast_db=os.path.splitext(file_ref)[0]
        outfile=os.path.splitext(self.args.file_faa)[0]+'.blast'
        if not os.path.isfile(blast_db+'.nhr'):
            print('\n\n###Build a new blast database\n')
            command="makeblastdb -dbtype prot -in {} -out {}".format(file_ref, blast_db)
            self.run_tool(command)
        #
        print('\n\n###Run proein-protein blast:\n')
        exe='blastp -evalue 1e-10 -num_threads 24 -outfmt 6 -num_alignments 5'
        command="{} -query {} -db {} -out {}".format(exe, \
                      self.args.file_faa, blast_db, outfile)
        self.run_tool(command)        
        #


#
    def run_blastn(self):
        if not os.path.isfile(self.args.blast_db+'.nhr'):
            print('\n\n###Build a new blast database\n')
            command="makeblastdb -dbtype nucl -in {} -out {}".format(\
                       self.args.file_fa, self.args.blast_db)
            self.run_tool(command)
        #
        print('\n\n###Run DNA-DNA blast:\n')
        exe='tblastn -evalue 0.001 -num_threads 24 -outfmt 7'
        blast_query='/mnt/rdedata12/yuan/AltriaSeq/bacteria_odb9/ancestral'
        command="{} -query {} -db {} -out {}".format(exe, blast_query,\
                          self.args.blast_db, self.args.blast_out)
        self.run_tool(command) 
        
        
    def run_blat(self, reads):
        fasta_file = self.args.dir_results + '/R2C2_temp_for_BLAT.fasta'
        input_file_fasta = open(fasta_file, 'w')
        for read in reads:
            input_file_fasta.write(">{}\n{}\n".format(read, reads[read][0]))
        input_file_fasta.close()
        #
        outfile=self.args.dir_results+'Splint_to_read_alignments.psl'
        os.system('blat -noHead -stepSize=1 -t=DNA q=DNA -minScore=15 \
              -minIdentity=10 %s %s %s' %(self.args.file_splint, fasta_file, outfile))        


#
    def run_nanoplot(self):
        print('\n\n\n###NanoPlot analyze sequencing quality\n')
        outdir=ab.basic().format_dir(self.args.dir_results+'seqQC')
        exe='NanoPlot -t 24 --minlength 100 --plots hex dot'
        command= "{} --fastq \"{}\" -o \"{}\" ".format(exe, self.args.file_fq, outdir)
        self.run_tool(command)

       
    def run_unicycler(self):
        print('\n\n\n###Unicycler assemble genome sequences\n')
        exe="unicycler --min_fasta_length 100 -t 24 --keep 2 --mode conservative"
        command= "{} -l \"{}\" -o \"{}\" ".format(exe, self.args.file_fq, self.args.dir_genome)
        self.run_tool(command)

    def run_quast(self):
        print('\n\n\n###summarize genome DNA assembly using the tool quast')
        outdir=ab.basic().format_dir(self.args.dir_results+'assessment')
        command= "quast.py -t 24 -o \"{}\" \"{}\" ".format(outdir, self.args.file_fa)
        self.run_tool(command)
        
    def run_prokka(self):
        if self.args.dir_annot.endswith('/'): 
            self.args.dir_annot=self.args.dir_annot[:-1]
        print('\n\n\n###Genome annotation using prokka\n')
        exe="prokka --kingdom Bacteria --force --rfam --cpus 24 --outdir"
        command= "{} \"{}\" \"{}\" ".format(exe, self.args.dir_annot, self.args.file_fa)
        self.run_tool(command)


#
    def run_mauve(self):
        #mauve_exe='/mnt/rdedata12/yuan/tar/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve'
        file_output=self.args.dir_mauve+'output'
        file_tree=self.args.dir_mauve+'guide_tree.txt'
        command="progressiveMauve --output={} --output-guide-tree={} --backbone-output={} {}".format(\
                      file_output, file_tree, self.args.file_backbone,  ' '.join(self.args.fa_files) )
        self.run_tool(command)
       

#run circos
    def run_circos(self):
        shell_command="circos -conf {} -outputdir {}".format(\
                self.args.dir_circos+'circos.conf',self.args.dir_circos)
        print("@@@@@@@@@@@@", shell_command)
        subprocess.Popen(shell_command, stdout=subprocess.PIPE, shell=True).stdout.read()
        


#compress and uncompress files using pigz or unpigz
#the parameter target could be file  or directories
    def run_pigz(self, target, compress=True):
        #
        if compress is True:
            print('Compress the files into *.gz')
            command="pigz {}".format(target)        
        else:
            print('Uncompress them if there are *.gz files in ', self.args.dir_rawdata)
            command="unpigz {}".format(target)
        #
        self.run_tool(command)
        #
#
    #def run_qc(self, rawfile):        
        #qc
        #if rawfile.endswith('.fastq') or rawfile.endswith('.fq'):
        #    shell_command="fastq_quality_filter -v -q 30 -p 80 -i \"{}\" -o \"{}\" ".format(\
        #    rawfile, self.args.dir_rawdata+os.path.basename(rawfile))
        #    #print "@@@@@@@@@@@@", shell_command
        #    subprocess.Popen(shell_command, stdout=subprocess.PIPE, shell=True).stdout.read()            
        
#        
    def run_hisat2(self, sample_name):
        
        #initiate 
        sample_dir=ab.basic().format_dir(self.args.dir_results+sample_name+'/')
        sample_prefix=sample_dir+sample_name
        sample_sam_file=sample_prefix+'.sam'
        novel_ss_file=sample_prefix+'_novel_ss.txt'
        summary_file=sample_prefix+'_summary.log'
        met_file=sample_prefix+'_met.log'
        sample_R1_str=self.args.sample_R1[sample_name]
        #print(sample_R1_str)
        
        #
        hisat2_basic="hisat2 -p 4 -q --dta --phred33 -x {} -S \"{}\"".format(\
                self.args.genome_index, sample_sam_file)
        hisat2_options = "{} --novel-splicesite-outfile {}  --summary-file {} --met-file {}".format(\
                          hisat2_basic, novel_ss_file, summary_file, met_file)
        if self.args.single_end is True:
            command = "{} -U \"{}\" ".format(hisat2_options, sample_R1_str)
        else:
            sample_R2_str=self.args.sample_R2[sample_name]
            #print(sample_R2_str)
            command = "{} -1 \"{}\" -2 \"{}\" ".format(hisat2_options,\
                           sample_R1_str, sample_R2_str)

        #print shell_command
        if not os.path.isfile(sample_sam_file):
            print('\n\n#######Sequence alignment', sample_name)
            self.run_tool(command)
 
        return 

    def run_samtools(self, sample_name, sample_dir=None):
        #initiate 
        if sample_dir is None:
            sample_dir=self.args.dir_results+sample_name
        dir_prefix=ab.basic().format_dir(sample_dir)+sample_name
        sample_sam_file=dir_prefix+'.sam'
        sample_bam_file=dir_prefix+'_sorted.bam'

        #convert sam to bam
        command="samtools sort -@ 24 -o \"{}\" \"{}\"; pigz {} ".format(\
                    sample_bam_file, sample_sam_file, sample_sam_file)
        if not os.path.isfile(sample_bam_file):
            print('Convert sam file into bam file: ', sample_name)
            self.run_tool(command)      
       
        #create RC table
        return 
        
    
    def run_stringtie(self, sample_name):
        #initiate 
        self.args.sample_dir=ab.basic().format_dir(self.args.dir_results+sample_name)
        self.args.sample_prefix=self.args.sample_dir+sample_name
        sample_bam_file=self.args.sample_prefix+'_sorted.bam'
        sample_gtf_file=self.args.sample_prefix+'.gtf'
                
        #assemble transcripts
        if not os.path.isfile(sample_gtf_file):
            command="stringtie -G \"{}\" -o \"{}\" -l {} \"{}\" ".format(\
                        self.args.genome_gtf_file, sample_gtf_file, sample_name, sample_bam_file)
            self.run_tool(command)
        
        #create RC table
        return 

    #merge transcripts
    def run_merge(self):
        #
        gtf_names=[self.args.dir_results+s+'/'+s+'.gtf' for s in self.args.sample_R1.keys()]
        mergelist_file=self.args.dir_results+'mergelist.txt'
        ab.basic().list_to_file(gtf_names, mergelist_file)

        #assemble transcripts
        self.args.merged_gtf_file=self.args.dir_results+'stringtie_merged.gtf'
        command="stringtie --merge -p 24 -G \"{}\" -o \"{}\" \"{}\" ".format(\
             self.args.genome_gtf_file, self.args.merged_gtf_file, mergelist_file)
        self.run_tool(command)
        
        #assemble transcripts
        outdir_prefix=ab.basic().format_dir(self.args.dir_results+'gffcompare')+'gffcmp'
        command="gffcompare -G -r \"{}\" -o {} \"{}\" ".format(\
             self.args.genome_gtf_file, outdir_prefix, self.args.merged_gtf_file)
        self.run_tool(command)      
        #os.system(shell_command)
#
    def run_ballgown(self, sample_name):
        bam_file = self.args.dir_results + sample_name+'/'+sample_name + '_sorted.bam'
        gtf_file = "{}ballgown/{}/{}.gtf".format(self.args.dir_results, sample_name, sample_name)
        if not os.path.isfile(gtf_file):
            command="stringtie -e -B -G \"{}\" -l {} -o \"{}\" \"{}\" ".format(\
                        self.args.genome_gtf_file, sample_name, gtf_file, bam_file)
            self.run_tool(command)

#nanop0ore
    def run_albacore(self):
        print('\n\n\n###Albacore convert fast5 to fastq\n')
        outdir=ab.basic().format_dir(self.args.dir_results+'fastq')
        exe='read_fast5_basecaller.py -t 24 -r -o fastq -f FLO-MIN106 -k SQK-PCS108'
        command= "{} -i \"{}\" -s \"{}\" ".format(exe, self.args.dir_fast5, outdir)
        self.run_tool(command)
        
        print('Concatnate fastq files into one file:',self.args.file_fq)
        dir_pass=outdir+'/workspace/pass/*fastq'
        command="cat {} > {}".format(dir_pass,self.args.file_fq)
        self.run_tool(command)

#nanopore
    def run_pychopper(self):
        barcodes_file='/mnt/rdedata12/yuan/github/pychopper/data/cdna_barcodes.fas'
        report_pdf=self.args.dir_results+'pychopper_report.pdf'
        full_length_fq=self.args.dir_results+'full_length.fq'
        un_fq=self.args.dir_results+'unclassified.fq'
        command="cdna_classifier.py -b {} -r {} -u {} {} {}".format(\
                    barcodes_file,report_pdf, un_fq, self.arg.file_fq, full_length_fq)
        self.run_tool(command)
        return full_length_fq
    
#nanopore
#query_file could fastq or fasta
    def run_minimap2_cdna(self,ref_file=None, bam_file=None):
        if ref_file is None:
            ref_file=self.args.genome_fa_file
        index_file=os.path.splitext(ref_file)[0]+'.mmi'
        #12-1:build index if it doesn't exist
        if not os.path.isfile(index_file):
            exe = "minimap2 -t 24 -k 14 -I 1000G"
            command = "{} -d {} {}".format(exe, index_file, ref_file)
            self.run_tool(command)
            
        #12-2:alignment
        if bam_file is None:
            bam_file=self.args.bam_file
        minimap2 = "minimap2 -t 24 -ax splice -uf"
        samtools = "samtools view -q 10 -F 2304 -Sb | samtools sort -@ 24"
        command="{} {} {} | {} -o {}; samtools index {}".format(minimap2,\
            index_file, self.args.query_file, samtools, bam_file, bam_file)
        self.run_tool(command)
        #
        return 1

        
#nanopore
    def run_pinfish(self):
        #12-3: convert BAM to GFF
        exe = self.args.dir_exe+"spliced_bam2gff -s -t 24"
        raw_gff=self.args.dir_pinfish+'raw_transcripts.gff'
        command = "{} -M {} > {}".format(exe, self.args.bam_file, raw_gff)
        self.run_tool(command)
        
        #12-4:cluster transcripts in GFF
        exe = self.args.dir_exe+"cluster_gff -p 1.0 -t 24 -c 10 -d 10 -e 30"
        cls_tab=self.args.dir_pinfish+'cluster_memberships.tsv'
        cls_gff=self.args.dir_pinfish+'clustered_transcripts.gff'
        command = "{} -a {} {} > {}".format(exe, cls_tab, raw_gff, cls_gff)
        self.run_tool(command)
        
        #12-5: collapse clustered read artifacts
        exe = self.args.dir_exe+"collapse_partials -d 5 -e 30 -f 5000"
        cls_gff_col = self.args.dir_pinfish+'clustered_transcripts_collapsed.gff'
        command = "{} {} > {}".format(exe, cls_gff, cls_gff_col)
        self.run_tool(command)

        #12-6: polish read clusters
        exe = self.args.dir_exe+"polish_clusters -t 24 -c 10"
        pol_trs = self.args.dir_pinfish+"polished_transcripts.fas"
        command="{} -a {} -o {} {}".format(exe, cls_tab, pol_trs, self.args.bam_file)
        self.run_tool(command)
        
        #
        return pol_trs
        
#nanopore
    def run_pinfish_polished(self):
        #12-8: convert BAM of polished transcripts to GFF
        exe = self.args.dir_exe+"spliced_bam2gff -s -t 24"
        gff=self.args.dir_pinfish+'polished_transcripts.gff'
        command = "{} -M {} > {}".format(exe, self.args.bam_file, gff)
        self.run_tool(command)   
        
        #12-9: collapse polished read artifacts
        exe = self.args.dir_exe+"collapse_partials -d 5 -e 30 -f 5000"
        gff_col = self.args.dir_pinfish+'polished_transcripts_collapsed.gff'
        command = "{} {} > {}".format(exe, gff, gff_col)
        self.run_tool(command)
        
        #12:10: Generate corrected transcriptome.
        output_fas = self.args.dir_pinfish+'corrected_transcriptome_polished_collapsed.fas'
        command = "gffread -w {} -g {} {}".format(output_fas, self.args.genome_fa_file, gff_col)
        self.run_tool(command)   
        
        #
        return gff_col
    
#compare gff fils
    def run_gffcompare(self):
        #  pinfish/polished_transcripts_collapsed.gff -o gff_compare
        exe="gffcompare -R -M -C -K"
        outdir=ab.basic().format_dir(self.args.dir_results+'gffcompare')+'gffcompare'
        gff_col = self.args.dir_pinfish+'polished_transcripts_collapsed.gff'
        command = "{} -r {} -o {} {}".format(exe, self.args.genome_gtf_file, outdir, gff_col)
        #print(command)
        self.run_tool(command)          

#nanopore: differential isofomrs 
#nanpore: demultiplexing
    def run_guppy_barcoder(self):
        outdir=ab.basic().format_dir(self.args.dir_results+'fastq')
        exe="guppy_barcoder --barcode_kits SQK-PBK004 -t 24 -q 0"
        command= "{} -i {} -s {}".format(exe, self.args.dir_fastq, outdir)
        self.run_tool(command)
        #
        #fastq files
        R1_files, R2_files, raw_files=ab.basic().seek_fq(outdir)
        barcode_fq=[x for x in raw_files if '/barcode' in x ]
        #print(len(barcode_fq))
        print('Concatnate fastq files into one file:',self.args.file_fq)
        command="cat {} > {}".format(' '.join(barcode_fq),self.args.file_fq)
        self.run_tool(command)                 

#query_file could fastq
    def run_minimap2_quant(self):
        index_file=self.args.transcripts_index+'.mmi'
        #12-1:build index if it doesn't exist
        if not os.path.isfile(index_file):
            command = "minimap2 -t 24 -I 1000G -d {} {}".format(\
                    index_file, self.args.ref_file)
            self.run_tool(command)
            
        #12-2:alignment
        minimap2 = "minimap2 -ax map-ont -t 24 -p 1.0 -N 100 {} {}".format(\
                    index_file, self.args.query_file)
        samtools = "samtools view -Sb | samtools sort -@ 24 -o {}".format(self.args.bam_file)
        command="{} | {} ; samtools index {}".format(minimap2,samtools, self.args.bam_file)
        self.run_tool(command)
        #
        return 1

##
    def run_salmon(self):
        salmon="salmon quant --noErrorModel -p 24 -l A"
        command="{} -t {} -a {} -o {}".format(salmon, self.args.transcripts_fa_file,\
                      self.args.bam_file, self.args.dir_sample)
        self.run_tool(command)        
       
        
######        
#mirdeep2 for miRNA
    def run_mirdeep2(self, sample_name):
        print('\n\n###Analyze {}: {}\n'.format(sample_name, self.args.sample_R1[sample_name]))
        outdir=ab.basic().format_dir(self.args.dir_results+sample_name)
        file_reads=outdir+'reads_collapsed.fa'
        file_arf=outdir+'reads_vs_refdb.arf'
        #mapper
        exe='mapper.pl {} -e -p {} -k {}'.format(self.args.sample_R1[sample_name], 
                       self.args.genome_index, self.args.adapter_3end)
        command = "cd {}; {} -h -i -j -n -l 18 -m -v -o 24 -s {} -t {}".format(
                        outdir, exe, file_reads, file_arf)
        self.run_tool(command)      
        
        #mirdeep2
        file_mature=self.args.dir_genome+'mirbase_mature_nta.fa'
        file_others=self.args.dir_genome+'mirbase_mature_others.fa'
        file_precursor=self.args.dir_genome+'mirbase_hairpin_nta.fa'
        log=outdir+'report.log'
        command = "cd {}; miRDeep2.pl {} {} {} {} {} {} 2>{}".format(outdir, file_reads, \
            self.args.genome_fa_file, file_arf, file_mature, file_others, file_precursor, log)
        self.run_tool(command) 
#end