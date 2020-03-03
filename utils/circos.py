#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 16:12:57 2019

@author: Tiezheng Yuan
"""
#
import subprocess
import os
#
import altriaseq.utils.biofile as bf
import altriaseq.utils.genome as ag
import altriaseq.utils.outer as ao


#############################################################################
#
class circos:
    def __init__(self,args):
        self.args=args
        self.lines=[]
#
    def collinear(self):
        self.args.dir_circos = self.args.dir_collinear
        #elicit collinear pairs
        collinear_dict=ag.genome(self.args).read_collinearity()
        
        #elicit contigs position
        contig_dict={}
        for name in self.args.genome_names:
            sub_args=self.args.sub_args[name]
            #print(vars(sub_args).keys())
            for ref in sub_args.ref_list:
                contig_dict[name+'_'+ref]={'genome':name, 'start':int(sub_args.ref_start[ref]),
                'end':int(sub_args.ref_end[ref])}
        #print(contig_dict) 
        
      
        #elicit gene position
        gff_dict=ag.genome(self.args).read_gff()
        #print(gff_dict.keys())
        
        #build link file
        f1_obj=open(self.args.dir_collinear+'link_forward.txt', 'wt')
        f2_obj=open(self.args.dir_collinear+'link_reverse.txt', 'wt')
        for key in collinear_dict:
            pairs=collinear_dict[key]['pairs']
            contig_A=collinear_dict[key]['A']
            genome_A=contig_dict[contig_A]['genome']
            contig_B=collinear_dict[key]['B']
            genome_B=contig_dict[contig_B]['genome']
            if not genome_A==genome_B:
                for gene_A, gene_B in pairs:
                    #print(genome_A, genome_B)
                    A_start,A_end=gff_dict[contig_A][gene_A]
                    A_start += contig_dict[contig_A]['start']
                    A_end += contig_dict[contig_A]['start']
                    B_start,B_end=gff_dict[contig_B][gene_B]
                    B_start += contig_dict[contig_B]['start']
                    B_end += contig_dict[contig_B]['start']
                    out=[genome_A,A_start,A_end,genome_B,B_start,B_end]
                    out=' '.join([str(x) for x in out])
                    #print(collinear_dict[key]['strand'])
                    if collinear_dict[key]['strand']=='plus':
                        f1_obj.write("{}\n".format(out))
                    else:
                        f2_obj.write("{}\n".format(out))
        f1_obj.close()
        f2_obj.close()
        #
        #build data files
        #karyotype.txt
        karyotype_files=[self.args.sub_args[x].file_karyotype for x in self.args.genome_names]
        self.args.dir_circos=self.args.dir_collinear
        outfile=self.args.dir_circos+'karyotype.txt'
        shell_command="cat {} > {} ".format(" ".join(karyotype_files),outfile)
        os.system(shell_command)
        #        #build *.conf
        self.export_conf(self.conf_circos(['plots']),'circos')
        self.export_conf(self.conf_ideogram(),'ideogram')
        self.export_conf(self.conf_ticks(),'ticks')
        self.export_conf(self.conf_plots_collinear(),'plots')
        #run circos
        ao.tool(self.args).run_circos()

               
#*.fa and *.gff files are available          
    def one_genome(self):
        print('\n### Visualize genome using Circos\n', self.args.chr_name)
        #retrieve from file_fa, file_gff, file_faa
        self.elicit_annotation()
        #print('###', self.args.file_fa)
        #print('###', self.args.file_gff)
        
        #build data files
        #karotype.txt
        self.karotype_file()
        #CDS file
        ag.genome(self.args).CDS_file()
        #GC file
        ag.genome(self.args).GC_file()
        
        #build *.conf
        self.export_conf(self.conf_circos(['highlights','plots']),'circos')
        self.export_conf(self.conf_ideogram(),'ideogram')
        self.export_conf(self.conf_ticks(),'ticks')
        self.export_conf(self.conf_plots(),'plots')
        self.export_conf(self.conf_highlights(),'highlights')
        
        #run circos
        ao.tool(self.args).run_circos()
        
#only *.fa available
    def simple_genome(self):
        print('\n\n\n### Visualize genome using Circos\n')
        #build data files
        #
        self.karotype_file()
        #
        ag.genome(self.args).GC_file()
        
        #build *.conf
        self.export_conf(self.conf_circos(),'circos')
        self.export_conf(self.conf_ideogram(),'ideogram')
        self.export_conf(self.conf_ticks(),'ticks')
        
        #run circos
        ao.tool(self.args).run_circos()


#        
    def multiple_genome(self):
        #combine data files
        data=['highlight_antisenseCDS', 'highlight_senseCDS',
              'karyotype', 'plot_GCcontent', 'plot_GCskew']
        data=[x+'.txt' for x in data]
        for file_name in data:
            file_str=" ".join([x+file_name for x in self.args.in_circos])
            outfile=self.args.dir_circos+file_name
            shell_command="cat {} > {} ".format(file_str,outfile)
            #print("@@@@@@@@@@@@", shell_command)
            subprocess.Popen(shell_command, stdout=subprocess.PIPE, shell=True).stdout.read()
        
        #copy conf files
        #build *.conf
        self.export_conf(self.conf_circos(['highlights','plots']),'circos')
        self.export_conf(self.conf_ideogram(),'ideogram')
        self.export_conf(self.conf_ticks(),'ticks')
        self.export_conf(self.conf_plots(),'plots')
        self.export_conf(self.conf_highlights(),'highlights')
        #update plots.conf
        self.append_conf(self.conf_links(),'plots')
        #buid links.txt
        ag.genome(self.args).read_backbone_file()
        
        #run circos
        ao.tool(self.args).run_circos()

# elicit annotations from DNA fasta files, AA fasta files, gff files
#type in (DNA, GFF, AA): args.file_fa, args.file_gff, args.file_faa
#args.dir_annot
    def elicit_annotation(self):
        #read fasta
        if self.args.file_fa:
            ref_dict,ref_list, ref_start, ref_end=bf.biofile(self.args).read_fa()
            self.args.ref_dict=ref_dict
            self.args.ref_list=ref_list
            self.args.ref_start=ref_start
            self.args.ref_end=ref_end
            #print(self.args.ref_list)
                
            #read gff
            #print(vars(self.args).keys())
            if 'file_gff' in vars(self.args).keys():
                pass
            else:
                files=[ x for x in os.listdir(self.args.dir_annot) if x.endswith('.gff')]
                if len(files)>0:
                    self.args.file_gff=self.args.dir_annot+files[0]
                    #print('%%%', self.args.file_gff)

            #read gff
            if 'file_faa' in vars(self.args).keys():
                pass
            else:
                files=[ x for x in os.listdir(self.args.dir_annot) if x.endswith('.faa')]
                if len(files)>0:
                    self.args.file_faa=self.args.dir_annot+files[0]
                    #print('%%%', self.args.file_gff)
                
        return self.args
            
#karotype file
    def karotype_file(self):
        #
        out_file=self.args.dir_circos+'karyotype.txt'
        out_obj = open(out_file, 'wt')
        print('Save into ', out_file)
        #chromosome 1
        out_obj.write('chr - {} {} {} {} black\n'.format(self.args.chr_name, \
            self.args.chr_name, self.args.ref_start['chr'],self.args.ref_end['chr']))
        #band
        col='white'
        for contig in self.args.ref_list[1:]:
            pos="{} {} {} {}".format('Band'+contig, 'Contig'+contig, \
                 self.args.ref_start[contig],self.args.ref_end[contig])
            #export contigs as band
            out_obj.write("band {} {} {}\n".format(self.args.chr_name, pos,col))
            #
            col='black' if col == 'white' else 'white'
            #print(pos)
        out_obj.close()

#
    def insert_tag(self, tag_name):
        self.lines.insert(0, "<{}>".format(tag_name) )
        self.lines.append("</{}>".format(tag_name))
        self.lines.append("") #add white line
    def insert_tag2(self, lines, tag_name):
        lines.insert(0, "<{}>".format(tag_name) )
        lines.append("</{}>".format(tag_name))
        lines.append("") #add white line        
        return lines
    
#conf_list=['highlights','plots']
    def conf_circos(self, conf_list=None):
        #karyotype
        line1=['karyotype = karyotype.txt',
               'chromosomes_units = 1000000',
               'chromosomes_display_defaults = yes','',
               #personal confs
               '<<include ideogram.conf>>', '<<include ticks.conf>>']
        line3=[# The remaining content is standard and required. It is imported
               # from default files in the Circos distribution.
               '<image>', 
               ' <<include /usr/local/bin/circos/etc/image.conf>>',
               '</image>','',
               # RGB/HSV color definitions, color lists, location of fonts, fill patterns.
               '<<include /usr/local/bin/circos/etc/colors_fonts_patterns.conf>>',
               # Debugging, I/O an dother system parameters
               '<<include /usr/local/bin/circos/etc/housekeeping.conf>>','']
        if conf_list is None:
            self.lines=line1+line3
        else:
            line2=[ '<<include '+x+'.conf>>' for x in conf_list]
            self.lines=line1+line2+line3
        return self.lines

    #
    def conf_ideogram(self):    
        self.lines=[' <spacing>','  default = 0.01r','  break = 0.25r',' </spacing>','',
                ' thickness = 80p', ' radius = 0.85r', '',
                ' show_label = yes', ' label_font = bold', ' label_with_tag = yes',
                ' label_radius = dims(ideogram,radius) + 0.05r',
                ' label_size = 48',' label_parallel = yes', ' label_case = upper',
                '', ' stroke_thickness = 3', ' stroke_color = black',' fill = yes',
                '',' show_bands = yes',' fill_bands = yes']
        self.insert_tag('ideogram')
        return self.lines
    
    def conf_ticks(self):
        lines_show=['show_ticks = yes','show_tick_labels = yes','show_grid = yes', '']
        lines_common=[' radius = dims(ideogram,radius_outer)',' multiplier = 1e-6',
                ' thickness = 4p',' size = 20p',' label_offset = 5p',' label_separation = 5p','']
        tick1=["  spacing= 0.1u","  color = blue","  show_label = yes",
               "  label_size = 30p", "  label_offset = 5p","  format = %.1f",
               "  grid = yes","  grid_color = black","  grid_start = 0.55r",
               "  grid_end = 0.95r","  grid_thickness = 1p", '']
        tick1=self.insert_tag2(tick1,'tick')
        tick2=["  spacing = 1u","  color = red","  show_label = yes",
               "  label_size = 60p","  label_offset = 5p","  format = %d",
               "  suffix	=Mb","  grid_start = 0.5r","  grid_end = 0.975r",
               "  rid_color = red","  grid_thickness = 4p","  grid = yes",'']
        tick2=self.insert_tag2(tick2,'tick')
        self.lines=lines_common+tick1+tick2
        self.insert_tag('ticks')
        self.lines=lines_show+self.lines
        return self.lines
        
    def conf_plots(self):
        self.lines=[
               #GC content",
               "<plot>", "  type = line", "  file = plot_GCcontent.txt",
               "  r1 = 0.85r","  r0 = 0.65r","  thickness = 5",
               "  max = 0.60",  "  min = 0", "  extend_bin = no",
               "  color = red", "  orientation = in",
               "</plot>", "",
               #GC skew",
               "<plot>", " type = line", " file = plot_GCskew.txt",
               " r1 = 0.65r"," r0 = 0.51r"," thickness = 5",
               " max = 0.49999999999999173"," min = -0.47826086956521324",
               " extend_bin = no"," color = red"," orientation = in","",
               " <rules>","  <rule>","  condition = var(value) < 0.0",
               "  color = blue","  </rule>"," </rules>","</plot>",]
        self.insert_tag('plots')
        return self.lines

    def conf_plots_collinear(self):
        self.conf_links(0.9)
        self.lines=self.insert_tag2([], 'plots')+self.lines
        return self.lines
    
    def conf_links(self, radius=0.5):
        ###genomic interaction",
        #
        link1=["  file = link_forward.txt","", "  bezier_radius = 0r", "  crest = 0.25",
               "  radius = "+str(radius)+"r", "  color     = red","  thickness = 2",""]
        link1=self.insert_tag2(link1,'link')
        #
        link2=["  file = link_reverse.txt","", "  bezier_radius = 0r", "  crest = 0.25",
               "  radius = "+str(radius)+"r", "  color     = blue","  thickness = 2",""]
        link2=self.insert_tag2(link2,'link')
        # make sure that the id field matches the required number-number format",
        rules=["   condition = var(id) =~ /(\d+)-(\d+)/",
               "   thickness = eval( my @match = \"var(id)\" =~ /(\d+)-(\d+)/; remap($match[0],1,100,1,10) )",
               "   z = eval( my @match = \"var(id)\" =~ /(\d+)-(\d+)/; $match[0] )",
                "   color = eval( my @match = \"var(id)\" =~ /(\d+)-(\d+)/; sprintf(\"spectral-9-div-%d_a%d\", remap($match[1],1,100,1,9), remap($match[1],1,100,5,1 ) ) )",
            ]
        rules=self.insert_tag2(rules,'rule')
        rules=self.insert_tag2(rules,'rules')
        #
        self.lines=link1+link2+rules
        self.insert_tag('links')
        return self.lines

    def conf_highlights(self):
        h1=["  init_counter = highlight:1", " file = highlight_senseCDS.txt",
               " fill_color = red"," r1 = 0.99r"," r0 = 0.93r"]
        h2=["  file = highlight_senseCDS.txt",
            " fill_color = blue"," r1 = 0.93r"," r0 = 0.86r"]        
        self.lines=self.insert_tag2(h1,'highlight')+self.insert_tag2(h2,'highlight')
        self.insert_tag('highlights')
        return self.lines
    
    def export_conf(self,lines, file_name):
        outfile=self.args.dir_circos+file_name+'.conf'
        f1_obj = open(outfile, 'wt')
        f1_obj.write('\n'.join(lines))
        f1_obj.close()
        return outfile
    
    def append_conf(self,lines, file_name):
        outfile=self.args.dir_circos+file_name+'.conf'
        print(outfile)
        f1_obj = open(outfile, 'a')
        f1_obj.write('\n'.join(lines))
        f1_obj.close()
        return outfile

