#!/usr/bin/Rscript
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)


#
args=commandArgs(trailingOnly = T)
data.attr<-list('dir_result'=args[1])#pass parameters
#data.attr<-list('dir_result'='/mnt/rdedata/2014/PaleYellow_Leaf/')
data.attr$dir_ballgown<-paste(data.attr$dir_result, 'ballgown/', sep='')
data.attr$dir_stat<-paste(data.attr$dir_result, 'stat/', sep='')
if (dir.exists(data.attr$dir_stat)==F) {  dir.create(data.attr$dir_stat)}

#read sample names
sample_info<-read.csv(paste(data.attr$dir_result, 'sample_info.csv',sep=''), header=F)
data.attr$sample_names<-as.character(sample_info[,1])
sample_path<-paste(data.attr$dir_ballgown, data.attr$sample_names, sep='')
names(sample_path)<-data.attr$sample_names
data.attr$samples<-sample_path
print(data.attr$sample_names)


#
bg<-ballgown(samples = data.attr$samples, pData=data.attr$phenotype)
data.attr$bg<-bg
bg_filt<-subset(bg, "rowVars(texpr(bg))>5", genomesubset=T)
data.attr$fpkm=texpr(bg_filt, 'FPKM')
rownames(data.attr$fpkm)<-texpr(bg_filt, 'all')$gene_name


#distribution of gene abundance
fpkm_filt<-log2(data.attr$fpkm+1)
png(paste(data.attr$dir_stat, 'FPKM_abundance_distribution.png', sep=''))
boxplot(fpkm_filt, las=2, ylab='log2 FPKM',
        main='Gene abundance across samples, FPKM>5')
dev.off()
print (dim(fpkm_filt))

# read counts at exon level 
rc1<-eexpr(data.attr$bg, meas='all')
head(rc1, n=10)
#merge by exon id
rc2<-merge(indexes(bg)$e2t, rc1, by=c('e_id'))
rc3<-rc2[,sapply(colnames(rc2), FUN=function(x){grepl('^rcount|t_id',x)})]
head(rc3,n=10)
dim(rc3)
#aggregate by transcripts id
rc4<-aggregate(rc3, by=list(rc3$t_id), FUN=sum)
rc4$t_id=rc4$Group.1
rc4<-rc4[,-1]
head(rc4)
dim(rc4)
#


#fpkm and cov table
ep<-texpr(bg, meas='all')
dim(ep);head(ep[1:12])
length(ep$t_id);length(unique(ep$t_id))
#combine
ep<-merge(ep, rc4, by='t_id')
head(ep)
dim(ep)
#sorted_ep<-ep[order(ep$FPKM.capture9,decreasing=T),]
#head(sorted_ep)
#export fpkm table
fpkm.df<-ep[,!grepl('^cov|^rcount',colnames(ep))]
head(fpkm.df)
write.csv(fpkm.df, paste(data.attr$dir_stat, 'transcripts_expression_level_FPKM.csv', sep=''), row.names=F, quote=F)
##export rc table
rc.df<-ep[,!grepl('^cov|^FPKM',colnames(ep))]
head(rc.df)
dim(rc.df)
write.csv(rc.df, paste(data.attr$dir_stat, 'transcripts_expression_level_RC.csv', sep=''), row.names=F, quote=F)


#exon id vs. transcripts id
#a=indexes(bg)$e2t
#head(a)
#length(unique(a$e_id))
#length(unique(a$t_id))

