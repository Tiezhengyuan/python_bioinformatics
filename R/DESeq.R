
#
####
get.data.attr<-function(args){
  data.attr<-list()
  # test if there is at least one argument: if not, return an error
  if (length(args)<2) {
    stop("At least two arguments must be supplied (input file).n", call.=FALSE)
  } else if (length(args)==2) {
    # default output dir
    args[3] = dirname(args[1])
  }
  
  #

  data.attr$rc_file=args[1]
  data.attr$de_file=args[2]
  data.attr$dir_results=ifelse(grepl('/$', args[3]), args[3], paste(args[3],'/',sep=''))
  
  

  #retrieve RC table
  data.attr$rc_df=read.table(data.attr$rc_file, header=T, row.names=1, check.names=F)
  data.attr$rc_df=round(data.attr$rc_df)
  #colnames(data.attr$rc_df)<-gsub('\\.', '-', colnames(data.attr$rc_df))
  dim(data.attr$rc_df)
  head(data.attr$rc_df)
  print('Sample names:')
  print(colnames(data.attr$rc_df))
  #
  
  return(data.attr)
}

####
get.DEpair<-function(data.attr){
  #extract pairs
  DEpair<-list()
  sample_names=colnames(data.attr$rc_df)
  
  #p<-read.csv('/mnt/rdedata12/yuan/AltriaSeq/examples/statistics/rawdata/nanopore_pairs.csv',header=T)
  p<-read.csv(data.attr$de_file,header=T)
  names=as.character(unlist(p))
  names=names[!is.na(names)]
  #names
  
  diff<-setdiff(names, sample_names)
  if(length(diff)>0){
    print(paste("Some sample names are not detected:", diff))
    stop(paste("keep consistency between sample names of", 
        data.attr$rc_file, "and those of", data.attr$de_file, sep=' '), call.=F)
  }else{
    A_name=colnames(p)[1]
    B_name=colnames(p)[2]
    DEpair[[A_name]]=as.character(p[!is.na(p[,1]),1])
    DEpair[[B_name]]=as.character(p[!is.na(p[,2]),2])
    #DEpair
  }
  
  return(DEpair)
}

####
DEG<-function(data.attr,DEpair){
  #import library
  library(DESeq)
  
  #colnames(data.attr$rc_df)
  sub_names=unlist(DEpair)
  countTable<-data.attr$rc_df[,sub_names]
  countTable<-countTable[apply(countTable,1, max)>5,]
  #head(countTable)
  condition<-sapply(names(DEpair), y=DEpair, FUN=function(x,y=y){
    #print(y[[x]])
    rep(x, each=length(y[[x]]))})
  condition<-factor(unlist(condition))
  Aname=names(AB)[1]
  Bname=names(AB)[2]
  #
  cds=newCountDataSet(countTable,condition)
  cds=estimateSizeFactors(cds)
  #sizeFactors(cds)
  if( length(DEpair[[1]])==1 | length(DEpair[[2]])==1 ){
    cds=estimateDispersions(cds,method='blind', sharingMode='fit-only')
  }else{
    cds=estimateDispersions(cds)
  }
  res=nbinomTest(cds, Aname,Bname)
  #head(res)
  
  #export df
  file_prefix=paste(data.attr$dir_results,Aname,'_', Bname, sep='')
  #head(countTable)
  nrc_df<-counts(cds,normalized=T)
  colnames(nrc_df)<-paste('normalized:',colnames(nrc_df),sep='')
  #head(nrc_df)
  comb_df<-cbind(countTable, nrc_df,res)
  #head(comb_df)
  csv_file=paste(file_prefix,'_DESeq.csv',sep='')
  print(csv_file)
  write.csv(comb_df,csv_file, row.names = T,quote=F)
  #
  #table(res$padj)
  resSig<-res[res$padj<0.1,]
  #resSig

  #plots
  pic_file<-paste(file_prefix,'_DESeq.png',sep='')
  print(pic_file)
  png(pic_file, width=300,height=110,res=300, units='mm')
    par(mfrow=c(1,3))
    #MA plot
    plotMA(res, main='M-A plot')
    #
    hist(res$pval, breaks=100, main='Historgram of p-values')
    #
    plot(x=0:6,y=0:6,col='white', xlab='',ylab='')
      text(3,4, paste('group1:', Aname,  ' num=', length(AB[[1]]),sep=''))
      text(3,3, paste('group2:', Bname,' num=', length(AB[[1]]),sep=''))
      text(3,2, paste('adjusted-pvalue threshold:', 'pvalue=0.1'))
      text(3,1, paste('significance DEGs:', nrow(resSig)))
  dev.off()

  #return(res)
}



##########################################################
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


#########################################################
#read data
data.attr<-get.data.attr(args)
#names(data.attr)

#sample pairs
AB<-get.DEpair(data.attr)
#AB

#differntial analysis
DEG(data.attr, AB)



