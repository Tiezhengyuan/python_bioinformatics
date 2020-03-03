
#
####
get.data.attr<-function(){
  
  #
  data.attr<-list()
  data.attr$rc_file='/mnt/rdedata4/miRNA/results/stat/matured_miRNA_RC.txt'
  data.attr$dir_results=paste(dirname(data.attr$rc_file), '/',sep='')
  data.attr$pheno_file='/mnt/rdedata4/miRNA/results/phenotype.csv'
  #
  data.attr$rc_df=read.table(data.attr$rc_file, header=T, row.names=1)
  data.attr$rc_df=round(data.attr$rc_df)
  colnames(data.attr$rc_df)<-gsub('\\.', '-', colnames(data.attr$rc_df))
  dim(data.attr$rc_df)
  head(data.attr$rc_df)
  #
  data.attr$pheno_df=read.csv(data.attr$pheno_file, header=T, row.names=1)
  head(data.attr$pheno_df)
  head(data.attr$pheno_df)
  
  return(data.attr)
}

####
get.DEpairs<-function(){
  #extract pairs
  DEpairs<-list()
  
  #1st pairs
  pair1<-c('Leaf-F6-G61-30-0hr-vs-Leaf-F6-G61-88-0hr','Leaf-F6-G61-30-0hr-vs-Leaf-F6-G61-30-72hr',
           'Leaf-F6-G61-30-0hr-vs-Leaf-F6-G61-76-0hr','Leaf-F6-G61-30-72hr-vs-Leaf-F6-G61-88-72hr',
           'Leaf-F6-G61-30-72hr-vs-Leaf-F6-G61-76-72hr','Leaf-F6-G61-30-72hr-vs-Leaf-F6-G61-93-72hr',
           'Leaf-F6-G61-76-0hr-vs-Leaf-F6-G61-88-0hr','Leaf-F6-G61-76-0hr-vs-Leaf-F6-G61-76-72hr',
           'Leaf-F6-G61-76-72hr-vs-Leaf-F6-G61-88-72hr','Leaf-F6-G61-76-72hr-vs-Leaf-F6-G61-93-72hr',
           'Leaf-F6-G61-88-0hr-vs-Leaf-F6-G61-88-72hr','Leaf-F6-G61-88-72hr-vs-Leaf-F6-G61-93-72hr',
           'Root-F6-G61-30-0hr-vs-Root-F6-G61-93-0hr','Root-F6-G61-30-0hr-vs-Root-F6-G61-88-0hr',
           'Root-F6-G61-30-0hr-vs-Root-F6-G61-30-72hr','Root-F6-G61-30-0hr-vs-Root-F6-G61-76-0hr',
           'Root-F6-G61-30-72hr-vs-Root-F6-G61-88-72hr','Root-F6-G61-30-72hr-vs-Root-F6-G61-76-72hr',
           'Root-F6-G61-30-72hr-vs-Root-F6-G61-93-72hr','Root-F6-G61-76-0hr-vs-Root-F6-G61-88-0hr',
           'Root-F6-G61-76-0hr-vs-Root-F6-G61-93-0hr','Root-F6-G61-76-0hr-vs-Root-F6-G61-76-72hr',
           'Root-F6-G61-76-72hr-vs-Root-F6-G61-88-72hr','Root-F6-G61-76-72hr-vs-Root-F6-G61-93-72hr',
           'Root-F6-G61-88-0hr-vs-Root-F6-G61-88-72hr','Root-F6-G61-88-0hr-vs-Root-F6-G61-93-0hr',
           'Root-F6-G61-88-72hr-vs-Root-F6-G61-93-72hr','Root-F6-G61-93-0hr-vs-Root-F6-G61-93-72hr',
           'Leaf-Chacoa-0hr-vs-Leaf-Chacoa-72hr','Leaf-BigCuban-0hr-vs-Leaf-BigCuban-72hr',
           'Root-Chacoa-0hr-vs-Root-Chacoa-72hr','Root-BigCuban-0hr-vs-Root-BigCuban-72hr',
           'Leaf-Chacoa-0hr-vs-Leaf-BigCuban-0hr','Leaf-Chacoa-72hr-vs-Leaf-BigCuban-72hr',
           'Root-Chacoa-0hr-vs-Root-BigCuban-0hr','Root-Chacoa-72hr-vs-Root-BigCuban-72hr')
  
  data.attr$pheno_df$pair1=paste(data.attr$pheno_df$TissueType,
              data.attr$pheno_df$Variety, data.attr$pheno_df$Condition,sep='-')
  for( p in pair1){
    DEpairs[[p]]<-list()
    pp=strsplit(p, '-vs-')[[1]]
    print(pp)
    A=as.character(pp[1])
    DEpairs[[p]][[A]]<-rownames(data.attr$pheno_df[data.attr$pheno_df$pair1==A,])
    B=as.character(pp[2])
    DEpairs[[p]][[B]]<-rownames(data.attr$pheno_df[data.attr$pheno_df$pair1==B,])
    print(DEpairs[[p]])
  }
  
  #2rd pairs
  pair2=c('Leaf-HighAlk-0hr-vs-Leaf-LowAlk-0hr','Root-HighAlk-72hr-vs-Root-LowAlk-72hr',
          'Root-HighAlk-0hr-vs-Root-HighAlk-72hr','Root-LowAlk-0hr-vs-Root-LowAlk-72hr',
          'Leaf-HighAlk-0hr-vs-Leaf-HighAlk-72hr','Leaf-LowAlk-0hr-vs-Leaf-LowAlk-72hr',
          'Leaf-HighAlk-72hr-vs-Leaf-LowAlk-72hr','Root-HighAlk-72hr-vs-Root-LowAlk-0hr')
  data.attr$pheno_df$pair2=paste(data.attr$pheno_df$TissueType, data.attr$pheno_df$Phenotype,
                                 data.attr$pheno_df$Condition,sep='-')
  head(data.attr$pheno_df,n=3)
  for( p in pair2){
    DEpairs[[p]]<-list()
    pp<-strsplit(p, '-vs-')[[1]]
    print(pp)
    A=as.character(pp[1])
    DEpairs[[p]][[A]]<-rownames(data.attr$pheno_df[data.attr$pheno_df$pair2==A,])
    B=as.character(pp[2])
    DEpairs[[p]][[B]]<-rownames(data.attr$pheno_df[data.attr$pheno_df$pair2==B,])
    print(DEpairs[[p]])
  }
  
  
  return(DEpairs)
}

####
DEG<-function(data.attr,AB){
  #colnames(data.attr$rc_df)
  sub_names<-unlist(AB)
  countTable<-data.attr$rc_df[, sub_names]
  countTable<-countTable[apply(countTable,1, max)>5,]
  #head(countTable)
  condition<-sapply(names(AB), y=AB, FUN=function(x,y=y){
    #print(y[[x]])
    rep(x, each=length(y[[x]]))})
  condition<-factor(unlist(condition))
  Aname=names(AB)[1]
  Bname=names(AB)[2]
  #
  cds=newCountDataSet(countTable,condition)
  cds=estimateSizeFactors(cds)
  #sizeFactors(cds)
  if( length(AB[[1]])==1 | length(AB[[2]])==1 ){
    cds=estimateDispersions(cds,method='blind', sharingMode='fit-only')
  }else{
    cds=estimateDispersions(cds)
  }
  res=nbinomTest(cds, Aname,Bname)

  #head(res)
  
  #export df
  #head(countTable)
  nrc_df<-counts(cds,normalized=T)
  colnames(nrc_df)<-paste('normalized:',colnames(nrc_df),sep='')
  #head(nrc_df)
  comb_df<-cbind(countTable, nrc_df,res)
  #head(comb_df)
  write.csv(comb_df,paste(data.attr$prefix,'_DEG.csv',sep=''),
            row.names = T,quote=F)
  #
  #table(res$padj)
  resSig<-res[res$padj<0.1,]
  #resSig
  
  #plots
  png(paste(data.attr$prefix,'_DEG.png',sep=''),
      width=300,height=110,res=300, units='mm')
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
  #
  return(res)
}
#########################################################
library(DESeq)


#read data
data.attr<-get.data.attr()
names(data.attr)

#sample pairs
DEpairs<-get.DEpairs()

#
for (name in names(DEpairs)){
  print(name)
  #print(DEpairs[[name]])
  #print("")
  data.attr$prefix=paste(data.attr$dir_results,name,sep='')
  AB=DEpairs[[name]]
  DEG(data.attr,AB)
}


