

## Download all series suppl data
 gse_relevant <- read.table("I:/get geo info/GSE_rel_list.txt",sep="\t",header=FALSE,row.names=NULL)
series.vector <- as.character(gse_relevant$V1)
num_series <- length(series.vector)
# 
# 
# setwd("I:/GEO_data/450k new")
# 
# library(doParallel)
# registerDoParallel(cores=8)
# foreach (i=1:num_series, .packages='GEOquery',.verbose=TRUE) %dopar% {
#   a <-getGEOSuppFiles(series.vector[i])  
# }

# set working directory where you want files to be stored
setwd("I:/GEO_data/450k new")

library(GEOquery)
## Alternatively, get data for one sample
#gse <- 'GSE32079' 
#a <-getGEOSuppFiles(gse)  

 
count<-1
 gse <-series.vector[count]; count = count+1

list.files(paste("I:/GEO_data/450k new/",gse,sep=''))

## load data into R
#change nrows to -1 to download all 450k
signals <- read.table(gzfile("I:/GEO_data/450k new/GSE29290/GSE29290_Matrix_Signal.txt.gz"),nrows=-1,header=TRUE,row.names=1,skip=0,sep='\t',dec = ",")

colnames(signals)


#remove suffixes from colnames
colnames(signals) <- gsub(".Signal_A","",colnames(signals))
colnames(signals) <- gsub(".Signal_B","",colnames(signals))
colnames(signals) <- gsub(".Detection","",colnames(signals))
colnames(signals) <- gsub(".Pval","",colnames(signals))

# locate relevant samples
samples.all <- colnames(signals)[seq(1,length(colnames(signals)),3)]



# get info about samples in series

series.info <- read.table("I:/get geo info/GSE_samples/joined/GSE29290.txt",sep='\t',nrows=-1,row.names=1,,header=TRUE,fill=TRUE,na.strings=c("NA","","0"),quote="\"",stringsAsFactors=FALSE)


relevant.samples.idx <- which(as.numeric(series.info$relevant)==1)

pheno <- series.info[relevant.samples.idx,]
num_samples <- length(series.info[,1])

## seperate tissues

 
# Important: This variable will change each time 
identifiers.all <- paste('Sample_',relevant.samples.idx,sep='')
 
 
 
unique.pheno <- unique(pheno[,c('tissue','cell_type','disease')])
unique.pheno.num <- length(unique.pheno[,1])
pheno.tables <- vector("list",unique.pheno.num)
group.names <- vector("list",unique.pheno.num)
identifiers.sub <-  vector("list",unique.pheno.num)
   
   
for (i in 1:unique.pheno.num){
  identifiers.sub[[i]] <- subset(identifiers.all,(pheno$tissue==unique.pheno[i,'tissue'] | is.na(pheno$tissue)) & 
                                   (pheno$cell_type==unique.pheno[i,'cell_type'] | is.na(pheno$cell_type)) & 
                                   pheno$disease==unique.pheno[i,'disease'])
  pheno.tables[[i]] <- subset(pheno,(tissue==unique.pheno[i,'tissue'] | is.na(tissue)) & 
                                (cell_type==unique.pheno[i,'cell_type'] | is.na(cell_type)) & 
                                disease==unique.pheno[i,'disease'])
  if (!is.na(unique.pheno$tissue[i])){
    group.names[[i]] <- paste(unique.pheno$tissue[i],'_',unique.pheno$disease[i],sep='')} else {
      group.names[[i]] <- paste(unique.pheno$cell_type[i],'_',unique.pheno$disease[i],sep='')
    }
}

setwd(gse)
 
## run analysis seperately for each tissue-disease group
library(RnBeads)
for (i in 1:unique.pheno.num){
  
relevant.samples.loc <- match(as.character(identifiers.sub[[i]]),samples.all)

# assign  unmethylated, methylated and pvalue matrices
U <- data.matrix(signals[,seq(1,(3*num_samples-2),3)])[,relevant.samples.loc]
M <- data.matrix(signals[,seq(2,(3*num_samples-1),3)])[,relevant.samples.loc]
p.values <- data.matrix(signals[,seq(3,(3*num_samples),3)])[,relevant.samples.loc]

#run rnbeads preprecossing

logger.start(fname=NA)
rnb.options(disk.dump.big.matrices=TRUE)
rnb.options(enforce.memory.management=TRUE)
rnb.set <- new('RnBeadRawSet',pheno.tables[[i]],U=U,M=M,p.values=p.values,useff=TRUE)

#parallel.setup(8)

#rnb.raw.set.greedy <- rnb.execute.greedycut(rnb.raw.set)
rnb.set <- rnb.execute.snp.removal(rnb.set)$dataset

rnb.set <- rnb.execute.normalization(rnb.set,method="bmiq",bgcorr.method="methylumi.lumi")

betas.table <- meth(rnb.set,row.names=TRUE)
pvalue.high <- which(dpval(rnb.set)>0.05,arr.ind=TRUE)
betas.table[pvalue.high[,'row'],pvalue.high[,'col']] <- NA
destroy(rnb.set)
write.table(betas.table,paste(gse,'_',group.names[[i]],'.txt',sep=''),sep='\t',col.names=NA)

}