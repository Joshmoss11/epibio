

## Download all series suppl data
# gse_relevant <- read.table("I:/get geo info/GSE_rel_list.txt",sep="\t",header=FALSE,row.names=NULL)
# series.vector <- as.character(gse_relevant$V1)
# num_series <- length(series.vector)
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
gse <- 'GSE32079' 
a <-getGEOSuppFiles(gse)  


## load data into R
#change nrows to -1 to download all 450k
signals <- read.table(gzfile("I:/GEO_data/450k new/GSE32079/GSE32079_non-normalized.txt.gz"),nrows=-1,header=TRUE,row.names=1,skip=0,sep='\t',dec = ".")




#remove suffixes from colnames
colnames(signals) <- gsub(".Signal_A","",colnames(signals))
colnames(signals) <- gsub(".Signal_B","",colnames(signals))
colnames(signals) <- gsub(".Detection","",colnames(signals))
colnames(signals) <- gsub(".Pval","",colnames(signals))

# locate relevant samples
samples.all <- colnames(signals)[seq(1,(3*num_samples-2),3)]

#samples.all <- gsub(".Signal_A","",colnames(signals)[seq(1,(3*num_samples-2),3)])


# get info about samples in series

series.info <- read.table("I:/get geo info/GSE_samples/joined/GSE32079.txt",sep='\t',row.names=1,header=TRUE)

# relevant.samples.idx <- as.numeric(series.info$relevant)
# relevant.samples.idx[is.na(relevant.samples.idx)] <- 0
# relevant.samples.idx <- relevant.samples.idx == 1

relevant.samples.idx <- which(as.numeric(series.info$relevant)==1)

pheno <- series.info[relevant.samples.idx,]
num_samples <- length(series.info[,1])

# seperate tissues
unique.pheno <- unique(pheno[,c('tissue','cell_type','disease')])
unique.pheno.num <- length(unique.pheno[,1])
pheno.tables <- vector("list",unique.pheno.num)
group.names <- vector("list",unique.pheno.num)
for (i in 1:unique.pheno.num){
  pheno.tables[[i]] <- subset(pheno,(tissue==unique.pheno[i,'tissue'] | is.na(tissue)) & 
                                (cell_type==unique.pheno[i,'cell_type'] | is.na(cell_type)) & 
                                disease==unique.pheno[i,'disease'])
  if (!is.na(pheno$tissue[i])){
    group.names[[i]] <- paste(pheno$tissue[i],'_',pheno$disease[i],sep='')} else {
      group.names[[i]] <- paste(pheno$cell_type[i],'_',pheno$disease[i],sep='')
    }
}
  
#                        
# tissues <- unique(as.character(pheno$tissue))
# cell_types <- unique(as.character(pheno$cell_type))
# diseases <- unique(as.character(pheno$disease))
# 
# if (is.na(tissues)){
#   tissues.num <- 0} else{
#     tissues.num <- length(tissues)
#   }
# if (is.na(cell_types)){
#   cell_types.num <- 0} else{
#     cell_types.num <- length(cell_types)
#   }
# 
# phenos.all <- vector("list",tissues.num+cell_types.num)



for (i in 1:unique.pheno.num){
  
relevant.samples.loc <- match(as.character(pheno.tables[[i]]$description),samples.all)

#relevant.samples.loc <- c(1,2,3,4)
#pheno <- pheno[1:4,]

# assign  unmethylated, methylated and pvalue matrices
U <- data.matrix(signals[,seq(1,(3*num_samples-2),3)])[,relevant.samples.loc]
M <- data.matrix(signals[,seq(2,(3*num_samples-1),3)])[,relevant.samples.loc]
p.values <- data.matrix(signals[,seq(3,(3*num_samples),3)])[,relevant.samples.loc]

#run rnbeads preprecossing
rnb.raw.set <- new('RnBeadRawSet',pheno,U=U,M=M,p.values=p.values,useff=FALSE)

logger.start(fname=NA)
#parallel.setup(8)
rnb.raw.set.greedy <- rnb.execute.greedycut(rnb.raw.set)
rnb.raw.set.greedy.snprem <- rnb.execute.snp.removal(rnb.raw.set)$dataset
rnb.set.norm <- rnb.execute.normalization(rnb.raw.set.greedy.snprem,method="bmiq",bgcorr.method="methylumi.lumi")
rnb.set.sexrem <- rnb.execute.sex.removal(rnb.set.norm)$dataset
rnb.execute.export.csv(rnb.set.sexrem, NA)
#meth.beta <- meth(rnb.set.sexrem)
