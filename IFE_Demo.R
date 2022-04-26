library(igraph)
library(doParallel)
library(foreach)
library(multtest)
library(stats)
library(stringr)
library(Matrix)
library(doSNOW)
library(Rcpp)

load("./Data/PwC_ReactomeKEGG.rda")
load("./Data/LRdb.rda")
load("./Data/mm2Hs.rda")

data<-read.csv(file = './Data/Mouse_IFE/IFE_express_normalized.csv',row.names = 1)
cell_ann<-read.csv(file = 'cell_ann.csv',row.names = 1)
celltype<-as.character(cell_ann$cell_type)
type_fre<-as.data.frame(table(celltype))
c.names <- as.character(type_fre$celltype)
cluster<-vector()
for (n in celltype) {
  cluster<-append(cluster,which(c.names == n))
}



parallel.core=20
cl<-makeSOCKcluster(parallel.core)
registerDoSNOW(cl)


source('./Program/R_info.R')
R_info(data = data,
      cluster = cluster,
      c.names = c.names,
      S.threshold =0,
      cell.prop = 0.2,
      LRdb=LRdb,
      mm2Hs = mm2Hs,
      PWC = PwC_ReactomeKEGG,
      species = "mus musculus")


source('./Program/LRscore.R')
LRscore(data = data,
        c.names = c.names,
        species = "mus musculus",
        celltype = celltype,
        Permutation.test = TRUE)


stopImplicitCluster()
stopCluster(cl)  

