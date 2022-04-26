#' @title separate_data
#' @description Segmentation of a cell type gene expression matrix based on receptor expression

#' @param Dat Gene expression matrix for a cell type
#' @param receptor single receptor
#' @param S.threshold segmentation threshold, i.e. receptor expression
#' 
#' @return segmented gene expression matrix

separate_data<-function(dat,receptor,S.threshold){
  dat<-as.matrix(dat)
  if(S.threshold==0){
    back_idx<-which(dat[receptor,]==0)
    exp_idx<-which(dat[receptor,]!=0) 
  }else{
    back_idx<-which(dat[receptor,]<=S.threshold)
    exp_idx<-which(dat[receptor,]>S.threshold) 
  }
  
  backgroup<-as.data.frame(dat[,back_idx],row.names = rownames(dat))
  colnames(backgroup)<-colnames(dat)[back_idx]
  expgroup<-as.data.frame(dat[,exp_idx],row.names = rownames(dat))
  colnames(expgroup)<-colnames(dat)[exp_idx]
  backgroup<-as(as.matrix(backgroup),"dgCMatrix")
  expgroup<-as(as.matrix(expgroup),"dgCMatrix")
  return(list(backgroup_ls=backgroup,expgroup_ls=expgroup))
}
  


