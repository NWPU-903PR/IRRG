#' @title Receptor information
#' @description calculate the impact scores of receptors on downstream genes in each cell type

#' @detail This function will save the relevant information such as 
#' the impact score of the receptors in each cell type in the local "R_info" folder

#' @param goi gene of interest (typically a receptor)
#' @param data a data frame of n rows (genes) and m columns (cells) of read or UMI counts
#' @param cluster a numeric vector of length m
#' @param c.names  cluster names
#' @param cell.prop a threshold, only the genes expressed in this proportion of
#' the cells of the coi will be taken into account
#' @param write_netdat Locally preserve gene expression information in individual receptor gene regulatory networks
#' @param write_netdat Save the adjacency matrix information in each receptor gene regulatory network locally
#' @param LRdb Ligand-Receptor database
#' @param PWC pathway database
#' @param species "homo sapiens" or "mus musculus"


R_info<-function(data,
                 cluster,
                 c.names,
                 cell.prop,
                 write_netdat=T,
                 write_net=T,
                 LRdb,
                 mm2Hs,
                 PWC,
                 species
                 ){
  R.list<-unique(as.character(LRdb$receptor))
  if(species =='mus musculus'){
    L.list<-unique(as.character(LRdb$ligand))
    Human<-mm2Hs[,2]
    Mouse<-data.frame(mm2Hs[,1],row.names = Human)
    R.list_t<-as.character(Mouse[R.list,])
    R.combine=cbind(R.list,R.list_t)
    
    L.list_t<-as.character(Mouse[L.list,])
    L.combine = cbind(L.list,L.list_t)
    
    R.list_t[which(is.na(R.list_t))]<-str_to_title(
      R.combine[which(is.na(R.combine[,"R.list_t"])),"R.list"]
    )
    L.list_t[which(is.na(L.list_t))]<-str_to_title(
      L.combine[which(is.na(L.combine[,"L.list_t"])),"L.list"]
    )
    R.combine<-cbind(R.list_t,R.list)
    L.combine<-cbind(L.list_t,L.list)
    
    colnames(R.combine)<-colnames(mm2Hs)
    colnames(L.combine)<-colnames(mm2Hs)
    
    mm2Hs<-unique(rbind(mm2Hs,R.combine))
    mm2Hs<-unique(rbind(mm2Hs,L.combine))
    
    R.list<-R.list_t
    }
  
  for (coi in c.names) {
    print(coi)
    print(paste('Finding receptor network in',coi,'......'))
    pb <- txtProgressBar(min = 1, max = length(R.list), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  
    net_data<-foreach(goi=R.list[1:length(R.list)],
                      .packages = "igraph",
                      .options.snow = opts
                
                      )%dopar%{
                        source("./Program/intra_network.R")
                        source("./Program/simplify_interactions.R")
                        intra_network(goi = goi,data=data,
                                      cluster = cluster,
                                      coi = coi,
                                      mm2Hs=mm2Hs,
                                      PwC_ReactomeKEGG=PWC,
                                      cell.prop = cell.prop,
                                      c.names = c.names,
                                      species = species,
                                      connected =TRUE)}
    close(pb)
    nd<-net_data#
    #-------------------------------------------------------#
    re_goi<-vector()
    for (i in seq_len(length(nd))) {
      if (is.null(nd[[i]]$net) || 
          is.null(nd[[i]]$data) ||
          dim(nd[[i]][["data"]])[1]<10){ 
        re_goi<-c(re_goi,i)
      }
    }
    if(length(re_goi)>0){
      R<-R.list[-re_goi]
      nd[re_goi]<-NULL
    }else{
      R<-R.list
    }
    
    
    net_list<-list()
    dat_list<-list()
    for (i in seq_len(length(nd))) {
      net_list[[i]]<-nd[[i]]$net
      dat_list[[i]]<-as(as.matrix(nd[[i]]$dat),"dgCMatrix")
    }
    #-------------------------------------------------------#
    
    if(write_netdat==T){
      R.net_datapath<-"Rnet_data"
      if (dir.exists(R.net_datapath)==FALSE){
        dir.create(R.net_datapath)
      }
      save(dat_list,file = paste('./',R.net_datapath,'/',coi,'.rda',sep = '')) 
    }
    
    sepdat_ls<-list()
    sepdat_ls<-foreach(dat=dat_list[1:length(dat_list)] ,
                       receptor=R[1:length(R)],
                       .packages = 'Matrix' 
    )%dopar%{
      source('./Program/separate_data.R')
      separate_data(dat,receptor,S.threshold=0)
      }
    
    
    backgroup_ls<-list()
    expgroup_ls<-list()
    for (i in seq_len(length(sepdat_ls))) {
      backgroup_ls[[i]]<-sepdat_ls[[i]]$backgroup_ls
      expgroup_ls[[i]]<-sepdat_ls[[i]]$expgroup_ls  
      
    }
    #-------------------------------------------------------#
   
    adj_ls<-list()
    adj_ls<-foreach(net=net_list[1:length(net_list)]
    )%dopar%{
      source('./Program/net2adj.R')
      net2adj(net)
      }
  
    #保存网络结构
    if(write_net==T){
      R.netpath<-"R_net"
      if (dir.exists(R.netpath)==FALSE){
        dir.create(R.netpath)
      }
      save(adj_ls,file = paste('./',R.netpath,'/',coi,'.rda',sep = ''))
    }
    #-------------------------------------------------------------------------------------#
    print('Caculating receptor score......')
    goi_index_V<-c()
    for (i in seq_len(length(R))) {
      netgene_name<-rownames(as.matrix(adj_ls[[i]]))
      index<-which(netgene_name%in%R[i])
      goi_index_V<-c(goi_index_V,index)
    }
    
    Rcpp::sourceCpp('./Program/cal_Rscore.cpp')
    
    now<-cal_Rscore(exp_ls = expgroup_ls,
                    back_ls = backgroup_ls,
                    adj_list = adj_ls,
                    goi_vec = goi_index_V)
    colnames(now)<-c('R_score_noGY','out','in','total','part')
    rownames(now)<-R
    now<-as.data.frame(now)
    R_noGY<-now$R_score_noGY
    tmp<-which(R_noGY==1)
    if(length(tmp)>0){
      Rscore<-(R_noGY-min(R_noGY[-tmp]))/(max(R_noGY[-tmp])-min(R_noGY[-tmp]))
      Rscore[tmp]<-1
    }else{
      Rscore<-(R_noGY-min(R_noGY))/(max(R_noGY)-min(R_noGY))
    }
    
    now<-cbind(now,Rscore)
    colnames(now)<-c('R_score_noGY','out','in','total','part','R_score')
    R.infopath = 'R_info'
    if (dir.exists(R.infopath)==FALSE){
      dir.create(R.infopath)
    }
    write.csv(now,file = paste('./',R.infopath,'/',coi,'_Rinfo.csv',sep = ''))
    
  }
}