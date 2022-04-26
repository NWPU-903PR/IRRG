#' @title net to adj
#' @description convert the network to an adjacency matrix

#' @param net list containing multiple receptor network structures

#' @return list of adjacency matrices for individual receptor gene regulatory networks
net2adj<-function(net){
  gene<-unique(c(net$a.gn,net$b.gn))
  m<-which(net[,'type']%in%c('complex','reaction'))
  change<-net[m,][,c(1,2)]
  changed<-change[,c(2,1)]
  net<-net[,c(1,2)]
  colnames(changed)<-colnames(net)
  utl<-rbind(net,changed)
  adj<-as.matrix(as_adj(graph_from_data_frame(utl)))
  adj<-adj[gene,gene]
  adj_mat<-as(as.matrix(adj),"dgCMatrix")
  return(adj_mat)
}

