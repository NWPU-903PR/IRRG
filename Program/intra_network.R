#' @title intra network
#' @description Building intracellular networks linked to receptors.
#'

#' @param goi gene of interest (typically a receptor)
#' @param data a data frame of n rows (genes) and m columns (cells) of read or UMI counts
#' @param cluster a numeric vector of length m
#' @param coi name of the cluster of interest
#' @param cell.prop a threshold, only the genes expressed in this proportion of the cells of the coi will be taken into account
#' @param c.names cluster names
#' @param species "homo sapiens" or "mus musculus"
#' @param connected a logical (if TRUE keeps only the genes connected to the goi)
#'
#' @return The function returns a list containing the gene regulatory network corresponding to each receptor 
#' and the expression information of the genes in the network


intra_network <- function(goi,data,cluster,coi,cell.prop,c.names=NULL,
                          species=c("homo sapiens","mus musculus"),
                          connected=FALSE,
                          mm2Hs,
                          PwC_ReactomeKEGG){

  if (is.null(c.names)==TRUE){
    c.names <- paste("cluster",seq_len(max(cluster)))
  }
  if (min(cluster)!=1){
    cluster <- cluster + 1 - min(cluster)
  }
  if (length(c.names)!=max(cluster) | sum(duplicated(c.names))>0 |
      grepl("/",paste(c.names,collapse =""))){
    stop("The length of c.names must be equal to the number of clusters and must
        contain no duplicates. The cluster names must not include special
        characters")
  }
  if (!is.element(coi,c.names)){
    stop(paste(coi,"must be included in c.names.","If c.names is not provided, it is set to cluster 1, cluster 2, ...,
        cluster N. WIth N the maximum number of clusters"))
  }
  #opar <- par()
  species <- match.arg(species)
  if (species=='mus musculus'){
    Mouse <- mm2Hs[,1]
    Human <- mm2Hs[,2]
    names(Human) <- Mouse
    names(Mouse) <- as.character(Human)
    m.names <- Human[rownames(data)]
    data <- subset(data,(!is.na(m.names)))
    m.names <- m.names[!is.na(m.names)]
    rownames(data) <- as.character(m.names)
    goi <- Human[goi]
  }
  
  col<-colnames(data)
  row<-rownames(data)
  data.tmp <- as.data.frame(data[,cluster==which(c.names==coi)])
  col<-col[cluster==which(c.names==coi)]
  row<-row[rowSums(data.tmp)>0]
  data.tmp <- as.data.frame(data.tmp[rowSums(data.tmp)>0,])
  
  colnames(data.tmp)<-col
  rownames(data.tmp)<-row
  
  
  # expressed genes
  good <- apply(data.tmp,1,function(x) sum(x>0)/ncol(data.tmp)>cell.prop)
  
  for (receptors in goi){
    if (!is.element(receptors,rownames(data.tmp))){
      a <- NULL 
      b<-NULL 
    }else {
      if (!receptors%in%unique(c(PwC_ReactomeKEGG$a.gn,PwC_ReactomeKEGG$b.gn))){
        a<-NULL
        b<-NULL
      }else{
        r_contain<-PwC_ReactomeKEGG[PwC_ReactomeKEGG$a.gn==receptors 
                                    |PwC_ReactomeKEGG$b.gn==receptors,]
        r_contain_gene<-unique(c(r_contain$a.gn,r_contain$b.gn))
        r_linked<-intersect(r_contain_gene,rownames(data.tmp))
        
        visible.genes <- unique(c(rownames(data.tmp)[good],receptors,r_linked))
        visible.n <- PwC_ReactomeKEGG[PwC_ReactomeKEGG$a.gn%in%visible.genes &
                                        PwC_ReactomeKEGG$b.gn%in%visible.genes,]
        if (!receptors%in%unique(c(visible.n$a.gn,visible.n$b.gn))){
          a<-NULL
          b<-NULL
        }else{
          red.visible.n <- simplify_interactions(visible.n,LRdb,autocrine = TRUE)
          net.n<-red.visible.n
          
          if (nrow(net.n)>0){
            add.net <- NULL
            nam=NULL
            siz=NULL
            
            
            
            net.tmp <- cbind(net.n[,seq_len(2)],type=net.n$type)
            net.f <- rbind(add.net,net.tmp)
            g.net <- graph_from_data_frame(net.f,directed=TRUE)
            g.net.tmp <- as.undirected(g.net)
            y <- shortest_paths(g.net.tmp,receptors,V(g.net.tmp))
            library(igraph)
            g.net.tmp <- graph_from_data_frame(net.f[,seq_len(2)],directed=FALSE)
            V(g.net.tmp)$status <- "pw.related"
            V(g.net.tmp)$status[
              which(unique(c(net.f$a.gn,net.f$b.gn)) %in% receptors)] =
              "gene.of.interest"
            V(g.net.tmp)$status[which(unique(c(net.f$a.gn,net.f$b.gn)) %in%
                                        net.f$a.gn[net.f$location=="extra"])] =
              "ligand"
            E(g.net.tmp)$int.type <- as.character(net.f$type)
            
            y <- unlist(lapply(y$vpath,function(x) length(x)))
            y <- y-1
            if (sum(y==-1)>0 & connected==FALSE){
              y[y==-1] <- 1
            }
            if (sum(y==-1)>0 & connected==TRUE){
              nam.tmp <- unique(c(net.f$a.gn,net.f$b.gn))[y==-1]
              net.f <- net.f[!net.f$a.gn%in%nam.tmp & !net.f$b.gn%in%nam.tmp,]
              y <- y[y!=-1]
            }
            data.tmp_list<-data.tmp
            if(species=='mus musculus'){
              net.f$a.gn = Mouse[net.f$a.gn]
              net.f$b.gn = Mouse[net.f$b.gn]
              rownames(data.tmp_list)<-as.character(Mouse[rownames(data.tmp_list)])
            }
            a <- net.f
            net_gene<-unique(c(net.f$a.gn,net.f$b.gn))
            b<-as.data.frame(data.tmp_list[net_gene,])
            colnames(b)<-col
            rownames(b)<-net_gene
          }
          else {
            a <- NULL
            b<-NULL
        }
      }
    }
  }
}
  return(list(net = a,data = b))
}