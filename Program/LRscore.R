#' @title LRscore
#' @description Calculate ligand-receptor interaction scores across all cell types

#' @detail This function will save all the ligand-receptor interaction information between cell types to the "LR_score_notest" folder, 
#' if "Permutation.test" is set to "TRUE", the interaction information of specific ligand-receptor pairs between each cell type 
#' will be saved in the "LR_score_test" folder.
#' 
#' @param data a data frame of n rows (genes) and m columns (cells) of read or UMI counts
#' @param c.names  cluster names
#' @param celltype List of cell types
#' @param species "homo sapiens" or "mus musculus". 
#' @param Permutation.test Logical parameter, whether to perform 
#' a permutation test to screen for specific ligand-receptor pairs
#' @param p_num Number of cell type permutations for permutation test
#' @param pVal_threshold  significance threshold


LRscore<-function(data,
                  c.names,
                  celltype,
                  species,
                  Permutation.test=TRUE,
                  p_num=1000,
                  pVal_threshold=0.05){
  Rinfo=1
  if (dir.exists('R_info')==FALSE){
    Rinfo=0
    print("There is no receptor score information yet, 
          you should execute the'R_info'function first; 
          The next calculation will only use the average expression information of the ligand and the receptor.")
  }
  LR_pair = LRdb[,c("ligand","receptor")]
  
  #配体受体库基因的改为小鼠对应的
  if(species =='mus musculus'){
    R.list<-unique(as.character(LRdb$receptor))
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
    
    Mouse<-as.data.frame(mm2Hs[,'Mouse gene name'])
    rownames(Mouse)<-mm2Hs[,'Gene name']
    LR_pair$ligand<-Mouse[as.character(LR_pair$ligand),]
    LR_pair$receptor<-Mouse[as.character(LR_pair$receptor),]
    
  }
  
  
  
  
  #============================================================================#
  #定义计算LRscore函数
  cal_LRscore <- function(coi.ligand,
                          coi.receptor,
                          R.info,
                          Rinfo_exist,
                          LR_pair,
                          cluster){
    data.ligand<-data[,cluster==which(c.names==coi.ligand)]
    data.receptor<-data[,cluster==which(c.names==coi.receptor)]
    # data.tmp <- data.tmp[rowSums(data.tmp)>0,]
    g.names<-rownames(data.ligand)
    cal.ligand <- function(m) {
      if (m %in% g.names) {
        m_expr<-mean(as.numeric(data.ligand[m,]))
      }else{
        m_expr<-0
      }
      return(m_expr)
    }
    cal.receptor <- function(m) {
      if (m %in% g.names) {
        m_expr<-mean(as.numeric(data.receptor[m,]))
      }else{
        m_expr<-0
      }
      return(m_expr)
    }
    loi.express<-lapply(as.character(LR_pair$ligand), cal.ligand)
    roi.express<-lapply(as.character(LR_pair$receptor), cal.receptor)
    if(Rinfo_exist){
      R.score<-list()
      for (R in LR_pair$receptor) {
        if(R%in%rownames(R.info)){
          R.score<-append(R.score,R.info[R,"R_score"])
        }else{
          R.score<-append(R.score,0)
        }
      }
      
      LR.info<-as.data.frame(cbind(as.numeric(loi.express),as.numeric(roi.express),as.numeric(R.score)))
      colnames(LR.info)<-c('loi.express','roi.express','R.score')
      
      LR.score <-as.data.frame(sqrt(LR.info$loi.express*LR.info$roi.express)*LR.info$R.score)
      colnames(LR.score)<-'LR.score'
      LR.express<-as.data.frame(sqrt(LR.info$loi.express*LR.info$roi.express))
      colnames(LR.express)<-'LR.express'
      
      res<-cbind(LR_pair,LR.info)
      res<-cbind(res,LR.express)
      res<-cbind(res,LR.score)
      
    }else{
      LR.info<-as.data.frame(cbind(as.numeric(loi.express),as.numeric(roi.express)))
      colnames(LR.info)<-c('loi.express','roi.express')
      
      LR.express<-as.data.frame(sqrt(LR.info$loi.express*LR.info$roi.express))
      colnames(LR.express)<-'LR.score'
      res<-cbind(LR_pair,LR.express)
    }
    return(res)
 
  }
  #-----------------------------------------------------------------------#
  res.path = 'LR_score_notest'
  if (dir.exists(res.path)==FALSE){
    dir.create(res.path)
  }
  for (coi.receptor in c.names) {
    if(Rinfo==1){
      R.info<-read.csv(file = paste('R_info/',coi.receptor,'_Rinfo.csv',sep = ''),row.names = 1)
      Rinfo_exist=TRUE
    }else{
      R.info=NA
      Rinfo_exist=FALSE
    }
      
    for (coi.ligand in c.names) {
      res<-cal_LRscore(coi.ligand = coi.ligand,
                       coi.receptor = coi.receptor,
                       R.info = R.info,
                       Rinfo_exist = Rinfo_exist,
                       LR_pair = LR_pair,
                       cluster = cluster)
      res<-res[which(res$LR.score!=0),]
      res<-res[order(-res$LR.score),]
      write.csv(res,file = paste('./',res.path,'/',coi.ligand,'_',coi.receptor,'.csv',sep = ''))
      cat(paste('The LRscore between ',coi.ligand,' and ',coi.receptor,' has been finished.\n',sep = ''))
      
    }
  }
  #---------------------------------------------------------------------------#
  if(Permutation.test==TRUE){
    res_sig.path = 'LR_score_test'
    if (dir.exists(res_sig.path)==FALSE){
      dir.create(res_sig.path)
    }
    
    for (coi.ligand in c.names) {
      for (coi.receptor in c.names) {
        res.true<-read.csv(file = paste('LR_score_notest/',coi.ligand,'_',coi.receptor,'.csv',sep = ''),row.names = 1)
        LR.express_test<-as.data.frame(matrix(nrow = nrow(res.true),ncol = p_num))
        loi.test<-(as.character(res.true$ligand))
        roi.test<-(as.character(res.true$receptor))
        permutation_test<-function(data,celltype,loi.test,roi.test){
          celltype.f<-sample(celltype)
          cluster.f<-vector()
          for (n in celltype.f) {
            cluster.f<-append(cluster.f,which(c.names == n))
          }
          data.test_l<-data[,cluster.f==which(c.names==coi.ligand)]
          data.test_r<-data[,cluster.f==which(c.names==coi.receptor)]
          
          data.test_l_means<-rowMeans(data.test_l)
          loi.means<-data.test_l_means[loi.test]
          data.test_r_means<-rowMeans(data.test_r)
          roi.means<-data.test_r_means[roi.test]
          lr_test<-as.numeric(sqrt(loi.means*roi.means))
          
          return(lr_test)
        }
        
        LR.express_test<-foreach(permute_idx=seq_len(p_num),
                                 .combine = cbind
        )%dopar%{
          permutation_test(data,celltype,loi.test,roi.test)}
        
        big_nums=as.data.frame(matrix(nrow = nrow(res.true),1))
        big_nums<-foreach(i=seq_len(nrow(res.true)),
                          .combine = rbind
        )%dopar%{
          length(which(LR.express_test[i,]>res.true[i,]$LR.express))}
        LR.express_pVal<-big_nums/p_num
        res.true<-cbind(res.true,LR.express_pVal)
        sig_idx<-which(big_nums<=p_num*pVal_threshold)
        res.true_sig<-res.true[sig_idx,]
        rownames(res.true_sig)<-seq_len(nrow(res.true_sig))

        res.true_sig<-res.true_sig[order(-res.true_sig$LR.score),]
        write.csv(res.true_sig,file = paste('./',res_sig.path,'/',coi.ligand,'_',coi.receptor,'.csv',sep = ''))
        cat(paste('The Permutation.test between ',coi.ligand,' and ',coi.receptor,' has been finished.\n',sep = ''))
      }
      
    }
  }
  
  info_total<-as.data.frame(matrix(nrow = length(c.names)^2,ncol = 4))
  i<-1
  for (ligand in c.names) {
    for (receptor in c.names) {
      if(Permutation.test==TRUE){
        res_utl_path<-paste('./',res_sig.path,'/',ligand,'_',receptor,'.csv',sep = '')
      }else{
        res_utl_path<-paste('./',res.path,'/',ligand,'_',receptor,'.csv',sep = '')
      }
      
      res_utl<-read.csv(file = res_utl_path,row.names = 1)
      info_total[i,]<-c(ligand,receptor,sum(res_utl$LR.score),nrow(res_utl))
      i<-i+1
    }
  }
  colnames(info_total)<-c('ligand','receptor','sum','total')
  
  summray_path<-'Celltype_communication_summary'
  if (dir.exists(summray_path)==FALSE){
    dir.create(summray_path)
  }
  write.csv(info_total,file = paste('./',summray_path,'/communication_info.csv',sep = ''))
  cat('The communication summary has been finished.\n')
  
}

  
  
  
  
  
  
  
  
  
  

