#' @title simplify_interactions

#' @param t the network to be simplified
#' @param lr ligand receptor interactions
#' @param autocrine a logical

#' @return t

simplify_interactions <- function(t,lr=NULL,autocrine=FALSE){
  
  mergeText <- function(a,b){
    if (length(a)==0)
      b
    else
      if (length(b)==0)
        a
    else
      paste(union(strsplit(a,';')[[1]],strsplit(b,';')[[1]]),collapse=';')
  }
  
  
  t$detailed.type <- t$type
  
  # reduce interaction types
  t$type[t$type %in% c('interacts-with','in-complex-with')] <- 'complex'
  t$type[t$type %in% c('chemical-affects','consumption-controlled-by',
                       'controls-expression-of','controls-phosphorylation-of',
                       'controls-production-of','controls-state-change-of',
                       'controls-transport-of',
                       'controls-transport-of-chemical')] <- 'control'
  t$type[t$type %in% c('catalysis-precedes','reacts-with',
                       'used-to-produce')] <- 'reaction'
  
  # merge duplicated interactions with same reduced type
  key <- paste(t$a.gn,t$b.gn,t$type,sep='|')
  dk <- which(duplicated(key))
  for (i in dk){
    jj <- setdiff(which(t$a.gn==t[i,"a.gn"] & t$b.gn==t[i,"b.gn"] &t$type==t[i,"type"]),i)
    for (j in jj){
      t[j,"detailed.type"] <- mergeText(t[j,"detailed.type"],t[i,"detailed.type"])
    }
  }
  if (length(dk)>0)
    t <- t[-dk,]
  
  # promote interaction types in case one interaction is still given for
  # multiple "reduced" interaction types
  key <- paste(t$a.gn,t$b.gn,sep='|')
  dk <- which(duplicated(key))
  for (i in dk){
    jj <- setdiff(which(t$a.gn==t[i,"a.gn"] & t$b.gn==t[i,"b.gn"]),i)
    for (j in jj){
      t[j,"detailed.type"] <- mergeText(t[j,"detailed.type"],t[i,
                                                               "detailed.type"])
      if (t$type[j]=='control'){
        if (t$type[i]=='reaction' || t$type[i]=='complex')
          t$type[j] <- t$type[i]
      }
    }
  }
  if (length(dk)>0)
    t <- t[-dk,]
  
  # eliminate reverse interactions
  l.key <- paste(t$a.gn,t$b.gn,sep='|')
  r.key <- paste(t$b.gn,t$a.gn,sep='|')
  rk <- which(r.key%in%l.key)
  to.remove <- NULL
  for (i in rk){
    j <- setdiff(which(t$a.gn==t[i,"b.gn"] & t$b.gn==t[i,"a.gn"]),i)
    if (j < i){
      t[j,"detailed.type"] <- mergeText(t[j,"detailed.type"],t[i,"detailed.type"])
      if (t$type[j]!='control'){
        if (t$type[i]=='control'){
          t$type[j] <- 'control'
          t$a.gn[j] <- t$a.gn[i]
          t$b.gn[j] <- t$b.gn[i]
        }
        else
          if (t$type[i]=='complex')
            t$type[j] <- 'complex'
      }
      to.remove <- c(to.remove,i)
    }
  }
  if (length(to.remove)>0)
    t <- t[-to.remove,]
  
  if (!autocrine && !is.null(lr)){
    # remove LR pairs
    bad <- (t$a.gn%in%lr$ligand & t$b.gn%in%lr$receptor) |
      (t$a.gn%in%lr$receptor & t$b.gn%in%lr$ligand)
    t <- t[!bad,]
  }
  
  return(t)
}