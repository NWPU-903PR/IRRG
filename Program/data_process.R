dat  = read.csv(file ='./Data/IFE_express.csv',row.names = 1)
calc_cpm <-
  function (expr_mat, Spikein=Spikein_idx) 
  {
    norm_factor <- colSums(expr_mat[-Spikein, ])
    return(t(t(expr_mat)/norm_factor) * 10^6)
  }
keep_feaure<-rowSums(dat>0)>5
dat_flitergene<-dat[keep_feaure,]
cpm_data<-calc_cpm(dat_flitergene)
norm_dat<-log2(cpm_data+1)
colsum<-colSums(dat_flitergene)
write.csv(norm_dat,file = "./IFE_norm_data.csv")