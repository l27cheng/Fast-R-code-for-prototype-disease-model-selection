
DMcmbn<-function(M1, M2)  
{   
   temp<-rep(1:M2, M1)
   temp1<-sort(rep(1:M1, M2),  method='quick')
   pair0<- cbind(temp1, temp)
   pair1 <-pair0[temp>temp1,]
   pair2 <-pair0[temp<temp1,] 
   return(list(pair1=pair1, pair2=pair2))
} 
#################################################################
getDummy<-function(dat)
 { m<-ncol(dat)
    Dummy <-matrix(0, nrow(dat), 3*m)
  for(i in 1:m) 
  { Dummy[, 3*i-2] <- (dat[,i]==0)
     Dummy[, 3*i-1] <- (dat[,i]==1)
     Dummy[, 3*i  ] <- (dat[,i]==2)
  } 
  return(Dummy)
 } 
 
##################################################################
MAFcal<-function(dummy)
{ n<-nrow(dummy) 
tab<-colSums(dummy)
maf<- (tab[seq(2, length(tab), 3) ] + tab[seq(3, length(tab), 3) ]*2)/n/2
return(maf)
} 
 
###################################################################################################################################
library(matrixStats) 
PTY_filter<-function(yy, Dummy, protodmall, protodm2all, MAFofdat, mafcut1, mafcut2, mm)   
 { 
 yy2<-cbind(1-yy, yy)
 N1<-colSums(yy2) 
 K<-nrow(protodmall)
 K2<-K*2
 N<-nrow(Dummy)

 M<-ncol(Dummy)/3
 pair0<-DMcmbn(M, M)
 pair1 <-pair0$pair1
 np<-nrow(pair1)
  
 #Dummy <- getDummy(dat)  
 
  print(gc())
 # sparseDummy<-Matrix(Dummy, sparse=T)
 # Dcase<-crossprod(sparseDummy[1:N1[2], ])
 # Dcont<-crossprod(sparseDummy[(N1[2]+1):N, ])
  
 Dcase<- eigenMapMatMult2(Dummy[1:N1[2], ])
 Dcont<- eigenMapMatMult2(Dummy[(N1[2]+1):N, ])
 rm(Dummy)
  print(gc())
intv<-rep(0, M)
for(i in 1:M) intv[i]<-MAFmatchsub(MAFofdat, i, mafcut1)
posall<- 6*(intv[pair1[,1]] -1) + intv[pair1[,2]]
  
	  tab<-matrix(0,9, np*2)
     for(i in 1:(M-1))  
      { i1<-3*i 
	    ii<-(i1-2):i1
		jj<-(i1+1):(3*M)
       tab1<-Dcase[ii,  jj]
       tab2<-Dcont[ii,  jj]
	   pos<-M*i-i*(i+1)/2
	   pos1<-(pos - M + i+1):pos
       tab[, pos1]<-c(tab1)
       tab[, pos1+np]<-c(tab2)
       }
   rm(Dcase, Dcont)
   print(gc())
   
   tab12<- matrix(0, K2, np*2)
   for(i in  1:ncol(protodmall))
   {
    pos<- which(posall==i) 
	pos<-c(pos, pos+np)
    protodm<-protodmall[ ,i]
	protodm2<-protodm2all[ ,i] 
    protodm2<-protodm2[!is.na(protodm2)]
    # tab12all[1:K, pos]<-mm[protodm,]%*%tab[,  pos]
    # tab12all[(K+1):(K+length(protodm2)), pos]<-mm[protodm2,]%*%tab[c(1,4,7,2,5,8,3,6,9),  pos]
     tab12[1:K, pos]<-eigenMatMult( mm[protodm,],  tab[,  pos])
	 tab12[(K+1):(K+length(protodm2)), pos] <- eigenMatMult( matrix(mm[protodm2,], ,9), tab[c(1,4,7,2,5,8,3,6,9),  pos])
   }
  rm(tab)
  print(gc())
 # tab12<- cbind(c(tab12all[,1:np]), c(tab12all[, (np+1):(np*2)]))
 # tab12<- cbind(tab12, N1[1]-tab12)  
 
  # N21<- tab12[, 1]+tab12[, 2]
  # N22<- N1[1]
  # temp<-abs(tab12[, 1]* (N22-tab12[,2]) - tab12[,2]* (N22-tab12[,1]))- N/2
  # temp[temp<0]<-0
  # chi2stat<- N*temp^2/N1mP/N21/N22
 
  tab121<-tab12[, 1:np]
  tab12<-tab12[, (np+1):(np*2)]

  temp<-abs(tab121* (N1[1]-tab12) - tab12* (N1[2]-tab121))- N/2
  tab121<- tab121+tab12
  print(gc())
  
  temp[temp<0]<-0
  tab12<- N/N1[1]/N1[2]*temp^2/tab121/(N- tab121)
  rm(temp)
  tab12[is.na(tab12)]<- -1
  chi2max<- colMaxs(tab12)
  chi2pos<- max.col(t(tab12), ties.method="first")

 rslt<-cbind(chi2max, chi2pos, posall)
 return(rslt)
#  return(chi2max)
}
 
 
#########################
MAFmatchsub<-function(MAFofdat, l1, mafcut1)
{
 mf1<-MAFofdat[l1] 
 temp<-which.min(abs(mafcut1-mf1))
 if(mafcut1[temp]<=mf1)  intv1<-temp+1 else  intv1<-temp   
 return(intv1)
}
##########################
PTY_filter2<-function(yy, Dummy1, Dummy2, protodmall, protodm2all, MAFofdat1, MAFofdat2, mafcut1, mafcut2, mm)   
 { 
 yy2<-cbind(1-yy, yy)
 N1<-colSums(yy2) 
 K<-nrow(protodmall)
 K2<-K*2
 N<-nrow(Dummy1)

 M1<-ncol(Dummy1)/3 
 M2<-ncol(Dummy2)/3 
 np<-M1*M2 
 # Dummy1 <- getDummy(dat1) 
 # Dummy2 <- getDummy(dat2) 
  
 Dcase<- eigenMapMatMult1(Dummy1[1:N1[2], ], Dummy2[1:N1[2], ])
 Dcont<- eigenMapMatMult1(Dummy1[(N1[2]+1):N, ], Dummy2[(N1[2]+1):N, ])
 rm(Dummy1, Dummy2)
  print(gc())
intv1<-rep(0, M1)
intv2<-rep(0, M2)
for(i in 1:M1) intv1[i]<-MAFmatchsub(MAFofdat1, i, mafcut1)
for(i in 1:M2) intv2[i]<-MAFmatchsub(MAFofdat2, i, mafcut1)
posall<- t(outer(6*(intv1-1),  intv2, "+"))
  
	  tab<-matrix(0,9, np*2)
     for(i in 1:M1)  
      { i1<-3*i 
	    ii<-(i1-2):i1 
       tab1<-Dcase[ii, ]
       tab2<-Dcont[ii, ] 
	   pos<- ((i-1)*M2+1):(i*M2)
       tab[, pos]<-c(tab1)
       tab[, pos+np]<-c(tab2)
       }
   rm(Dcase, Dcont)
   print(gc())
   
   tab12<- matrix(0, K2, np*2)
   for(i in  1:ncol(protodmall))
   {
    pos<- which(posall==i) 
	pos<-c(pos, pos+np)
    protodm<-protodmall[ ,i]
	protodm2<-protodm2all[ ,i] 
    protodm2<-protodm2[!is.na(protodm2)]
     tab12[1:K, pos]<-eigenMatMult( mm[protodm,],  tab[,  pos])
	 tab12[(K+1):(K+length(protodm2)), pos] <- eigenMatMult( mm[protodm2,], tab[c(1,4,7,2,5,8,3,6,9),  pos])
   }
  rm(tab)
  print(gc()) 
  
  tab121<-tab12[, 1:np]
  tab12<-tab12[, (np+1):(np*2)]

  temp<-abs(tab121* (N1[1]-tab12) - tab12* (N1[2]-tab121))- N/2
  tab121<- tab121+tab12
  print(gc())
  
  temp[temp<0]<-0
  tab12<- N/N1[1]/N1[2]*temp^2/tab121/(N- tab121)
  rm(temp)
  tab12[is.na(tab12)]<- -1
  chi2max<- colMaxs(tab12)
 # chi2pos<- max.col(t(tab12), ties.method="first")

  #rslt<-cbind(chi2max, chi2pos, c(posall))
  return(chi2max)
}
 
###################################################################################################################################
   
MDR1_filter<-function(yy,  dat)   
 { 
  yy2<-cbind(1-yy, yy)
  N1<- colSums(yy2) 
 thresh<-N1[2]/N1[1]
 N<-nrow(dat)
 M<-ncol(dat) 
 pair0<-DMcmbn(M, M)
 pair1 <-pair0$pair1
 np<-nrow(pair1)   
 Dummy <- getDummy(dat)  
 Dcase<- eigenMapMatMult2(Dummy[1:N1[2], ])
 Dcont<- eigenMapMatMult2(Dummy[(N1[2]+1):N, ])
 DM<-(Dcase/Dcont>=thresh)*1
 
  tab12<-matrix(0,9, np*2)
  cndM<-matrix(0,9,np)
     for(i in 1:(M-1))  
      { i1<-3*i 
	    ii<-(i1-2):i1
		jj<-(i1+1):(3*M)
       tab1<-Dcase[ii,  jj]
       tab2<-Dcont[ii,  jj]
	   pos<-M*i-i*(i+1)/2
	   pos1<-(pos - M + i+1):pos
       tab12[, pos1]<-c(tab1)
       tab12[, pos1+np]<-c(tab2)
	   cndM[, pos1]<-DM[ii, jj]
       }  
   
  tab12[, 1:np]<- cndM*tab12[,  1:np]
  tab12[, (np+1):(2*np)]<- cndM*tab12[,  (np+1):(2*np)]
  tab12<-colSums(tab12, na.rm=T)
  
  tab121<-tab12[1:np]
  tab12<-tab12[ (np+1):(np*2)]

  temp<-abs(tab121* (N1[1]-tab12) - tab12* (N1[2]-tab121))- N/2
  tab121<- tab121+tab12
  #temp[temp<0]<-0
  tab12<- N/N1[1]/N1[2]*temp^2/tab121/(N- tab121)
	   
  return(tab12)
}
 
############################

MDR1_filter2<-function(yy,  dat1, dat2 )   
 { 
  yy2<-cbind(1-yy, yy)
  N1<- colSums(yy2) 
 thresh<-N1[2]/N1[1]
 N<-nrow(dat1)
 M1<-ncol(dat1) 
 M2<-ncol(dat2)   
 np<-M1*M2 
 Dummy1 <- getDummy(dat1) 
 Dummy2 <- getDummy(dat2) 
  rm(dat1, dat2)
 Dcase<- eigenMapMatMult1(Dummy1[1:N1[2], ], Dummy2[1:N1[2], ])
 Dcont<- eigenMapMatMult1(Dummy1[(N1[2]+1):N, ], Dummy2[(N1[2]+1):N, ])
 DM<-(Dcase/Dcont> thresh)*1
 rm(Dummy1, Dummy2)
 
  tab12<-matrix(0,9, np*2)
  cndM<-matrix(0,9,np)
     for(i in 1:M1)  
      { i1<-3*i 
	    ii<-(i1-2):i1 
       tab1<-Dcase[ii,  ]
       tab2<-Dcont[ii,  ]
	   pos<-((i-1)*M2+1):(i*M2)
       tab12[, pos]<-c(tab1)
       tab12[, pos+np]<-c(tab2)
	   cndM[, pos]<-DM[ii,  ]
       }    
     rm(Dcase, Dcont)
	 
  tab12[, 1:np]<- cndM*tab12[,  1:np]
  tab12[, (np+1):(2*np)]<- cndM*tab12[,  (np+1):(2*np)]
  tab12<-colSums(tab12, na.rm=T)
  
  tab121<-tab12[1:np]
  tab12<-tab12[ (np+1):(np*2)]

  temp<-abs(tab121* (N1[1]-tab12) - tab12* (N1[2]-tab121))- N/2
  tab121<- tab121+tab12
  temp[temp<0]<-0
  tab12<- N/N1[1]/N1[2]*temp^2/tab121/(N- tab121)

  return(tab12)
}
 
###################################################################################################################################

MDR9cutPTY_filter2<-function(yy, dat1, dat2, protodmall, protodm2all,   mafcut1, mafcut2, mm)
{
 yy2<-cbind(1-yy, yy)
 N1<-colSums(yy2) 
 N<-nrow(dat1)
 M1<-ncol(dat1) 
 M2<-ncol(dat2) 
 np<-M1*M2 
 
 Dummy1 <- getDummy(dat1) 
 Dummy2 <- getDummy(dat2) 
  rm(dat1, dat2)
 Dcase<- eigenMapMatMult1(Dummy1[1:N1[2], ], Dummy2[1:N1[2], ])
 Dcont<- eigenMapMatMult1(Dummy1[(N1[2]+1):N, ], Dummy2[(N1[2]+1):N, ]) 
 
 MAFofdat1<-MAFcal(Dummy1[yy==0,])  
 MAFofdat2<-MAFcal(Dummy2[yy==0,])  
 rm(Dummy1, Dummy2)
 K<-nrow(protodmall)
 K2<-K*2
intv1<-rep(0, M1)
intv2<-rep(0, M2)
for(i in 1:M1) intv1[i]<-MAFmatchsub(MAFofdat1, i, mafcut1)
for(i in 1:M2) intv2[i]<-MAFmatchsub(MAFofdat2, i, mafcut1)
posall<- t(outer(6*(intv1-1),  intv2, "+"))
  
 tab<-matrix(0,9, np*2)
     for(i in 1:M1)  
      { i1<-3*i 
	    ii<-(i1-2):i1 
       tab1<-Dcase[ii, ]
       tab2<-Dcont[ii, ] 
	   pos<- ((i-1)*M2+1):(i*M2)
       tab[, pos]<-c(tab1)
       tab[, pos+np]<-c(tab2)
       }  
  rm(Dcase, Dcont, tab1, tab2)	    
  
  ratio<-tab[ ,1:np]/tab[ , (np+1):(2*np)]
  cndM<-(ratio> sum(yy==1)/sum(yy==0))*1
  tab12<-tab0<-tab
  tab12<-cbind( cndM, cndM)*tab12
  rm(cndM)
  tab12<-colSums(tab12, na.rm=T) 
  tab121<-tab12[1:np]
  tab12<-tab12[ (np+1):(np*2)]
  temp<-abs(tab121* (N1[1]-tab12) - tab12* (N1[2]-tab121))- N/2
  tab121<- tab121+tab12
  rm(tab12)
  temp[temp<0]<-0
  chi2_mdr<- N/N1[1]/N1[2]*temp^2/tab121/(N- tab121)
rm(temp, tab121)

#9cut
  ratio[is.na(ratio)]<-1
  orderofRatio<-apply(ratio, 2,  order)
  rm( ratio )
      for(i in 1:np)  
       {  
	     tab[, i]<- cumsum(tab[orderofRatio[,i], i])
		 tab[, i+np]<- cumsum(tab[orderofRatio[,i], i+np])   
       } 
	   rm(orderofRatio)
  tab<-tab[-9, ]     
  tab121<-tab[ ,1:np]
  tab<-tab[ ,(np+1):(np*2)] 
  temp<-abs(tab121* (N1[1]-tab) - tab* (N1[2]-tab121))- N/2
  tab121<- tab121+tab
  temp[temp<0]<-0
  tab<- N/N1[1]/N1[2]*temp^2/tab121/(N- tab121)
  rm(temp, tab121)
  tab[is.na(tab)]<- -1
  chi2_rs<- colMaxs(tab)
rm(tab)

#pty     
  tab12<- matrix(0, K2, np*2)
   for(i in  1:ncol(protodmall))
   {
    pos<- which(posall==i) 
	pos<-c(pos, pos+np)
    protodm<-protodmall[ ,i]
	protodm2<-protodm2all[ ,i] 
    protodm2<-protodm2[!is.na(protodm2)]
     tab12[1:K, pos]<-eigenMatMult( mm[protodm,],  tab0[,  pos])
	 tab12[(K+1):(K+length(protodm2)), pos] <- eigenMatMult( mm[protodm2,], tab0[c(1,4,7,2,5,8,3,6,9),  pos])
   }
  rm(tab0)
  print(gc()) 
  tab121<-tab12[, 1:np]
  tab12<-tab12[, (np+1):(np*2)]
  temp<-abs(tab121* (N1[1]-tab12) - tab12* (N1[2]-tab121))- N/2
  tab121<- tab121+tab12
  print(gc())
  
  temp[temp<0]<-0
  tab12<- N/N1[1]/N1[2]*temp^2/tab121/(N- tab121)
  rm(temp)
  tab12[is.na(tab12)]<- -1
  chi2_pty<- colMaxs(tab12)
 
return(cbind(chi2_mdr, chi2_rs, chi2_pty))
} 

 
###################################################################################################################################
  
funPrtMDR1fwd<-function(yy, dat, protodmall, protodm2all, Cnt, Sm1, mafcut1, mafcut2, mm)
{
 yy2<-cbind(1-yy, yy)
 N1<-colSums(yy2) 
  thresh<-N1[2]/N1[1]
 K<-nrow(protodmall)
 K2<-K*2
 N<-nrow(dat)
 M<-ncol(dat) 
 pair0<-DMcmbn(M, M)
 pair1 <-pair0$pair1
 np<-nrow(pair1)   
   
 Dummy <- getDummy(dat) 
 rm(dat) 
 Dcase<- eigenMapMatMult2(Dummy[1:N1[2], ])
 Dcont<- eigenMapMatMult2(Dummy[(N1[2]+1):N, ])
 DM<-(Dcase/Dcont>=thresh)*1 
  MAFofdat<-MAFcal(Dummy[yy==0,])  
 rm(Dummy)
 
intv<-rep(0, M)
for(i in 1:M) intv[i]<-MAFmatchsub(MAFofdat, i, mafcut1)
posall<- 6*(intv[pair1[,1]] -1) + intv[pair1[,2]]
	   
   tab<-matrix(0,9, np*2)
     cndM<-matrix(0,9,np)
     for(i in 1:(M-1))  
      { i1<-3*i 
	    ii<-(i1-2):i1
		jj<-(i1+1):(3*M)
       tab1<-Dcase[ii,  jj]
       tab2<-Dcont[ii,  jj]
	   pos<-M*i-i*(i+1)/2
	   pos1<-(pos - M + i+1):pos
       tab[, pos1]<-c(tab1)
       tab[, pos1+np]<-c(tab2)
	   cndM[, pos1]<-DM[ii, jj]
       }
   rm(Dcase, Dcont)
	   
   tab12<- matrix(0, K2, np*2)
   for(i in  1:ncol(protodmall))
   {
    pos<- which(posall==i) 
	pos<-c(pos, pos+np)
    protodm<-protodmall[ ,i]
	protodm2<-protodm2all[ ,i] 
    protodm2<-protodm2[!is.na(protodm2)]
	
    sm1<-Sm1[[i]]
	cnt<-Cnt[, i]
    dm1 <- as.numeric(names(sort(sm1[which(is.element(as.numeric(colnames(sm1)), cnt)), dm]))[1:3])
	protodm<-c(protodm, dm1)

    tab12[1:K, pos]<-eigenMatMult( mm[protodm,],  tab[,  pos])
    tab12[(K+1):(K+length(protodm2)), pos] <- eigenMatMult( matrix(mm[protodm2,], ,9), tab[c(1,4,7,2,5,8,3,6,9),  pos])
   }

  tab[, 1:np]<- cndM*tab[,  1:np]
  tab[, (np+1):(2*np)]<- cndM*tab[,  (np+1):(2*np)]
  tab<-colSums(tab, na.rm=T)
  tab12<-rbind(tab, tab12)
    rm(tab)  
   
  # tab1<-tab[1:np]
  # tab<-tab[ (np+1):(np*2)]
  # temp<-abs(tab1* (N1[1]-tab) - tab* (N1[2]-tab1))- N/2
  # tab1<- tab1+tab
  # #temp[temp<0]<-0
  # tab<- N/N1[1]/N1[2]*temp^2/tab1/(N- tab1)
  
   tab121<-tab12[, 1:np]
  tab12<-tab12[, (np+1):(np*2)]

  temp<-abs(tab121* (N1[1]-tab12) - tab12* (N1[2]-tab121))- N/2
  tab121<- tab121+tab12
  
  temp[temp<0]<-0
  tab12<- N/N1[1]/N1[2]*temp^2/tab121/(N- tab121)
  rm(temp)
  tab12[is.na(tab12)]<- -1
  chi2max<- colMaxs(tab12)
  chi2pos<- max.col(t(tab12), ties.method="first")

 rslt<-cbind(chi2max, chi2pos, posall)
 return(rslt)
} 

###################################################################################################################################

MDR9cut_filter2<-function(yy, dat1, dat2)
{
 yy2<-cbind(1-yy, yy)
 N1<-colSums(yy2) 
  thresh<-N1[2]/N1[1]
 N<-nrow(dat1)
 M1<-ncol(dat1) 
 M2<-ncol(dat2) 
 np<-M1*M2 
 Dummy1 <- getDummy(dat1) 
 Dummy2 <- getDummy(dat2) 
  rm(dat1, dat2)
 Dcase<- eigenMapMatMult1(Dummy1[1:N1[2], ], Dummy2[1:N1[2], ])
 Dcont<- eigenMapMatMult1(Dummy1[(N1[2]+1):N, ], Dummy2[(N1[2]+1):N, ]) 
 
  rm(Dummy1, Dummy2)
  tab<-matrix(0,9, np*2)
 # cndM<-matrix(0,9,np)
     for(i in 1:M1)  
      { i1<-3*i 
	    ii<-(i1-2):i1 
       tab1<-Dcase[ii,  ]
       tab2<-Dcont[ii,  ]
	   pos<-((i-1)*M2+1):(i*M2) 
       tab[, pos]<-c(tab1)
       tab[, pos+np]<-c(tab2)
	  # cndM[, pos]<-ratio[ii, ]
       }  
  rm(Dcase, Dcont, tab1, tab2)	   
  ratio<-tab[ ,1:np]/tab[ , (np+1):(2*np)]
  cndM<-(ratio>thresh)*1
  tab12<-tab
  tab12<-cbind( cndM, cndM)*tab12
  rm(cndM)
  tab12<-colSums(tab12, na.rm=T)
  
  tab121<-tab12[1:np]
  tab12<-tab12[ (np+1):(np*2)]
  temp<-abs(tab121* (N1[1]-tab12) - tab12* (N1[2]-tab121))- N/2
  tab121<- tab121+tab12
  temp[temp<0]<-0
  tab12<- N/N1[1]/N1[2]*temp^2/tab121/(N- tab121)

 
  ratio[is.na(ratio)]<-1
  orderofRatio<-apply(ratio, 2,  order)
  rm( ratio )
 
      for(i in 1:np)  
       {  
	     tab[, i]<- cumsum(tab[orderofRatio[,i], i])
		 tab[, i+np]<- cumsum(tab[orderofRatio[,i], i+np])   
       } 
	   rm(orderofRatio)
	tab<-tab[-9, ]     
  tab121<-tab[ ,1:np]
  tab<-tab[ ,(np+1):(np*2)]

  temp<-abs(tab121* (N1[1]-tab) - tab* (N1[2]-tab121))- N/2
  tab121<- tab121+tab
  temp[temp<0]<-0
  tab<- N/N1[1]/N1[2]*temp^2/tab121/(N- tab121)
  rm(temp)
  tab[is.na(tab)]<- -1
 chi2max<- colMaxs(tab)

return(cbind(tab12, chi2max))
} 
  
 
###################################################################################################################################

fun9cut_filter2<-function(yy, dat1, dat2)
{
 yy2<-cbind(1-yy, yy)
 N1<-colSums(yy2) 
 N<-nrow(dat1)
 M1<-ncol(dat1) 
 M2<-ncol(dat2) 
 np<-M1*M2 
 Dummy1 <- getDummy(dat1) 
 Dummy2 <- getDummy(dat2) 
  rm(dat1, dat2)
 Dcase<- eigenMapMatMult1(Dummy1[1:N1[2], ], Dummy2[1:N1[2], ])
 Dcont<- eigenMapMatMult1(Dummy1[(N1[2]+1):N, ], Dummy2[(N1[2]+1):N, ]) 
 ratio<-Dcase/Dcont
  rm(Dummy1, Dummy2)
  tab12<-matrix(0,9, np*2)
  cndM<-matrix(0,9,np)
     for(i in 1:M1)  
      { i1<-3*i 
	    ii<-(i1-2):i1 
       tab1<-Dcase[ii,  ]
       tab2<-Dcont[ii,  ]
	   pos<-((i-1)*M2+1):(i*M2) 
       tab12[, pos]<-c(tab1)
       tab12[, pos+np]<-c(tab2)
	   cndM[, pos]<-ratio[ii, ]
       }  
  cndM[is.na(cndM)]<-1
  orderofRatio<-apply(cndM, 2,  order)
  rm(Dcase, Dcont, ratio, cndM, tab1, tab2)
 # tab12<-matrix(0, 9, np*2)
      for(i in 1:np)  
       {  
	     tab12[, i]<- cumsum(tab12[orderofRatio[,i], i])
		 tab12[, i+np]<- cumsum(tab12[orderofRatio[,i], i+np])   
       } 
	   rm(orderofRatio)
	tab12<-tab12[-9, ]     
  tab121<-tab12[ ,1:np]
  tab12<-tab12[ ,(np+1):(np*2)]

  temp<-abs(tab121* (N1[1]-tab12) - tab12* (N1[2]-tab121))- N/2
  tab121<- tab121+tab12
  temp[temp<0]<-0
  tab12<- N/N1[1]/N1[2]*temp^2/tab121/(N- tab121)
  tab12[is.na(tab12)]<- -1
 chi2max<- colMaxs(tab12)

return(chi2max)
}

############################################
fun9cut_filter<-function(yy, dat)
{
 yy2<-cbind(1-yy, yy)
 N1<-colSums(yy2) 
 N<-nrow(dat)
 M<-ncol(dat) 
 pair0<-DMcmbn(M, M)
 pair1 <-pair0$pair1
 np<-nrow(pair1)
 
  Dummy <- getDummy(dat) 
  #rm(dat)
  Dcase<- eigenMapMatMult2(Dummy[1:N1[2], ])
  Dcont<- eigenMapMatMult2(Dummy[(N1[2]+1):N, ])
  ratio<-Dcase/Dcont
  #rm(Dummy)
  
  tab12<-matrix(0,9, np*2)
  cndM<-matrix(0,9,np)
     for(i in 1:(M-1))  
      { i1<-3*i 
	    ii<-(i1-2):i1
		jj<-(i1+1):(3*M)
       tab1<-Dcase[ii,  jj]
       tab2<-Dcont[ii,  jj]
	   pos<-M*i-i*(i+1)/2
	   pos1<-(pos - M + i+1):pos
       tab12[, pos1]<-c(tab1)
       tab12[, pos1+np]<-c(tab2)
	   cndM[, pos1]<-ratio[ii, jj]
       }  
  cndM[is.na(cndM)]<-1
  orderofRatio<-apply(cndM, 2,  order)
  
       for(i in 1:np)  
       {  
	     tab12[, i]<- cumsum(tab12[orderofRatio[,i], i])
		 tab12[, i+np]<- cumsum(tab12[orderofRatio[,i], i+np])   
       } 
	tab12<-tab12[-9, ]     
  tab121<-tab12[ ,1:np]
  tab12<-tab12[ ,(np+1):(np*2)]

  temp<-abs(tab121* (N1[1]-tab12) - tab12* (N1[2]-tab121))- N/2
  tab121<- tab121+tab12
  temp[temp<0]<-0
  tab12<- N/N1[1]/N1[2]*temp^2/tab121/(N- tab121)
  tab12[is.na(tab12)]<- -1
 chi2max<- colMaxs(tab12)

return(chi2max)
}
 
#######################