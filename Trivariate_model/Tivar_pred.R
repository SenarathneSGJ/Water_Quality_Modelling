pred=function(TR.para,Fac.loc,dist.mat,z_new){
  
  Ey1.p = TR.para[1]+TR.para[2]*Fac.loc[,1]+TR.para[3]*Fac.loc[,2]+TR.para[4]*Fac.loc[,3]+TR.para[25]*z_new[,1]+TR.para[26]*z_new[,2]+TR.para[27]*z_new[,3]
  Ey2.p = TR.para[5]+TR.para[6]*Fac.loc[,1]+TR.para[7]*Fac.loc[,2]+TR.para[8]*Fac.loc[,3]+TR.para[28]*z_new[,1]+TR.para[29]*z_new[,2]+TR.para[30]*z_new[,3]
  Ey3.p = TR.para[9]+TR.para[10]*Fac.loc[,1]+TR.para[11]*Fac.loc[,2]+TR.para[12]*Fac.loc[,3]+TR.para[31]*z_new[,1]+TR.para[32]*z_new[,2]+TR.para[33]*z_new[,3]
  
  
  r21 <- 1/(1+exp(-TR.para[14]))
  r11 <- exp(TR.para[13])*r21
  
  Gdist1 <- exp(-1*(dist.mat/r21))
  Gd1 <- matrix(ncol=ncol(Gdist1),nrow=nrow(Gdist1))
  Gd1[Gdist1!=1] <- r11*Gdist1[Gdist1!=1]
  Gd1[Gdist1==1] <- r11+r01
  
  r22 <- 1/(1+exp(-TR.para[16]))
  r12 <- exp(TR.para[15])*r22
  
  Gdist2 <- exp(-1*(dist.mat/r22))
  Gd2 <- matrix(ncol=ncol(Gdist2),nrow=nrow(Gdist2))
  Gd2[Gdist2!=1] <- r12*Gdist2[Gdist2!=1]
  Gd2[Gdist2==1] <- r12+r02
  
  r23 <- 1/(1+exp(-TR.para[18]))
  r13 <- exp(TR.para[17])*r23
  
  Gdist3 <- exp(-1*(dist.mat/r23))
  Gd3 <- matrix(ncol=ncol(Gdist3),nrow=nrow(Gdist3))
  Gd3[Gdist3!=1] <- r13*Gdist3[Gdist3!=1]
  Gd3[Gdist3==1] <- r13+r03
  
  var.Y1= exp(TR.para[19])
  var.Y2= exp(TR.para[20])
  var.Y3= exp(TR.para[21])
  
  cov.Y1 = diag(rep(var.Y1,ncol(Gd1)))+Gd1
  cov.Y2 = diag(rep(var.Y2,ncol(Gd2)))+Gd2
  cov.Y3 = diag(rep(var.Y3,ncol(Gd3)))+Gd3
  
  sdY1 = sqrt(diag(cov.Y1))
  sdY2 = sqrt(diag(cov.Y2))
  sdY3 = sqrt(diag(cov.Y3))
  
  RR12 = 1/(1+exp(-TR.para[22]))
  C.Y1Y2 = RR12*sdY1*sdY2
  
  RR13 = 1/(1+exp(-TR.para[23]))
  C.Y1Y3 = RR13*sdY1*sdY3
  
  RR23 = 1/(1+exp(-TR.para[24]))
  C.Y2Y3 = RR23*sdY2*sdY3
  
  Cmat.Y1Y2=diag(C.Y1Y2)
  Cmat.Y1Y3=diag(C.Y1Y3)
  Cmat.Y2Y3=diag(C.Y2Y3)
  
  Cov.mat1 = cbind(cov.Y1,Cmat.Y1Y2,Cmat.Y1Y3)
  Cov.mat2 = cbind(t(Cmat.Y1Y2),cov.Y2,Cmat.Y2Y3)
  Cov.mat3 = cbind(t(Cmat.Y1Y3),t(Cmat.Y2Y3),cov.Y3)
  Cov.mat = rbind(Cov.mat1,Cov.mat2,Cov.mat3)
  
  N.obs= 3*length(Ey1.p)/nrow(Cov.mat)
  l.Y1=nrow(Cov.mat)/3
  
  Pred.Y =matrix(ncol=3,nrow=length(Ey1.p))
  for(i in 1:N.obs){
    StrtY = l.Y1*(i-1)+1
    EndY = l.Y1*i
    
    Ey =  c(Ey1.p[StrtY:EndY],Ey2.p[StrtY:EndY],Ey3.p[StrtY:EndY])
    Yi = rmvnorm(1,Ey,Cov.mat)
    Y.mat=matrix(exp(Yi),ncol=3)
    Pred.Y[StrtY:EndY,]=Y.mat
  }
  
  return(Pred.Y)
}