pred=function(TR.para,Fac.loc,dist.mat,z_new){
  
  Ey1.p = TR.para[1]+TR.para[2]*Fac.loc[,1]+TR.para[3]*Fac.loc[,2]+TR.para[4]*Fac.loc[,3]+TR.para[16]*z_new[,1]+TR.para[17]*z_new[,2]+TR.para[18]*z_new[,3]
  Ey2.p = TR.para[5]+TR.para[6]*Fac.loc[,1]+TR.para[7]*Fac.loc[,2]+TR.para[8]*Fac.loc[,3]+TR.para[19]*z_new[,1]+TR.para[20]*z_new[,2]+TR.para[21]*z_new[,3]
  
  r21 <- 1/(1+exp(-TR.para[10]))
  r11 <- exp(TR.para[9])*r21
  
  Gdist1 <- exp(-1*(dist.mat/r21))
  Gd1 <- matrix(ncol=ncol(Gdist1),nrow=nrow(Gdist1))
  Gd1[Gdist1!=1] <- r11*Gdist1[Gdist1!=1]
  Gd1[Gdist1==1] <- r11+r01
  
  r22 <- 1/(1+exp(-TR.para[12]))
  r12 <- exp(TR.para[11])*r22
  
  
  Gdist2 <- exp(-1*(dist.mat/r22))
  Gd2 <- matrix(ncol=ncol(Gdist2),nrow=nrow(Gdist2))
  Gd2[Gdist2!=1] <- r12*Gdist2[Gdist2!=1]
  Gd2[Gdist2==1] <- r12+r02
  
  var.Y1= exp(TR.para[13])
  var.Y2= exp(TR.para[14])
  
  cov.Y1 = diag(rep(var.Y1,ncol(Gd1)))+Gd1
  cov.Y2 = diag(rep(var.Y2,ncol(Gd2)))+Gd2
  
  sdY1 = sqrt(diag(cov.Y1))
  sdY2 = sqrt(diag(cov.Y2))
  
  RR = 1/(1+exp(-TR.para[15]))
  C.Y1Y2 = RR*sdY1*sdY2
  
  Cmat.Y1Y2=diag(C.Y1Y2)
  
  Cov.mat1 = cbind(cov.Y1,Cmat.Y1Y2)
  Cov.mat2 = cbind(t(Cmat.Y1Y2),cov.Y2)
  Cov.mat = rbind(Cov.mat1,Cov.mat2)
  
  N.obs= 2*length(Ey1.p)/nrow(Cov.mat)
  l.Y1=nrow(Cov.mat)/2
  
  Pred.Y =matrix(ncol=2,nrow=length(Ey1.p))
  for(i in 1:N.obs){
    StrtY = l.Y1*(i-1)+1
    EndY = l.Y1*i
    
    Ey =  c(Ey1.p[StrtY:EndY],Ey2.p[StrtY:EndY])
    Yi = rmvnorm(1,Ey,Cov.mat)
    Y.mat=matrix(exp(Yi),ncol=2)
    Pred.Y[StrtY:EndY,]=Y.mat
  }
  
  return(Pred.Y)
}