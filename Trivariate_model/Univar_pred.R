pred2=function(TR.para,Fac.loc,dist.mat,z_new){
  
  Ey1.p = TR.para[1]+TR.para[2]*Fac.loc[,1]+TR.para[3]*Fac.loc[,2]+TR.para[4]*Fac.loc[,3]+TR.para[8]*z_new[,1]+TR.para[9]*z_new[,2]+TR.para[10]*z_new[,3]
  
  r21 <- 1/(1+exp(-TR.para[6]))
  r11 <- exp(TR.para[5])*r21
  
  Gdist1 <- exp(-1*(dist.mat/r21))
  Gd1 <- matrix(ncol=ncol(Gdist1),nrow=nrow(Gdist1))
  Gd1[Gdist1!=1] <- r11*Gdist1[Gdist1!=1]
  Gd1[Gdist1==1] <- r11+r01
  
  var.Y1= exp(TR.para[7])
  
  cov.Y1 = diag(rep(var.Y1,ncol(Gd1)))+Gd1
  
  N.obs= length(Ey1.p)/nrow(cov.Y1)
  l.Y1=nrow(cov.Y1)
  
  Pred.Y =matrix(ncol=1,nrow=length(Ey1.p))
  for(i in 1:N.obs){
    StrtY = l.Y1*(i-1)+1
    EndY = l.Y1*i
    
    Ey =  c(Ey1.p[StrtY:EndY])
    Yi = rmvnorm(1,Ey,cov.Y1)
    Y.mat=exp(Yi)
    Pred.Y[StrtY:EndY,]=Y.mat
  }
  
  return(Pred.Y)
}