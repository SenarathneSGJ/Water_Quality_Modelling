NLog_Like=function(dataX,samp_no,dataY,para,dist.Y,z){
  
  Ey1=para[1]+para[2]*dataX[,1]+para[3]*dataX[,2]+para[4]*dataX[,3]+para[16]*z[,1]+para[17]*z[,2]+para[18]*z[,3]
  Ey2=para[5]+para[6]*dataX[,1]+para[7]*dataX[,2]+para[8]*dataX[,3]+para[19]*z[,1]+para[20]*z[,2]+para[21]*z[,3]
  
  r21 <- 1/(1+exp(-para[10]))
  r21[r21==0]=0.00001
  r11 <- exp(para[9])*r21
  r11[r11>exp(600)] <- exp(600)
  
  Gdist1 <- exp(-1*(dist.Y/r21))
  Gd1 <- matrix(ncol=ncol(Gdist1),nrow=nrow(Gdist1))
  Gd1[Gdist1!=1] <- r11*Gdist1[Gdist1!=1]
  Gd1[Gdist1==1] <- r11+r01
  
  r22 <- 1/(1+exp(-para[12]))
  r22[r22==0]=0.00001
  r12 <- exp(para[11])*r22
  r12[r12==Inf] <- exp(600)
  
  
  Gdist2 <- exp(-1*(dist.Y/r22))
  Gd2 <- matrix(ncol=ncol(Gdist2),nrow=nrow(Gdist2))
  Gd2[Gdist2!=1] <- r12*Gdist2[Gdist2!=1]
  Gd2[Gdist2==1] <- r12+r02
  
  var.Y1= exp(para[13])
  var.Y2= exp(para[14])
  
  cov.Y1 = diag(rep(var.Y1,ncol(Gd1)))+Gd1
  cov.Y2 = diag(rep(var.Y2,ncol(Gd2)))+Gd2
  
  sdY1 = sqrt(diag(cov.Y1))
  sdY2 = sqrt(diag(cov.Y2))
  
  RR = 1/(1+exp(-para[15]))
  C.Y1Y2 = RR*sdY1*sdY2
  
  Cmat.Y1Y2=diag(C.Y1Y2)
  
  Cov.mat1 = cbind(cov.Y1,Cmat.Y1Y2)
  Cov.mat2 = cbind(t(Cmat.Y1Y2),cov.Y2)
  Cov.mat = rbind(Cov.mat1,Cov.mat2)
  
  N.obs= ceiling(2*nrow(dataY)/nrow(Cov.mat))
  l.Y1=nrow(Cov.mat)/2
  
  Lk_Y =c()
  for(i in 1:N.obs){
    if(i==N.obs){
      StrtY = l.Y1*(i-1)+1
      EndY = nrow(dataY)
      
      Y = log(c(dataY[StrtY:EndY,1],dataY[StrtY:EndY,2]))
      Ey =  c(Ey1[StrtY:EndY],Ey2[StrtY:EndY])
      Lk_Y[i] = dmvnorm(Y,Ey,Cov.mat[c(samp_no,samp_no+22),c(samp_no,samp_no+22)],log=T)
      
    }else{
      StrtY = l.Y1*(i-1)+1
      EndY = l.Y1*i
      
      Y = log(c(dataY[StrtY:EndY,1],dataY[StrtY:EndY,2]))
      Ey =  c(Ey1[StrtY:EndY],Ey2[StrtY:EndY])
      Lk_Y[i] = dmvnorm(Y,Ey,Cov.mat,log=T)
    }
  }
  
  log_likelihood <- sum(Lk_Y)
  Neg_log_like <- -1*log_likelihood
  Neg_log_like[Neg_log_like>1e+100]=1e+100
  Neg_log_like[Neg_log_like< -1e+100]= -1e+100
  print(Neg_log_like)
  return(Neg_log_like)
  
}