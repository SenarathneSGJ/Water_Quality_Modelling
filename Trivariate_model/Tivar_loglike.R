NLog_Like=function(dataX,samp_no,dataY,para,dist.Y,z){
  
  Ey1=para[1]+para[2]*dataX[,1]+para[3]*dataX[,2]+para[4]*dataX[,3]+para[25]*z[,1]+para[26]*z[,2]+para[27]*z[,3]
  Ey2=para[5]+para[6]*dataX[,1]+para[7]*dataX[,2]+para[8]*dataX[,3]+para[28]*z[,1]+para[29]*z[,2]+para[30]*z[,3]
  Ey3=para[9]+para[10]*dataX[,1]+para[11]*dataX[,2]+para[12]*dataX[,3]+para[31]*z[,1]+para[32]*z[,2]+para[33]*z[,3]
  
  r21 <- 1/(1+exp(-para[14]))
  r21[r21==0]=0.00001
  r11 <- exp(para[13])*r21
  r11[r11>exp(600)] <- exp(600)
  
  Gdist1 <- exp(-1*(dist.Y/r21))
  Gd1 <- matrix(ncol=ncol(Gdist1),nrow=nrow(Gdist1))
  Gd1[Gdist1!=1] <- r11*Gdist1[Gdist1!=1]
  Gd1[Gdist1==1] <- r11+r01
  
  r22 <- 1/(1+exp(-para[16]))
  r22[r22==0]=0.00001
  r12 <- exp(para[15])*r22
  r12[r12==Inf] <- exp(600)
  
  Gdist2 <- exp(-1*(dist.Y/r22))
  Gd2 <- matrix(ncol=ncol(Gdist2),nrow=nrow(Gdist2))
  Gd2[Gdist2!=1] <- r12*Gdist2[Gdist2!=1]
  Gd2[Gdist2==1] <- r12+r02
  
  r23 <- 1/(1+exp(-para[18]))
  r23[r23==0]=0.00001
  r13 <- exp(para[17])*r23
  r13[r13==Inf] <- exp(600)
  
  Gdist3 <- exp(-1*(dist.Y/r23))
  Gd3 <- matrix(ncol=ncol(Gdist3),nrow=nrow(Gdist3))
  Gd3[Gdist3!=1] <- r13*Gdist3[Gdist3!=1]
  Gd3[Gdist3==1] <- r13+r03
  
  var.Y1= exp(para[19])
  var.Y2= exp(para[20])
  var.Y3= exp(para[21])
  
  cov.Y1 = diag(rep(var.Y1,ncol(Gd1)))+Gd1
  cov.Y2 = diag(rep(var.Y2,ncol(Gd2)))+Gd2
  cov.Y3 = diag(rep(var.Y3,ncol(Gd3)))+Gd3
  
  sdY1 = sqrt(diag(cov.Y1))
  sdY2 = sqrt(diag(cov.Y2))
  sdY3 = sqrt(diag(cov.Y3))
  
  RR12 = 1/(1+exp(-para[22]))
  C.Y1Y2 = RR12*sdY1*sdY2
  
  RR13 = 1/(1+exp(-para[23]))
  C.Y1Y3 = RR13*sdY1*sdY3
  
  RR23 = 1/(1+exp(-para[24]))
  C.Y2Y3 = RR23*sdY2*sdY3
  
  Cmat.Y1Y2=diag(C.Y1Y2)
  Cmat.Y1Y3=diag(C.Y1Y3)
  Cmat.Y2Y3=diag(C.Y2Y3)
  
  Cov.mat1 = cbind(cov.Y1,Cmat.Y1Y2,Cmat.Y1Y3)
  Cov.mat2 = cbind(t(Cmat.Y1Y2),cov.Y2,Cmat.Y2Y3)
  Cov.mat3 = cbind(t(Cmat.Y1Y3),t(Cmat.Y2Y3),cov.Y3)
  Cov.mat = rbind(Cov.mat1,Cov.mat2,Cov.mat3)
  
  N.obs= ceiling(3*nrow(dataY)/nrow(Cov.mat))
  l.Y1=nrow(Cov.mat)/3
  
  Lk_Y =c()
  for(i in 1:N.obs){
    if(i==N.obs){
      StrtY = l.Y1*(i-1)+1
      EndY = nrow(dataY)
      
      Y = log(c(dataY[StrtY:EndY,1],dataY[StrtY:EndY,2],dataY[StrtY:EndY,3]))
      Ey =  c(Ey1[StrtY:EndY],Ey2[StrtY:EndY],Ey3[StrtY:EndY])
      Lk_Y[i] = dmvnorm(Y,Ey,Cov.mat[c(samp_no,samp_no+22,samp_no+44),c(samp_no,samp_no+22,samp_no+44)],log=T)
      
    }else{
      StrtY = l.Y1*(i-1)+1
      EndY = l.Y1*i
      
      Y = log(c(dataY[StrtY:EndY,1],dataY[StrtY:EndY,2],dataY[StrtY:EndY,3]))
      Ey =  c(Ey1[StrtY:EndY],Ey2[StrtY:EndY],Ey3[StrtY:EndY])
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