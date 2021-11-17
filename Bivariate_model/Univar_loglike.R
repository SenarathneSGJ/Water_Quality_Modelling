NLog_Like_Res1=function(dataX,samp_no,dataY,para,dist.Y,z){
  
  Ey1=para[1]+para[2]*dataX[,1]+para[3]*dataX[,2]+para[4]*dataX[,3]+para[8]*z[,1]+para[9]*z[,2]+para[10]*z[,3]
  
  r21 <- 1/(1+exp(-para[6]))
  r21[r21==0]=0.00001
  r11 <- exp(para[5])*r21
  r11[r11>exp(600)] <- exp(600)
  
  Gdist1 <- exp(-1*(dist.Y/r21))
  Gd1 <- matrix(ncol=ncol(Gdist1),nrow=nrow(Gdist1))
  Gd1[Gdist1!=1] <- r11*Gdist1[Gdist1!=1]
  Gd1[Gdist1==1] <- r11+r01
  
  var.Y1= exp(para[7])
  Cov.mat = diag(rep(var.Y1,ncol(Gd1)))+Gd1
  
  N.obs= ceiling(nrow(dataY)/nrow(Cov.mat))
  l.Y1=nrow(Cov.mat)
  
  Lk_Y =c()
  for(i in 1:N.obs){
    if(i==N.obs){
      StrtY = l.Y1*(i-1)+1
      EndY = nrow(dataY)
      
      Y = log(c(dataY[StrtY:EndY,1]))
      Ey =  c(Ey1[StrtY:EndY])
      Lk_Y[i] = dmvnorm(Y,Ey,Cov.mat[samp_no,samp_no],log=T)
      
    }else{
      StrtY = l.Y1*(i-1)+1
      EndY = l.Y1*i
      
      Y = log(c(dataY[StrtY:EndY,1]))
      Ey =  c(Ey1[StrtY:EndY])
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