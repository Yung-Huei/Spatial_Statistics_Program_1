library(spam);library(dotCall64);library(grid);library(maps);library(geoR)
library(fields);library(MASS);library(autoFRK);library(psych);library(SpatialTools)
load("0 02.RData")
truebeta <- c(1,1,1)
I <- diag(1,100); In <- diag(1,3)

beta_gls <- matrix(,n,3)
biasgls <- matrix(,n,3); vargls <- matrix(,n,3); msegls <- c()
#beta_ls<-matrix(,n,3)
#biasls <- matrix(,n,3); varls <- matrix(,n,3); msels <- c()

max_K <- matrix(,n,1)
K <- matrix(,n,2)
Kauto <- matrix(,n,3) #beta, K
biasKauto <- matrix(,n,3); varKauto <- matrix(,n,3); mseKauto <- c()
Kauto2 <- matrix(,n,3) #beta, K
biasKauto2 <- matrix(,n,3); varKauto2 <- matrix(,n,3); mseKauto2 <- c()
Kauto3 <- matrix(,n,3) #beta, K
biasKauto3 <- matrix(,n,3); varKauto3 <- matrix(,n,3); mseKauto3 <- c()
Kauto4 <- matrix(,n,3) #beta, K
biasKauto4 <- matrix(,n,3); varKauto4 <- matrix(,n,3); mseKauto4 <- c()
#Kauto5 <- matrix(,n,3) #beta, K
#biasKauto5 <- matrix(,n,3); varKauto5 <- matrix(,n,3); mseKauto5 <- c()

t1 <- proc.time()
u <- 1
for(u in 1:n){
  set.seed(u+11)
  #data------------
   bination <- cbind(grids,t(rbind(x1[u,],x2[u,],y[u,],z[u,])))
   Data <- rbind(bination[sample_100[u,1:100],])
   X <- cbind(1,Data[,3:4])

  #LSE------------
  # lmm<-lm(Data[,5]~Data[,3:4])
  # beta_ls[u,]<-lmm$coefficients
  # biasls[u,] <- beta_ls[u,]-truebeta

 #GLS------------
  #c(sigmasq[1],lambdavector[1]);c(4,p)
   gls <- likfit(data=Data[,5],coords=Data[,1:2],trend=trend.spatial(~Data[,3:4]),ini=c(4,p),
            fix.nugget=F,fix.kappa=F,cov.model="matern",lik.method ="ML",message=F)   
   beta_gls[u,] <- gls$beta
   Xbeta_gls <- gls$beta[1]+gls$beta[2]*X[,2]+gls$beta[3]*X[,3]     #X*hat_beta_gls
   
   biasgls[u,] <- beta_gls[u,]-truebeta
   msegls[u] <- sum((beta_gls[u,]-truebeta)^2)

   ytilde_w <- autoFRK(Data=Data[,5], loc=Data[,1:2], mu=Xbeta_gls)
   
#hat_K[u,4] <- ncol(ytilde_w$G)
   max_K[u,] <- ncol(ytilde_w$G)
    A <- function(par){
     Data_k <- matrix(,100,2)
     Data_k[,1] <- predict(autoFRK(Data=Data[,3], loc=Data[,1:2], mu=0, Kseq=par) )$pred.value
     Data_k[,2] <- predict(autoFRK(Data=Data[,4], loc=Data[,1:2], mu=0, Kseq=par) )$pred.value
     if(prod(Data_k[,1])==0 || prod(Data_k[,2])==0){
     return(c(Inf,Inf,Inf,Inf))
     }else{
     fit <- likfit(data=Data[,5],coords=Data[,1:2],trend=trend.spatial(~Data_k),ini=c(4,p),
            fix.nugget=F,fix.kappa=F,cov.model="matern",lik.method ="ML",message=F) 
     return(c( fit$beta,sqrt(sum((gls$beta - fit$beta)^2))))
     }
   }
oo <- sapply(3: max_K[u,], function(i){A(i)})
a <- quantile(oo[4,],seq(0.25,0.35,0.1))
K[u,1] <- min(which(oo[4,-(1:3)]<=a[1]))+5
K[u,2] <- min(which(oo[4,-(1:3)]<=a[2]))+5
#K[u,3] <- min(which(oo[4,-(1:3)]<=a[3]))+5
#K[u,4] <- min(which(oo[4,-(1:3)]<=a[4]))+5
#K[u,5] <- min(which(oo[4,-(1:3)]<=a[5]))+5

#plot(x=c(4:max_K[u,]),y=oo[4,1:93],type="b",col=8)
#a <- 3
#hat_Kvalue <- 0
#repeat{
#  diff <- matrix(,max_K[u,]-a-2,a)
#  for(j in 1:a){
#    if(j==a){
#    diff[,j] <- abs(oo[4,-(1:j)]-oo[4,-((ncol(oo)-j+1):ncol(oo))])
#    }else{
#    diff[,j] <- abs(oo[4,-(1:j)]-oo[4,-((ncol(oo)-j+1):ncol(oo))])[-(1:(a-j))]
#    }
#  }
#  o2 <- sapply(1:nrow(diff),function(j){sum(diff[j,])})
#  hat_K <- max( which(o2==min(o2)) )+2+a
#print(hat_K)
#  if(hat_K==hat_Kvalue){
#  break
#  }else{
#hat_Kvalue <- hat_K
#  a <- a+1
#print(a)
#  }
#}
#K[u,] <- hat_K-a
 #K=auto------------
   B <- function(par){

     x1_w <- autoFRK(Data=Data[,3], loc=Data[,1:2], mu=0, Kseq=par)
     x2_w <- autoFRK(Data=Data[,4], loc=Data[,1:2], mu=0, Kseq=par)
     ytilde_w <- autoFRK(Data=Data[,5], loc=Data[,1:2], mu=Xbeta_gls, Kseq=par)
     
     hat_x1 <- predict(x1_w)$pred.value
     hat_x2 <- predict(x2_w)$pred.value
     hat_X <- cbind(0,hat_x1,hat_x2)                                  #hat_X

     hat_ytilde_w <- predict(ytilde_w)$pred.value                     #hat_W*
     cov_ytilde_w <- ytilde_w$G %*% ytilde_w$M %*% t(ytilde_w$G)      #FMF
     Q <- solve(In - solve(t(X) %*% X) %*% t(X) %*% (hat_X))

     P <- X %*% solve(t(X) %*% X) %*% t(X)
     cov_rsr_y <- (I-P) %*% cov_ytilde_w %*% t(I-P) + ytilde_w$s*I
     b_rsr <- solve(t(X) %*% solve(cov_rsr_y) %*% X) %*% t(X) %*% solve(cov_rsr_y) %*% Data[,5]
     beta_FRK <- Q %*% (b_rsr - solve(t(X) %*% X) %*% t(X) %*% (hat_ytilde_w + (hat_X) %*% gls$beta)) 
     return( c(beta_FRK) )
   }
   Kauto[u,] <- B(K[u,1])
   biasKauto[u,] <- Kauto[u,1:3]-truebeta
   mseKauto[u] <- sum((Kauto[u,1:3]-truebeta)^2)  

   Kauto2[u,] <- B(K[u,2])
   biasKauto2[u,] <- Kauto2[u,1:3]-truebeta
   mseKauto2[u] <- sum((Kauto2[u,1:3]-truebeta)^2)

#   Kauto3[u,] <- B(K[u,3])
#   biasKauto3[u,] <- Kauto3[u,1:3]-truebeta
#   mseKauto3[u] <- sum((Kauto3[u,1:3]-truebeta)^2)

#   Kauto4[u,] <- B(K[u,4])
#   biasKauto4[u,] <- Kauto4[u,1:3]-truebeta
#   mseKauto4[u] <- sum((Kauto4[u,1:3]-truebeta)^2)

#   Kauto5[u,] <- B(K[u,5])
#   biasKauto5[u,] <- Kauto5[u,1:3]-truebeta
#   mseKauto5[u] <- sum((Kauto5[u,1:3]-truebeta)^2)
  #replace90
  T <- 20
  beta_FRK <- matrix(,T,3)
  I2 <- diag(1,90)
  for(i in 1:T){
   subData <- Data[sample(1:100,90),]
   subX <- cbind(1,subData[,3:4])
   gls <- likfit(data=subData[,5],coords=subData[,1:2],trend=trend.spatial(~subData[,3:4]),ini=c(4,p),
            fix.nugget=F,fix.kappa=F,cov.model="matern",lik.method ="ML",message=F)   
   Xbeta_gls <- gls$beta[1]+gls$beta[2]*subX[,2]+gls$beta[3]*subX[,3]     #X*hat_beta_gls
   ytilde_w <- autoFRK(Data=subData[,5], loc=subData[,1:2], mu=Xbeta_gls)
     x1_w <- autoFRK(Data=subData[,3], loc=subData[,1:2], mu=0, Kseq=ncol(ytilde_w$G))
     x2_w <- autoFRK(Data=subData[,4], loc=subData[,1:2], mu=0, Kseq=ncol(ytilde_w$G))
     
     hat_x1 <- predict(x1_w)$pred.value
     hat_x2 <- predict(x2_w)$pred.value
     hat_X <- cbind(0,hat_x1,hat_x2)                                  #hat_X

    hat_ytilde_w <- predict(ytilde_w)$pred.value                     #hat_W*
     cov_ytilde_w <- ytilde_w$G %*% ytilde_w$M %*% t(ytilde_w$G)      #FMF
     Q <- solve(In - solve(t(subX) %*% subX) %*% t(subX) %*% (hat_X))

     P <- subX %*% solve(t(subX) %*% subX) %*% t(subX)
     cov_rsr_y <- (I2-P) %*% cov_ytilde_w %*% t(I2-P) + ytilde_w$s*I2
     b_rsr <- solve(t(subX) %*% solve(cov_rsr_y) %*% subX) %*% t(subX) %*% solve(cov_rsr_y) %*% subData[,5]
     beta_FRK[i,] <- Q %*% (b_rsr - solve(t(subX) %*% subX) %*% t(subX) %*% (hat_ytilde_w + (hat_X) %*% gls$beta)) 
  }
   Kauto3[u,] <- sapply(1:3,function(j){mean(beta_FRK[,j])})
   biasKauto3[u,] <- Kauto3[u,1:3]-truebeta
   mseKauto3[u] <- sum((Kauto3[u,1:3]-truebeta)^2)

  T <- 50
  beta_FRK <- matrix(,T,3)
  I2 <- diag(1,90)
  for(i in 1:T){
   subData <- Data[sample(1:100,90),]
   subX <- cbind(1,subData[,3:4])
   gls <- likfit(data=subData[,5],coords=subData[,1:2],trend=trend.spatial(~subData[,3:4]),ini=c(4,p),
            fix.nugget=F,fix.kappa=F,cov.model="matern",lik.method ="ML",message=F)   
   Xbeta_gls <- gls$beta[1]+gls$beta[2]*subX[,2]+gls$beta[3]*subX[,3]     #X*hat_beta_gls
   ytilde_w <- autoFRK(Data=subData[,5], loc=subData[,1:2], mu=Xbeta_gls)
     x1_w <- autoFRK(Data=subData[,3], loc=subData[,1:2], mu=0, Kseq=ncol(ytilde_w$G))
     x2_w <- autoFRK(Data=subData[,4], loc=subData[,1:2], mu=0, Kseq=ncol(ytilde_w$G))
     
     hat_x1 <- predict(x1_w)$pred.value
     hat_x2 <- predict(x2_w)$pred.value
     hat_X <- cbind(0,hat_x1,hat_x2)                                  #hat_X

    hat_ytilde_w <- predict(ytilde_w)$pred.value                     #hat_W*
     cov_ytilde_w <- ytilde_w$G %*% ytilde_w$M %*% t(ytilde_w$G)      #FMF
     Q <- solve(In - solve(t(subX) %*% subX) %*% t(subX) %*% (hat_X))

     P <- subX %*% solve(t(subX) %*% subX) %*% t(subX)
     cov_rsr_y <- (I2-P) %*% cov_ytilde_w %*% t(I2-P) + ytilde_w$s*I2
     b_rsr <- solve(t(subX) %*% solve(cov_rsr_y) %*% subX) %*% t(subX) %*% solve(cov_rsr_y) %*% subData[,5]
     beta_FRK[i,] <- Q %*% (b_rsr - solve(t(subX) %*% subX) %*% t(subX) %*% (hat_ytilde_w + (hat_X) %*% gls$beta)) 
  }
   Kauto4[u,] <- sapply(1:3,function(j){mean(beta_FRK[,j])})
   biasKauto4[u,] <- Kauto4[u,1:3]-truebeta
   mseKauto4[u] <- sum((Kauto4[u,1:3]-truebeta)^2)
   print(u)
}
proc.time()-t1
w <- n


mean(K[1:w,1])
sd(K[1:w,1])/sqrt(w)
sapply(1:3, function(i){round(mean(Kauto[1:w,i]),4)})
sapply(1:3, function(i){round(sd(Kauto[1:w,i])/sqrt(w),4)})

sapply(1:3, function(i){round(mean(biasKauto[1:w,i]),4)})

v <- 1
for(v in 1:w){
  varKauto[v,] <- sapply(1:3, function(i){(Kauto[v,i]-mean(Kauto[1:w,i]))^2 })
}
vv<-sapply(1:3, function(i){round(mean(varKauto[1:w,i]),4)})
vv
TF <-sapply(1:w, function(v){
  L <- Kauto[v,1:3]-1.96*sqrt(vv)
  U <- Kauto[v,1:3]+1.96*sqrt(vv)
  L< c(beta0,beta1,beta2) & c(beta0,beta1,beta2)<U
})
table(TF[1,])/w
table(TF[2,])/w
table(TF[3,])/w
round(mean(mseKauto[1:w]),4)
round(sd(mseKauto[1:w])/sqrt(w),4)
#=============================================================================
mean(K[1:w,2])
sd(K[1:w,2])/sqrt(w)
sapply(1:3, function(i){round(mean(Kauto2[1:w,i]),4)})
sapply(1:3, function(i){round(sd(Kauto2[1:w,i])/sqrt(w),4)})

sapply(1:3, function(i){round(mean(biasKauto2[1:w,i]),4)})

v <- 1
for(v in 1:w){
  varKauto2[v,] <- sapply(1:3, function(i){(Kauto2[v,i]-mean(Kauto2[1:w,i]))^2 })
}
vv<-sapply(1:3, function(i){round(mean(varKauto2[1:w,i]),4)})
vv
TF <-sapply(1:w, function(v){
  L <- Kauto2[v,1:3]-1.96*sqrt(vv)
  U <- Kauto2[v,1:3]+1.96*sqrt(vv)
  L< c(beta0,beta1,beta2) & c(beta0,beta1,beta2)<U
})
table(TF[1,])/w
table(TF[2,])/w
table(TF[3,])/w
round(mean(mseKauto2[1:w]),4)
round(sd(mseKauto2[1:w])/sqrt(w),4)

#======================================================================================

#mean(K[1:w,3])
#sd(K[1:w,3])/sqrt(w)
sapply(1:3, function(i){round(mean(Kauto3[1:w,i]),4)})
sapply(1:3, function(i){round(sd(Kauto3[1:w,i])/sqrt(w),4)})

sapply(1:3, function(i){round(mean(biasKauto3[1:w,i]),4)})

v <- 1
for(v in 1:w){
  varKauto3[v,] <- sapply(1:3, function(i){(Kauto3[v,i]-mean(Kauto3[1:w,i]))^2 })
}
vv<-sapply(1:3, function(i){round(mean(varKauto3[1:w,i]),4)})
vv
TF <-sapply(1:w, function(v){
  L <- Kauto3[v,1:3]-1.96*sqrt(vv)
  U <- Kauto3[v,1:3]+1.96*sqrt(vv)
  L< c(beta0,beta1,beta2) & c(beta0,beta1,beta2)<U
})
table(TF[1,])/w
table(TF[2,])/w
table(TF[3,])/w
round(mean(mseKauto3[1:w]),4)
round(sd(mseKauto3[1:w])/sqrt(w),4)

#======================================================================================
#mean(K[1:w,4])
#sd(K[1:w,4])/sqrt(w)
sapply(1:3, function(i){round(mean(Kauto4[1:w,i]),4)})
sapply(1:3, function(i){round(sd(Kauto4[1:w,i])/sqrt(w),4)})

sapply(1:3, function(i){round(mean(biasKauto4[1:w,i]),4)})

v <- 1
for(v in 1:w){
  varKauto4[v,] <- sapply(1:3, function(i){(Kauto4[v,i]-mean(Kauto4[1:w,i]))^2 })
}
vv<-sapply(1:3, function(i){round(mean(varKauto4[1:w,i]),4)})
vv
TF <-sapply(1:w, function(v){
  L <- Kauto4[v,1:3]-1.96*sqrt(vv)
  U <- Kauto4[v,1:3]+1.96*sqrt(vv)
  L< c(beta0,beta1,beta2) & c(beta0,beta1,beta2)<U
})
table(TF[1,])/w
table(TF[2,])/w
table(TF[3,])/w
round(mean(mseKauto4[1:w]),4)
round(sd(mseKauto4[1:w])/sqrt(w),4)



#======================================================================================

#sapply(1:3, function(i){round(mean(beta_ls[1:w,i]),4)})
#sapply(1:3, function(i){round(sd(beta_ls[1:w,i])/sqrt(w),4)})

#sapply(1:3, function(i){round(mean(biasls[1:w,i]),4)})

#v <- 1
#for(v in 1:w){
#varls[v,] <- sapply(1:3, function(i){(beta_ls[v,i]-mean(beta_ls[1:w,i]))^2 })
#}
#vv<-sapply(1:3, function(i){round(mean(varls[1:w,i]),4)})
#vv
#TF <-sapply(1:w, function(v){
#L <- beta_ls[v,1:3]-1.96*sqrt(vv)
#U <- beta_ls[v,1:3]+1.96*sqrt(vv)
#L< c(beta0,beta1,beta2) & c(beta0,beta1,beta2)<U
#})
#table(TF[1,])/w
#table(TF[2,])/w
#table(TF[3,])/w

#round(mean(msels[1:w]),4)
#round(sd(msels[1:w])/sqrt(w),4)

#=============
sapply(1:3, function(i){round(mean(beta_gls[1:w,i]),4)})
sapply(1:3, function(i){round(sd(beta_gls[1:w,i])/sqrt(w),4)})

sapply(1:3, function(i){round(mean(biasgls[1:w,i]),4)})

v <- 1
for(v in 1:w){
vargls[v,] <- sapply(1:3, function(i){(beta_gls[v,i]-mean(beta_gls[1:w,i]))^2 })
}
vv<-sapply(1:3, function(i){round(mean(vargls[1:w,i]),4)})
vv
TF <-sapply(1:w, function(v){
L <- beta_gls[v,1:3]-1.96*sqrt(vv)
U <- beta_gls[v,1:3]+1.96*sqrt(vv)
L< c(beta0,beta1,beta2) & c(beta0,beta1,beta2)<U
})
table(TF[1,])/w
table(TF[2,])/w
table(TF[3,])/w

round(mean(msegls[1:w]),4)
round(sd(msegls[1:w])/sqrt(w),4)

#setEPS()
#postscript("081.eps")
#MSE <- c(mean(mseKauto[1:w]),mean(mseKauto2[1:w]),mean(mseKauto3[1:w]),mean(mseKauto4[1:w]),mean(mseKauto5[1:w]))
#MSEv <- c(sd(mseKauto[1:w]),sd(mseKauto2[1:w]),sd(mseKauto3[1:w]),sd(mseKauto4[1:w]),sd(mseKauto5[1:w]))
#L <- MSE-1.96*MSEv/sqrt(w)
#U <- MSE+1.96*MSEv/sqrt(w)
#plot(x=seq(0.15,0.55,0.1),y=MSE,xlab="quantile",type="b",lty=1,col=2,xaxt="n",ylim=c(0.46,0.60),main=expression(paste(rho," ","=0"," , ",phi," ","=0.2")))
#lines(x=seq(0.15,0.55,0.1),y=L,type="b",lty=2,col=4)
#lines(x=seq(0.15,0.55,0.1),y=U,type="b",lty=2,col=4)
#axis(1,seq(0.15,0.55,0.1))
#dev.off()
save.image("002result.RData")

