library(MASS)
library(ggplot2)

f <- function(x){1+sin(20*x)-cos(30*x)}
xx <- seq(0.00001,1,length.out=300)

plot(xx,sapply(xx,f))

#Taken from Plumlee
phif =function(tv,p){
  lambdan = 4/(pi^2*(2*seq(0,p-1)+1)^2)
  3*sin(outer(tv,1/sqrt(lambdan)))*outer(rep(1,length(tv)),sqrt(lambdan))
}

xxdown <- seq(0.00001,1,length.out=ndown)
plot(xxdown,sapply(xxdown,f),ylab='f(t)',xlab='t')
lines(xx,sapply(xx,f))

p <- 30
phiptr <- phif(xxdown,p)
xxnew <- sort(c(xx,xxdown))
phiptest <- phif(xxnew,p)
puedoinvphitr<-solve(phiptr%*%t(phiptr))

fdata <- sapply(xxdown,f)

nsamps <- 1000
data <- mvrnorm(n = nsamps, mu= phiptest%*%t(phiptr)%*%puedoinvphitr%*%fdata , Sigma=phiptest%*%t(phiptest) -  phiptest%*%t(phiptr)%*%puedoinvphitr%*%phiptr%*%t(phiptest))

plot(xx,sapply(xx,f),type="l",xlab='x',ylab="f")
for(i in 1:100){
  lines(xxnew,data[i,],col="red",lty=2)
}
lines(xx,sapply(xx,f),type="l")
points(xxdown,sapply(xxdown,f),pch=16)


quantiles_data <- apply(data,2, quantile, probs = c(0.01,0.99))

ggplot() + geom_ribbon(aes(x=xxnew,ymin=quantiles_data[1,],ymax=quantiles_data[2,]),alpha=0.1, colour="red",fill="red")+
  geom_line(aes(x=xxnew,y=sapply(xxnew,f)), colour="black") + 
  geom_line(aes(xxnew,phiptest%*%t(phiptr)%*%puedoinvphitr%*%fdata),colour='blue') + 
  geom_point(aes(xxdown,sapply(xxdown,f)),size=2)+ 
  xlab("x") +
  ylab("f")



corr <- function(tv1,tv2){
  K1 <-matrix(tv1,nrow = length(tv1), ncol=length(tv2))
  K2 <-matrix(tv2,nrow = length(tv1), ncol=length(tv2),byrow=TRUE)
  K<-9/2*pmin(K1,K2)
  K
}

predictGP <-function(t_test, t_tr, f_tr){# obtain the fitted info
  gammahat <- 0
  
  # obtain the final covariance structure
  
  corr0 <- corr(t_tr, t_tr)  # c(t_tr, t_tr)
  kappa0 <- corr0 
  kappa0_inv <- solve(kappa0)
  kappa0_inv_f0 <- kappa0_inv%*%(as.vector(f_tr) - gammahat)
  corr_Tv <- corr(t_test, t_tr)# c(t_test, t_tr)
  kappa_Tv <- corr_Tv# c(t_test, t_tr)
  
  corr_toTv <- corr(t_tr, t_test)# c(t_tr, t_test)
  kappa_toTv <- corr_toTv#  c(., .)
  
  gamma_vec <- (gammahat + kappa_Tv%*%kappa0_inv_f0)
  kappa_Tv%*%kappa0_inv_f0
  gamma_vec
  corr_TvTv <- corr(t_test, t_test)# c(t_test, t_test)
  dkappa <- corr_TvTv
  kappa_vec <- pmax(dkappa -kappa_Tv%*%kappa0_inv%*%(kappa_toTv),10**-13)
  return(list(pred_mean=gamma_vec, pred_var=kappa_vec))}

t_tr = as.matrix(xxdown )
t_test = as.matrix(sort(c(t_tr,xx) ))
f_tr = sapply(xxdown,f)
predinfo <- predictGP(t_test,t_tr, f_tr)

nsamps <- 100000
data <- mvrnorm(n = nsamps, mu= predinfo$pred_mean , Sigma=predinfo$pred_var)
quantiles_data <- apply(data,2, quantile, probs = c(0.01,0.99))

pl<- ggplot() + 
  geom_ribbon(aes(x=t_test,ymin=quantiles_data[1,],ymax=quantiles_data[2,]),alpha=0.1, colour="red",fill="red") + 
  geom_line(aes(x=t_test,y=sapply(t_test,f)), colour="black") + 
  geom_line(aes(t_test,predinfo$pred_mean),colour='blue') + 
  geom_point(aes(xxdown,sapply(xxdown,f)),size=2)+ 
  xlab("x") +ylab("f")

pl

# approximations are not the same after the end points -- why?
# Gaussian process very inaccurate before the first end point and when 0 is included in the test set. 
# Baynesian Regression Problem is fine
# Gammahat does not seem matter as according to Plumlee. Not sure why.
#TODO: Try other covariance functions. How do the results change.
#TOOD: why is (9/2)*min(t,t') the covariance function of infinite sines
p <- 1000
phiptr <- phif(xxdown,p)
xxnew <- sort(c(xx,xxdown))
phiptest <- phif(xxnew,p)
puedoinvphitr<-solve(phiptr%*%t(phiptr))

fdata <- sapply(xxdown,f)

nsamps <- 100000
data <- mvrnorm(n = nsamps, mu= phiptest%*%t(phiptr)%*%puedoinvphitr%*%fdata , Sigma=phiptest%*%t(phiptest) -  phiptest%*%t(phiptr)%*%puedoinvphitr%*%phiptr%*%t(phiptest))
quantiles_data1 <- apply(data,2, quantile, probs = c(0.01,0.99))

pl + geom_ribbon(aes(x=xxnew,ymin=quantiles_data1[1,],ymax=quantiles_data1[2,]),alpha=0.2, colour="yellow",fill="yellow")

