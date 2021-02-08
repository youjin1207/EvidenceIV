rm(list = ls());
# install.packages("ivreg", dependencies = TRUE)
library("MASS")
library(dplyr)
library("ivreg")
library(senstrat)
trts = c("11","10","01","00");
trt.mat = rbind(c(1,1,0,0), c(1,0,1,0));
rownames(trt.mat) = c("Z1","Z2")
colnames(trt.mat)=trts; 
nTrt = length(trts);

nRep=1000;
ni = 50; 
I = 50; 

n = ni*I;
p2 = .4; 
delta = .14; 
lambda2 = .1; 
epsilon = 1 

p1.0 = 1/3 - delta*p2
p1.1 = p1.0 +delta;
block = rep(1:I, each = ni);
 

p.mat = c(); 

for(Rep in 1:nRep){
  if(Rep %%10==0) print(Rep)
  
  Z2 = rbinom(n, 1, p2);
  Z1 = rep(NA,n); 
  eta = rnorm(n, sd=sqrt(.06));
  n.1 = sum(Z2==1);
  
  
  Z1[Z2==1] = rbinom(n.1, 1, p1.1);
  Z1[Z2==0] = rbinom(n - n.1, 1, p1.0);
  D = sapply( .3*Z1+.25*Z2 + eta, function(x) rbinom(1,1,max(0,min(1,x)))); 
  
  trt =paste(Z1,Z2,sep="") 
  
  R = rnorm(n, mean=0,sd=epsilon) + lambda2*Z2; 
  
  dat = data.frame(Z1,Z2,D,R,block, trt);
  
  lm.2sls = summary(ivreg(R ~ D | Z1 + Z2, data = dat)); 
  p.2sls = (lm.2sls$coefficients)[2,4];
  
  ## Reinforced design (RD) 
  sc.1 <-hodgeslehmann(y=dat$R,z=dat$Z1,st=dat$block, align="hl")
  res.1  = senstrat(sc=sc.1, z=dat$Z1, st=dat$block, gamma=1);
  p.1  = res.1$Result["P-value"];
  
  
  sc.2  <-hodgeslehmann(y=dat$R,z=dat$Z2,st=paste(dat$block,dat$Z1,sep=""), align="hl")
  res.2   = senstrat(sc=sc.2 , z=dat$Z2, st=paste(dat$block,dat$Z1,sep=""), gamma=1);
  p.2   = res.2$Result["P-value"];

  sc.D   <-hodgeslehmann(y=dat$R,z=dat$D,st=paste(dat$block,dat$Z1, dat$Z2,sep=""), align="hl")
  res.D    = senstrat(sc=sc.D  , z=dat$D, st=paste(dat$block,dat$Z1, dat$Z2,sep=""), gamma=1);
  p.D    = res.D $Result["P-value"];
  
  p.rd = c(p.1, p.2, p.D);
  
  # Balanced block design (BB) 
  
  dat.b = c(); 
  
  for( i in 1:I){
    temp= dat[dat$block==i,];
    nb = min(table(temp$trt))
    if(nb > 1){
      for(t  in trts){
        dat.b = rbind(dat.b, temp[temp$trt==t,][1:nb,])
      }
    }else{
      if(nb == 1){
        for(t  in trts){
          nb.t = sum(temp$trt==t); 
          if(nb.t>1){
            dat.b = rbind(dat.b, temp[temp$trt==t,][1:nb,])
          }else{
            dat.b = rbind(dat.b, temp[temp$trt==t,])
          }
        }
      }
    }
    
  }
  
  sc.1 <-hodgeslehmann(y=dat.b$R,z=dat.b$Z1,st=dat.b$block, align="hl")
  res.1  = senstrat(sc=sc.1, z=dat.b$Z1, st=dat.b$block, gamma=1);
  p.1  = res.1$Result["P-value"];
  
  
  sc.2  <-hodgeslehmann(y=dat.b$R,z=dat.b$Z2,st=paste(dat.b$block,dat.b$Z1,sep=""), align="hl")
  res.2   = senstrat(sc=sc.2 , z=dat.b$Z2, st=paste(dat.b$block,dat.b$Z1,sep=""), gamma=1);
  p.2   = res.2$Result["P-value"];
  
  sc.D   <-hodgeslehmann(y=dat.b$R,z=dat.b$D,st=paste(dat.b$block,dat.b$Z1, dat.b$Z2,sep=""), align="hl")
  res.D    = senstrat(sc=sc.D  , z=dat.b$D, st=paste(dat.b$block,dat.b$Z1, dat.b$Z2,sep=""), gamma=1);
  p.D    = res.D $Result["P-value"];
  
  p.bb = c(p.1, p.2, p.D);
  
  p.mat = rbind(p.mat, c(p.2sls, p.rd, p.bb))

}

res = c(lambda2, epsilon, n, ni, i, apply(p.mat, 2, function(x) mean(x<0.05))); 
names(res) = c("lambda2", "epsilon", "n", "ni", "i", "p.2sls", "p.rd.1", "p.rd.2", "p.rd.D", "p.bb.1", "p.bb.2", "p.bb.D" )

save(res, file = paste("lambda2_", 10*lambda2, "_ep_", 10*epsilon, ".RData", sep="")) 

