rm(list = ls());

library(dplyr);
# install.packages("senstrat")
library("lmtest")
library("sandwich")
library("xtable")
library("sensitivitymv")


#################
library(senstrat)

if.bb.ms=FALSE;
alpha = 0.05;
sigma = 0.5; 
beta = 0.3; 
 

# design parameter ... 
design.mat = rbind( c(1,1,1,1), 
                    c(1,2,1,2),
                    c(1,4,1,4),
                    c(1,2,2,4),
                    c(2,2,2,2),
                    #c(2,4,2,4),
                    c(4,4,4,4) ); 
colnames(design.mat) = c("00","01","10","11");

N = 1440;
K = 3; 
gammas = c(1,1.05, 1.1, 1.15); 
ngam = length(gammas);

nRep=1000;

res.mat = c(); 

nDesign = nrow(design.mat)

for (i in 1:nDesign){
  design = design.mat[i,];
  print(design)
  m = sum(design); 
  B  = N/m; 
  
  nPos.mat = matrix(0, nrow=length(gammas), ncol = K+1); 
  colnames(nPos.mat) = c("z1","z2","z3","combined")
  nPos.bb.ms.mat = matrix(0, nrow=length(gammas), ncol=2*(K+1)); 
  colnames(nPos.bb.ms.mat) = c("z1","z2","z3","combined", "z1.bb.ms","z2.bb.ms","z3","combined")
  
  for(rep in 1:nRep){
    if(rep%%100==0){print(rep)}; 
    
    z1 = rep( c(rep(0, design[1]+design[2]), rep(1, design[3]+design[4])), B); 
    z2 = rep( c(rep(0, design[1]), rep(1, design[2]), rep(0, design[3]), rep(1,design[4])), B); 
    # sanity.check
    # temp = cbind(z1,z2); temp[1:12,]
    eta = rnorm(N,mean=0,sd = sqrt(0.06))
    epsilon = rnorm(N,mean=0,sd = sigma)
    
    xi = .3*z1 + .25*z2 + eta
    z3 = sapply( xi, function(x) rbinom(1,1, max(0, min(1, x)))); 
    Y = beta*z3 + .1*z2 + epsilon;
    cor(Y,z1)
    cor(Y,z2)
    cor(Y,z3)
    
    st.bb = rep(1:B, each = m); 
    st.bb.3 = as.factor(paste(z1,z2,sep="_"));
    
    
    
    sc.bb.1 <- hodgeslehmann(y=Y, z=z1, st=st.bb, align="hl");
    sc.bb.2 <- hodgeslehmann(y=Y, z=z2, st=st.bb, align="hl");
    sc.bb.3 <- hodgeslehmann(y=Y, z=z3, st=st.bb.3, align="hl");
    
    if(!if.bb.ms){
      st.bb.ms.1 = as.factor(paste(st.bb, z2, sep="_"));
      st.bb.ms.2 = as.factor(paste(st.bb, z1, sep="_"));
      sc.bb.ms.1 <- hodgeslehmann(y=Y, z=z1, st=st.bb.ms.1, align="hl");
      sc.bb.ms.2 <- hodgeslehmann(y=Y, z=z2, st=st.bb.ms.2, align="hl");
      p.bb = matrix(NA, nrow=ngam, ncol = (K+1));
      colnames(p.bb) = c("z1", "z2", "z3", "c"); 
      rownames(p.bb) = gammas; 
      for(gam in gammas){
        p1.bb = senstrat(sc=sc.bb.1, z=z1, st=st.bb,   gamma = gam)$Result["P-value"];
        p2.bb = senstrat(sc=sc.bb.2, z=z2, st=st.bb,   gamma = gam)$Result["P-value"];
        p3.bb = senstrat(sc=sc.bb.3, z=z3, st=st.bb.3, gamma = gam)$Result["P-value"];
        
        pc.bb = truncatedP( sort(c(p1.bb,p2.bb,p3.bb))[2:3],trunc = 0.2); 
        p.bb[as.character(gam), ] = c(p1.bb,p2.bb,p3.bb, pc.bb)
      }
      pos = ifelse(p.bb <= alpha,1,0);
      # pos.combined = ifelse(apply(p.bb, 1, function(x) max(x)<= sqrt(alpha)),1,0);
      # nPos.mat = nPos.mat + cbind(pos,pos.combined)
      nPos.mat = nPos.mat + pos; 
    }else{
      st.bb.ms.1 = as.factor(paste(st.bb, z2, sep="_"));
      st.bb.ms.2 = as.factor(paste(st.bb, z1, sep="_"));
      sc.bb.ms.1 <- hodgeslehmann(y=Y, z=z1, st=st.bb.ms.1, align="hl");
      sc.bb.ms.2 <- hodgeslehmann(y=Y, z=z2, st=st.bb.ms.2, align="hl");
      p.bb = matrix(NA, nrow=ngam, ncol = 2*K-1);
      colnames(p.bb) = c("z1", "z2", "z3", "z1.bb.ms", "z2.bb.ms"); 
      rownames(p.bb) = gammas; 
      for(gam in gammas){
        p1.bb = senstrat(sc=sc.bb.1, z=z1, st=st.bb,   gamma = gam)$Result["P-value"];
        p2.bb = senstrat(sc=sc.bb.2, z=z2, st=st.bb,   gamma = gam)$Result["P-value"];
        p3.bb = senstrat(sc=sc.bb.3, z=z3, st=st.bb.3, gamma = gam)$Result["P-value"];
        
        p1.bb.ms = senstrat(sc=sc.bb.ms.1, z=z1, st=st.bb.ms.1,   gamma = gam)$Result["P-value"];
        p2.bb.ms = senstrat(sc=sc.bb.ms.2, z=z2, st=st.bb.ms.2,   gamma = gam)$Result["P-value"];
        p.bb[as.character(gam), ] = c(p1.bb,p2.bb,p3.bb, p1.bb.ms, p2.bb.ms)
      }
      
    }
  }
  
  rPos.mat = nPos.mat/nRep; 
  res.mat = cbind(res.mat, rPos.mat);
}
colnames(res.mat) = rep(c("Z1", "Z2", "D", "C"), nDesign)
rownames(res.mat) = gammas

fn = paste("power_analysis_sigma_", 10*sigma, "_beta", 10*beta,"_", nRep, "rep", sep="")

sink(
  paste( fn, ".txt", sep="")); 
xtable(round(res.mat,3)); 
sink()

save(res.mat, 
     file = paste( fn, ".RData", sep=""))

res.mat
