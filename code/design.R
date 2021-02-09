rm(list = ls());
library(MASS)
library(dplyr)
library(senstrat)
#==========================
#------------------------------
# Bikram 2020
#------------------------------
VEwilcoxon<-function (y, z, st, tau = 0) { 
  ymiss<-is.na(y)
  y<-y[!ymiss]
  z<-z[!ymiss]
  st<-st[!ymiss]
  if (is.factor(st))   st <- as.integer(st)
  ust <- sort(unique(st))
  nst <- length(ust)
  if (tau != 0)     y <- y - tau * z
  sc <- rep(NA, length(y))
  for (i in 1:nst) {
    who <- (st == ust[i])
    yi <- y[who]
    vi <- rank(yi)/(length(yi)+1)
    sc[who] <- vi
  }
  list(sc=sc,z=z,st=st)
} 

#==========================
#------------------------------
# Generate data from null model
#------------------------------
generateData.fixedZ1Z2.list <- function(var.epsilon = var.epsilon,#var.eta = var.eta, 
                                        # Response model 
                                        #alpha = alpha, beta = beta,
                                        # IV on response --- potential violation of the exclusion restriction 
                                        lambda1 = 0, lambda2 = 0.1,
                                        # IV on exposure --- if weak?
                                        #phi1 = .2, phi2 = .25, nu=nu,
                                        ## Design parameter
                                        N.trt = N.trt, 
                                        trts = c("11","10","01","00"),
                                        trt.mat = trt.mat
){
  if(length(N.trt)==1){
    N.trt=rep(N.trt,4);
    names(N.trt) = trts; 
  }
  
  data.list = list(); 
  for (trt in trts){
    #EE = generateEE(N=N.trt[trt], rho=rho, var.epsilon = var.epsilon, var.eta = var.eta); 
    epsilon = rnorm(N.trt[trt],mean=0,sd = sqrt(var.epsilon))
    # eta = EE[,"eta"];
    # xi = nu + phi1*trt.mat["Z1",trt] + phi2*trt.mat["Z2",trt] + eta
    # p = sapply(1:N.trt[trt], function(x) max(0, min(1,xi[x])));
    # Z3 = sapply(1:N.trt[trt], function(x) rbinom(1,1,p[x]) );
    
    #R = alpha + lambda1*trt.mat["Z1",trt] + lambda2*trt.mat["Z2",trt] + beta*Z3+epsilon;
    R = lambda2*trt.mat["Z2",trt] + epsilon
    
    data = cbind(trt.mat["Z1",trt], trt.mat["Z2",trt], R);
    colnames(data) = c("Z1","Z2","R");
    data.list[[trt]] = data
  }
  
  return(data.list)
}

#################
trts = c("11","10","01","00");
trt.mat = rbind(c(1,1,0,0), c(1,0,1,0));
rownames(trt.mat) = c("Z1","Z2")
colnames(trt.mat)=trts; 
nTrt = length(trts);

# design parameter ... 
design.mat = rbind( c(1,1,1,1), 
                    c(1,2,1,2),
                    c(1,4,1,4),
                    c(1,2,2,4),
                    c(2,2,2,2),
                    #c(2,4,2,4),
                    c(4,4,4,4),
                    #c(2,8,2,8),
                    #c(2,4,4,8),
                    c(1,1,1,2),
                    c(1,1,1,3)); 


N = 3600;
nRep=10^4; #^4; 
# if(N < 500){
#   nRep = 10^4;
# }else{
#   nRep = 10^3
# }

nDesign = nrow(design.mat); 
design.id.vec = apply(design.mat, 1, function(x) paste(x,collapse ="-"));
rownames(design.mat) = design.id.vec; 

ni.vec = apply(design.mat, 1, sum);
nBlock.vec = N/ni.vec; 
block.id.mat = apply( cbind(nBlock.vec, ni.vec), 1, function(x) rep(1:x[1], x[2]) );
colnames(block.id.mat)= design.id.vec; 

temp = apply(design.mat,1,function(x) x/sum(x));
N.trt.mat = t(temp*N);
colnames(N.trt.mat)= trts; 
rownames(N.trt.mat) = design.id.vec; 

N.trt.max = apply(N.trt.mat,2,max)



lambda1 = 0; lambda2 = 0.1; var.epsilon = 0; #var.eta = 0.06;

p.nRep.hl =  matrix(NA, nrow=nRep, ncol = nDesign);
colnames(p.nRep.hl) = design.id.vec; 

p.nRep.ve = p.nRep.hl;


for(Rep in 1:nRep){
  print(Rep); 
  data.max.list = generateData.fixedZ1Z2.list(var.epsilon = var.epsilon, #,var.eta = var.eta, 
                                              # Response model 
                                              #alpha = alpha, beta = beta,
                                              # IV on response --- potential violation of the exclusion restriction 
                                              lambda1 = lambda1, lambda2 = lambda2,
                                              # IV on exposure --- if weak?
                                              #phi1 = .2, phi2 = .25, nu=nu,
                                              ## Design parameter
                                              N.trt = N.trt.max,
                                              trts = c("11","10","01","00"), trt.mat=trt.mat 
  );
  
  
  for ( design in design.id.vec){
    N.trt = N.trt.mat[design,];
    
    ## generate the data matrix arranged in the order of 00-01-10-11
    data = c();
    for (trt in trts){
      data = rbind(data, data.max.list[[trt]][1:N.trt[trt],]);
    }
    block = block.id.mat[,design];
    
    sc.hl <-hodgeslehmann(y=data[,"R"],z=data[,"Z1"],st=block, align="hl")
    res.hl = senstrat(sc=sc.hl, z=data[,"Z1"],st=block, gamma=1)
    p.nRep.hl[Rep, design] = res.hl$Result["P-value"];
    
    sc.ve <-VEwilcoxon(y=data[,"R"],z=data[,"Z1"],st=block, tau = 0)$sc
    res.ve = senstrat(sc=sc.ve, z=data[,"Z1"],st=block, gamma=1)
    p.nRep.ve[Rep, design] = res.ve$Result["P-value"];
  }
}

print(var.epsilon)
print(lambda2)
write.table(
  p.nRep.hl,paste("asymValidityBBD_lambda2_", lambda2, "_10sigma",10*sqrt(var.epsilon), "_N", N, "_hl_nRep",Rep,".txt",sep=""),quote=F,sep="\t", row.names=F, col.names=T)
write.table(
  p.nRep.ve,paste("asymValidityBBD_lambda2_", lambda2, "_10sigma",10*sqrt(var.epsilon), "_N", N, "_ve_nRep",Rep,".txt",sep=""),quote=F,sep="\t", row.names=F, col.names=T)

apply(p.nRep.hl, 2, function(x) mean(x<0.05)); 
apply(p.nRep.ve, 2, function(x) mean(x<0.05))
