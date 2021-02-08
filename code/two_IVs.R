library(MASS)
library(MatchIt)
library(senstrat)
library(sensitivitymv)
## data generating model
generateData.oneStratum <- function(seed = seed, # for(rep in 1:nRep) seed = rep; 
                                    ## situation parameter
                                    p1 = p1, p2 = p2, delta= delta, #dependence between Z1 and Z2
                                    rho=rho, # Unobserved covariate; rho \neq 0 --> unobserved confounding
                                    var.epsilon = var.epsilon, var.eta = var.eta, 
                                    # Response model 
                                    alpha = alpha, beta = beta,
                                    # IV on response --- potential violation of the exclusion restriction 
                                    lambda1 = 0, lambda2 = 0.1,
                                    # IV on exposure --- weak/strong IV
                                    phi1 = .2, phi2 = .25, nu=nu,
                                    ## Design parameter
                                    N = N # number of units in each stratum 
){
  
  set.seed(seed)
  #p2.1 = P(z2=1|z1=1); p2.0 = P(z2=1|z1=0)
  p2.0 = p2 - delta*p1;
  p2.1 = p2.0 + delta
  
  # sanity check 
  #p2.0*(1-p1)+p2.1*p1 ==  p2;
  
  
  Z1 = rbinom(N, 1, p1); 
  n1.1 = sum(Z1);
  n1.0 = N - n1.1;
  
  Z2.1 = rbinom(n1.1,1,p2.1);
  Z2.0 = rbinom(n1.0,1,p2.0);
  
  Z2 = rep(NA,N);
  Z2[Z1==1] = Z2.1;
  Z2[Z1==0] = Z2.0; 
  
  # sanity check
  #r2 = mean(Z2); 
  temp = mvrnorm(n=N, mu=c(0,0), Sigma = cbind(c(1,rho),c(rho,1)), tol = 1e-6, empirical = FALSE)
  epsilon = sqrt(var.epsilon)*temp[,1];
  eta = sqrt(var.eta)*temp[,2];
  
  xi = nu + phi1*Z1 + phi2*Z2 + eta; 
  p = sapply(1:N, function(x) max(0, min(1,xi[x])));
  
  #print(p)
  Z3 = sapply(1:N, function(x) rbinom(1,1,p[x]) );
  
  R = alpha + lambda1*Z1 + lambda2*Z2 + beta*Z3+epsilon;
  trt = apply(cbind(Z1,Z2),1,function(x) paste(x,collapse=""))
  data = cbind(Z1, Z2, Z3, R, p);
  rownames(data) = trt; 
  colnames(data) = c("Z1","Z2","Z3","R","p");
  return(list(data=data,N.trt = table(trt), N1 = sum(Z1))); 
}

### parameter
ni = 50; 
I = 50; 
alpha = 0; beta = 0; 
nu = 0;
rho = 0.3; var.epsilon = 0.25; var.eta = 0.06;
##  depending on the beta values and lambda2 values
nested_beta00_lambda01 = nested_beta01_lambda01 = nested_beta02_lambda01 = nested_beta03_lambda01 = nested_beta04_lambda01 = 
nested_beta00_lambda00 = nested_beta01_lambda00 = nested_beta02_lambda00 = nested_beta03_lambda00 = nested_beta04_lambda00  = list()

## for beta = 0.0, 0.1, 0.2, 0.3, 0.4
for(Rep in 1:1000){
  set.seed(Rep)
  beta = 0.0 # beta = 0.0, 0.1, 0.2, 0.3, 0.4
  lambda2 = 0.1 # lambda2 = 0.0, 0.1
  deltas = seq(0.1, 0.9, 0.2)
  p.mutual1.nRep = p.mutual2.nRep = p.kruskal.nRep = rep(NA, length(deltas))
  
  for(d in 1:length(deltas)){
    p1 = 0.8; p2 = deltas[d]*p1 
    data.total = list(); 
    
    N.trt.s.mat = matrix(0, nrow=I, ncol=4); 
    colnames(N.trt.s.mat) = trts; 
    rownames(N.trt.s.mat) = 1:I; 
    
    ni.1.s = rep(NA,I);
    names(ni.1.s) = 1:I; 
    
    Z1.ub = Z2.ub = c(); 
    R.ub = c(); 
    for(s in 1:I){
      temp = generateData.oneStratum(seed = Rep*I+s, 
                                     p1 = p1, p2 = p2, delta= deltas[d], 
                                     rho=rho, var.epsilon = var.epsilon,var.eta = var.eta, 
                                     alpha = alpha, beta = beta,
                                     lambda1 = 0, lambda2 = lambda2,
                                     phi1 = 0.5, phi2 = 0.5, nu=nu,
                                     N = ni); 
      Z1.ub = c(Z1.ub, temp$data[,"Z1"]);
      Z2.ub = c(Z2.ub, temp$data[,"Z2"]);
      R.ub = c(R.ub, temp$data[,"R"]); 
      data.total[[s]] = temp$data;
      N.trt.s.mat[s,names(temp$N.trt)] = temp$N.trt; 
    }
      
    ## mutual stratification ##
    Ymat1 = Treat1 = c()
    index = 0
    for(s in 1:length(data.total)){
      dat = data.frame(data.total[[s]])
      if(sum(table(dat$Z1, dat$Z2)!=0) == 3){
      match1 = matchit(Z1 ~ Z2, method = "exact", data = dat)
      nj = length(table(match1$subclass))
      for(j in 1:nj){
        index = index + 1
        treated = dat$R[which(match1$subclass == j & dat$Z1 == 1)]
        control = dat$R[which(match1$subclass == j & dat$Z1 == 0)]
        Ymat1 = rbind(Ymat1, cbind( c(treated, control), rep(index, sum(match1$subclass == j, na.rm = TRUE))))
        Treat1 = c(Treat1, c(rep(1,length(treated)), rep(0, length(control)) ) )
      }
    }
    }
    sampleSize.mutual1.nRep = nrow(Ymat1)
    sc.mutual1 = hodgeslehmann(y= Ymat1[,1], z= Treat1, 
                              st= Ymat1[,2], align="hl")
    res.mutual1 = senstrat(sc=sc.mutual1, z= Treat1,
                      st= Ymat1[,2], gamma=1)
    p.mutual1.nRep[d] = res.mutual1$Result["P-value"];
    
    Ymat2 = Treat2 = c()
    index = 0
    for(s in 1:length(data.total)){
      dat = data.frame(data.total[[s]])
      if(sum(table(dat$Z1, dat$Z2)!=0) == 3){
      match2 = matchit(Z2 ~ Z1, method = "exact", data = dat)
      nj = length(table(match2$subclass))
      for(j in 1:nj){
        index = index + 1
        treated = dat$R[which(match2$subclass == j & dat$Z2 == 1)]
        control = dat$R[which(match2$subclass == j & dat$Z2 == 0)]
        Ymat2 = rbind(Ymat2, cbind( c(treated, control), rep(index, sum(match2$subclass == j, na.rm = TRUE))))
        Treat2 = c(Treat2, c(rep(1,length(treated)), rep(0, length(control)) ) )
      }
    }
    }
    sampleSize.mutual1.nRep = nrow(Ymat2)
    sc.mutual2 = hodgeslehmann(y= Ymat2[,1], z= Treat2, 
                              st= Ymat2[,2], align="hl")
    res.mutual2 = senstrat(sc=sc.mutual2, z= Treat2,
                      st= Ymat2[,2],gamma=1)
    p.mutual2.nRep[d] = res.mutual2$Result["P-value"];


    ## single instrument ##
    Ymat_star = Treat_star = c()
    index = 0
    for(s in 1:length(data.total)){
      dat = data.frame(data.total[[s]])
      dat$Zstar = ifelse(dat$Z1 == 0, 0, ifelse(dat$Z2 == 1, 2, 1))
      index = index + 1
      treated1 = dat$R[which(dat$Zstar == 1)]
      treated2 = dat$R[which(dat$Zstar == 2)]
      control = dat$R[which(dat$Zstar == 0)]
      Ymat_star = rbind(Ymat_star, cbind( c(treated1, treated2, control), rep(NA, nrow(dat))))
      Treat_star = c(Treat_star, c(rep(1,length(treated1)), rep(2, length(treated2)), rep(0, length(control)) ) )
  
    }
    sampleSize.mutual1.nRep = nrow(Ymat_star)
    fit = kruskal.test(Ymat_star[,1] ~ Treat_star) # no stratification
    p.kruskal.nRep[d] = fit$p.value

  nested_beta00_lambda01[[ii]] = list(p.mutual1.nRep = p.mutual1.nRep,
        p.mutual2.nRep = p.mutual2.nRep, p.kruskal.nRep = p.kruskal.nRep)
  }
}

### Figure 1 ####
## read beta = 0.0, 0.1, 0.2, 0.3, 0.4.
dat = nested_beta00_lambda00
pval.mutual1 = pval.mutual2 = pval.kruskal = matrix(NA, length(dat), 5)
for(i in 1:length(dat)){
  pval.mutual1[i,] = dat[[i]][[1]]
  pval.mutual2[i,] = dat[[i]][[2]]
  pval.kruskal[i,] = dat[[i]][[3]]
}
beta00_lambda00_mutual1 = colMeans(pval.mutual1 <= 0.05)
beta00_lambda00_mutual2 = colMeans(pval.mutual2 <= 0.05)
beta00_lambda00_mutual = colMeans(pval.mutual1 <= 0.05^(1/2) & pval.mutual2 <= 0.05^(1/2))
beta00_lambda00_kruskal = colMeans(pval.kruskal <= 0.05)

dat = nested_beta00_lambda01
pval.mutual1 = pval.mutual2 = pval.kruskal = matrix(NA, length(dat), 5)
for(i in 1:length(dat)){
  pval.mutual1[i,] = dat[[i]][[1]]
  pval.mutual2[i,] = dat[[i]][[2]]
  pval.kruskal[i,] = dat[[i]][[3]]
}
beta00_lambda01_mutual1 = colMeans(pval.mutual1 <= 0.05)
beta00_lambda01_mutual2 = colMeans(pval.mutual2 <= 0.05)
beta00_lambda01_mutual = colMeans(pval.mutual1 <= 0.05 & pval.mutual2 <= 0.05)
beta00_lambda01_kruskal = colMeans(pval.kruskal <= 0.05)

tmp.pvals = matrix(NA, nrow = length(dat), 5)
for(i in 1:length(dat)){
  for(j in 1:5){
    tmp.pvals[i,j] = truncatedP(c(pval.mutual1[i,j], pval.mutual2[i,j]), 1)
  }
}
beta00_lambda01_combine = colMeans(tmp.pvals <= 0.05)

## power plot ## 
## lambda = 0.0, 0.1 and delta = 0.3, 0.5, 0.7
pdf("Figure/lambda00_delta03.pdf",  width = 7, height = 6)
par(mfrow = c(1,1),   mar = c(5,5,5,3),  cex.lab = 1.5, 
    cex.main = 2.0, cex.axis = 1.5, tcl = 0.5, lwd = 2)
beta.points = seq(0, 0.4, 0.1)
plot(beta.points, rep(NA, length(beta.points)),
     xlim = c(0, 0.44), ylim = c(0, 1),
     xlab = expression(beta^"*"), ylab = "Type I error/Power",
     main = expression(paste(delta, " = 0.3, ", lambda[2], " = 0.0" )))
## s_bar = 2, gamma = 1
points(beta.points, c(beta00_lambda00_combine[2], beta01_lambda00_combine[2], 
                      beta02_lambda00_combine[2],
                      beta03_lambda00_combine[2],
                      1),
       type = "b", cex = 1, col = "black", pch = 20, lty = 1, lwd = 2)
## s_bar = 2, gamma = 1.05
points(beta.points, c(beta00_lambda00_kruskal[2], beta01_lambda00_kruskal[2], 
                      beta02_lambda00_kruskal[2], beta03_lambda00_kruskal[2],
                      beta04_lambda00_kruskal[2]),
       type = "b", cex = 1, col = "coral2", pch = 20, lty = 2, lwd = 2)

abline(h = 0.05, lty = 2, col = "grey")
abline(v = 0.0, lty = 2, col = "grey")
#legend(x = 0.15, y = 0.3, c(expression(paste("Mutual Stratification (", nu, "=2)" )), 
#                           "Kruskal-Wallis Test"), 
#       col = c("black", "coral2"),
#       cex = 1.3, lty = c(1, 2), bty = "n", lwd = 2)
dev.off()
