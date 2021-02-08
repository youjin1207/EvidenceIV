### library ###
library(MASS)
library(MatchIt)
library(senstrat)
library(sensitivitymv)
## data generating model
generateData.oneStratum.five <- function(seed = seed, # 
                         ## situation parameter
                         p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5,
                         delta= delta, #dependence between Z1 and Z2
                         rho=rho, # Unobserved covariate; rho \neq 0 --> unobserved confounding
                         var.epsilon = var.epsilon, var.eta = var.eta, 
                         # Response model 
                         alpha = alpha, beta = beta,
                         # IV on response --- potential violation of the exclusion restriction 
                         lambda1 = 0, lambda2 = 0, lambda3 = 0, lambda4 = 0, lambda5 = 0,
                         # IV on exposure --- weak/strong IV
                         phi1 = 0.2, phi2 = 0.2, phi3 = 0.2, phi4 = 0.2, phi5 = 0.2, nu=nu,
                         ## Design parameter
                         N = N # number of units in each stratum 
){
  
  set.seed(seed)
  #p2.1 = P(z2=1|z1=1); p2.0 = P(z2=1|z1=0)
  p2.0 = p2 - delta*p1; p2.1 = p2.0 + delta
  p3.0 = p3 - delta*p2; p3.1 = p3.0 + delta
  p4.0 = p4 - delta*p3; p4.1 = p4.0 + delta
  p5.0 = p5 - delta*p4; p5.1 = p5.0 + delta
  
  Z1 = rbinom(N, 1, p1); 
  Z2 = rep(NA,N); Z3 = rep(NA,N); Z4 = rep(NA,N); Z5 = rep(NA,N);
  n1.1 = sum(Z1); n1.0 = N - n1.1; Z2.1 = rbinom(n1.1,1,p2.1); Z2.0 = rbinom(n1.0,1,p2.0);
  Z2[Z1==1] = Z2.1;Z2[Z1==0] = Z2.0;

  n2.1 = sum(Z2); n2.0 = N - n2.1; Z3.1 = rbinom(n2.1,1,p3.1); Z3.0 = rbinom(n2.0,1,p3.0);
  Z3[Z2==1] = Z3.1; Z3[Z2==0] = Z3.0; 

  n3.1 = sum(Z3); n3.0 = N - n3.1; Z4.1 = rbinom(n3.1,1,p4.1); Z4.0 = rbinom(n3.0,1,p4.0);
  Z4[Z3==1] = Z4.1; Z4[Z3==0] = Z4.0;
  
  n4.1 = sum(Z4); n4.0 = N - n4.1; Z5.1 = rbinom(n4.1,1,p5.1); Z5.0 = rbinom(n4.0,1,p5.0);
  Z5[Z4==1] = Z5.1; Z5[Z4==0] = Z5.0;
  

  temp = mvrnorm(n=N, mu=c(0,0), Sigma = cbind(c(1,rho),c(rho,1)), tol = 1e-6, empirical = FALSE)
  epsilon = sqrt(var.epsilon)*temp[,1];
  eta = sqrt(var.eta)*temp[,2];
  
  xi = nu + phi1*Z1 + phi2*Z2 + phi3*Z3 + phi4*Z4 + phi5*Z5 + eta; 
  p = sapply(1:N, function(x) max(0, min(1,xi[x])));
  
  D = sapply(1:N, function(x) rbinom(1,1,p[x]) );
  
  R = alpha + lambda1*Z1 + lambda2*Z2 + lambda3*Z3 + lambda4*Z4 + lambda5*Z5 + beta*D  + epsilon;
  trt = apply(cbind(Z1,Z2,Z3,Z4,Z5),1,function(x) paste(x,collapse=""))
  data = cbind(Z1, Z2, Z3, Z4, Z5, D, R, p);
  rownames(data) = trt; 
  colnames(data) = c("Z1","Z2","Z3", "Z4", "Z5", "D", "R","p");
  return(list(data=data,N.trt = table(trt), N1 = sum(Z1))); 
}

## data generating model
ni = 200; 
I = 50; 
alpha = 0; 
nu= 0;
rho = 0.3; var.epsilon = 0.25; var.eta = 0.06;

##  depending on the beta values and lambda3 values
nested_five_beta00 = nested_five_beta01 = nested_five_beta02 = nested_five_beta03 = nested_five_beta04 = 
nested_five_beta00_lambda01 = nested_five_beta01_lambda01 = nested_five_beta02_lambda01 = 
nested_five_beta03_lambda01 = nested_five_beta04_lambda01 = list()

## for beta = 0.0, 0.1, 0.2, 0.3, 0.4
for(Rep in 1:1000){
  set.seed(Rep)
  beta = 0.0; # beta = 0.0, 0.1, 0.2, 0.3, 0.4
  deltas = seq(0.3, 0.7, 0.1)
  p.mutual1.nRep = p.mutual2.nRep =  p.mutual3.nRep  = 
   p.mutual4.nRep =  p.mutual5.nRep  = p.kruskal.nRep = rep(NA, length(deltas))
  for(d in 1:length(deltas)){
     p1 = 0.8; p2 = deltas[d]*p1; p3 = deltas[d]*p2; p4 = deltas[d]*p3; p5 = deltas[d]*p4;
     phi1 = phi2 = phi3 = phi4 = phi5 = 0.1
     lambda1 = lambda2 = lambda3 = lambda4 = lambda5 = 0
  
      data.total = list(); 
      N.trt.s.mat = matrix(0, nrow=I, ncol=4); 
      colnames(N.trt.s.mat) = trts; 
      rownames(N.trt.s.mat) = 1:I; 

      ni.1.s = rep(NA,I);
      names(ni.1.s) = 1:I; 
      Z1.ub = Z2.ub = c(); 
      R.ub = c(); 
      for(s in 1:I){
        temp = generateData.oneStratum.five(seed = Rep*I+s, 
                                       p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5,
                                       delta= deltas[d], 
                                       rho=rho, var.epsilon = var.epsilon,var.eta = var.eta, 
                                       alpha = alpha, beta = beta,
                                       lambda1 = 0, lambda2 = lambda2, lambda3 = lambda3, lambda4 = lambda4,
                                       phi1 = phi1, phi2 = phi2, phi3 = phi3, phi4 = phi4, phi5 = phi5, 
                                       nu=nu,
                                       N = ni); 
        Z1.ub = c(Z1.ub, temp$data[,"Z1"]);
        Z2.ub = c(Z2.ub, temp$data[,"Z2"]);
        R.ub = c(R.ub, temp$data[,"R"]); 
        data.total[[s]] = temp$data;
      }

      ## mutual stratificiation ##
      ###
      Ymat1 = Treat1 = c()
      index = 0
      for(s in 1:length(data.total)){
        dat = data.frame(data.total[[s]])
        match1 = tryCatch(matchit(Z1 ~ Z2, method = "exact", data = dat), error=function(err) NA); 
        if(!is.na( match1 )){
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
                        st= Ymat1[,2],gamma=1)
      p.mutual1.nRep[d] = res.mutual1$Result["P-value"];
      ###
      Ymat2 = Treat2 = c()
      index = 0
      for(s in 1:length(data.total)){
        dat = data.frame(data.total[[s]])
        match2 = tryCatch(matchit(Z2 ~ Z1 + Z3, method = "exact", data = dat), error=function(err) NA); 
        if(!is.na( match2 )){
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
    ###
    Ymat3 = Treat3 = c()
    index = 0
    for(s in 1:length(data.total)){
      dat = data.frame(data.total[[s]])   
      match3 = tryCatch(matchit(Z3 ~ Z2 + Z4, method = "exact", data = dat), error=function(err) NA); 
      if(!is.na( match3 )){
      nj = length(table(match3$subclass))
      for(j in 1:nj){
        index = index + 1
        treated = dat$R[which(match3$subclass == j & dat$Z3 == 1)]
        control = dat$R[which(match3$subclass == j & dat$Z3 == 0)]
        Ymat3 = rbind(Ymat3, cbind( c(treated, control), rep(index, sum(match3$subclass == j, na.rm = TRUE))))
        Treat3 = c(Treat3, c(rep(1,length(treated)), rep(0, length(control)) ) )
      }
    }
    }
    sampleSize.mutual1.nRep = nrow(Ymat3)
    sc.mutual3 = hodgeslehmann(y= Ymat3[,1], z= Treat3, 
                              st= Ymat3[,2], align="hl")
    res.mutual3 = senstrat(sc=sc.mutual3, z= Treat3,
                      st= Ymat3[,2],gamma=1)
    p.mutual3.nRep[d] = res.mutual3$Result["P-value"];
    ###
    Ymat4 = Treat4 = c()
    index = 0
    for(s in 1:length(data.total)){
      dat = data.frame(data.total[[s]])   
      match4 = tryCatch(matchit(Z4 ~ Z3 + Z5, method = "exact", data = dat), error=function(err) NA); 
      if(!is.na( match4 )){
      nj = length(table(match4$subclass))
      for(j in 1:nj){
        index = index + 1
        treated = dat$R[which(match4$subclass == j & dat$Z4 == 1)]
        control = dat$R[which(match4$subclass == j & dat$Z4 == 0)]
        Ymat4 = rbind(Ymat4, cbind( c(treated, control), rep(index, sum(match4$subclass == j, na.rm = TRUE))))
        Treat4 = c(Treat4, c(rep(1,length(treated)), rep(0, length(control)) ) )
      }
    }
    }
    sampleSize.mutual1.nRep = nrow(Ymat4)
    sc.mutual4 = hodgeslehmann(y= Ymat4[,1], z= Treat4, 
                              st= Ymat4[,2], align="hl")
    res.mutual4 = senstrat(sc=sc.mutual4, z= Treat4,
                      st= Ymat4[,2],gamma=1)
    p.mutual4.nRep[d] = res.mutual4$Result["P-value"];
    ###
    Ymat5 = Treat5 = c()
    index = 0
    for(s in 1:length(data.total)){
      dat = data.frame(data.total[[s]])   
      match5 = tryCatch(matchit(Z5 ~ Z4, method = "exact", data = dat), error=function(err) NA); 
      if(!is.na( match5 )){
      nj = length(table(match5$subclass))
      for(j in 1:nj){
        index = index + 1
        treated = dat$R[which(match5$subclass == j & dat$Z5 == 1)]
        control = dat$R[which(match5$subclass == j & dat$Z5 == 0)]
        Ymat5 = rbind(Ymat5, cbind( c(treated, control), rep(index, sum(match5$subclass == j, na.rm = TRUE))))
        Treat5 = c(Treat5, c(rep(1,length(treated)), rep(0, length(control)) ) )
      }
    }
    }
    sampleSize.mutual1.nRep = nrow(Ymat5)
    sc.mutual5 = hodgeslehmann(y= Ymat5[,1], z= Treat5, 
                              st= Ymat5[,2], align="hl")
    res.mutual5 = senstrat(sc=sc.mutual5, z= Treat5,
                      st= Ymat5[,2],gamma=1)
    p.mutual5.nRep[d] = res.mutual5$Result["P-value"];

    ## single instrument ##
    Ymat_star = Treat_star = c()
    index = 0
    for(s in 1:length(data.total)){
      dat = data.frame(data.total[[s]])
      dat$Zstar = ifelse(dat$Z1 == 0, 0, ifelse(dat$Z2 == 0, 1,
      ifelse(dat$Z3 == 0, 2, ifelse(dat$Z4 == 0, 3,  ifelse(dat$Z5 == 0, 4, 5)))))
      index = index + 1
      treated1 = dat$R[which(dat$Zstar == 1)]
      treated2 = dat$R[which(dat$Zstar == 2)]
      treated3 = dat$R[which(dat$Zstar == 3)]
      treated4 = dat$R[which(dat$Zstar == 4)]
      treated5 = dat$R[which(dat$Zstar == 5)]
      control = dat$R[which(dat$Zstar == 0)]
      Ymat_star = rbind(Ymat_star, cbind( c(treated1, treated2, treated3, treated4, treated5, control), rep(NA, nrow(dat))))
      Treat_star = c(Treat_star, c(rep(1,length(treated1)), rep(2, length(treated2)), 
        rep(3,length(treated3)), rep(4, length(treated4)), rep(5,length(treated5)), rep(0, length(control)) ) )
    }
    sampleSize.mutual1.nRep = nrow(Ymat_star)
    fit = kruskal.test(Ymat_star[,1] ~ Treat_star) #
    p.kruskal.nRep[d] = fit$p.value
  }

  nested_five_beta00[[ii]] = return(list(p.mutual1.nRep = p.mutual1.nRep,
        p.mutual2.nRep = p.mutual2.nRep,  p.mutual3.nRep = p.mutual3.nRep,  p.mutual4.nRep = p.mutual4.nRep,
        p.mutual5.nRep = p.mutual5.nRep,
        p.kruskal.nRep = p.kruskal.nRep))

}



## Figure 2 ## 
## beta = 0.0, 0.1, 0.2, 0.3, 0.4
dat = nested_five_beta00
pval.mutual1 = pval.mutual2 = pval.mutual3 = pval.mutual4 = 
  pval.mutual5 = pval.kruskal = matrix(NA, length(dat), 5)
for(i in 1:length(dat)){
  pval.mutual1[i,] = dat[[i]][[1]]
  pval.mutual2[i,] = dat[[i]][[2]]
  pval.mutual3[i,] = dat[[i]][[3]]
  pval.mutual4[i,] = dat[[i]][[4]]
  pval.mutual5[i,] = dat[[i]][[5]]
  pval.kruskal[i,] = dat[[i]][[6]]
}
five_beta00_lambda00_mutual = colMeans(pval.mutual1 <= 0.05^(1/5) & pval.mutual2 <= 0.05^(1/5) & 
           pval.mutual3 <= 0.05^(1/5) & pval.mutual4 <= 0.05^(1/5) & 
           pval.mutual5 <= 0.05^(1/5))
five_beta00_lambda00_kruskal = colMeans(pval.kruskal <= 0.05)
tmp.pvals = matrix(NA, nrow = length(dat), 5)
for(i in 1:length(dat)){
  for(j in 1:5){
    tmp = c(pval.mutual1[i,j], pval.mutual2[i,j],
            pval.mutual3[i,j], pval.mutual4[i,j],
            pval.mutual5[i,j])
    tmp = tmp[order(tmp)][1:5]
    tmp.pvals[i,j] = truncatedP(tmp, 1)
  }
}
five_beta00_lambda00_combine = colMeans(tmp.pvals <= 0.05)


dat = nested_five_beta00_lambda01
pval.mutual1 = pval.mutual2 = pval.mutual3 = pval.mutual4 = 
  pval.mutual5 = pval.kruskal = matrix(NA, length(dat), 5)
for(i in 1:length(dat)){
  pval.mutual1[i,] = dat[[i]][[1]]
  pval.mutual2[i,] = dat[[i]][[2]]
  pval.mutual3[i,] = dat[[i]][[3]]
  pval.mutual4[i,] = dat[[i]][[4]]
  pval.mutual5[i,] = dat[[i]][[5]]
  pval.kruskal[i,] = dat[[i]][[6]]
}
five_beta00_lambda01_mutual = colMeans(pval.mutual1 <= 0.05^(1/4) & pval.mutual2 <= 0.05^(1/4) & 
           pval.mutual3 <= 0.05^(1/4) & pval.mutual4 <= 0.05^(1/4) & 
           pval.mutual5 <= 0.05^(1/4))
five_beta00_lambda01_kruskal = colMeans(pval.kruskal <= 0.05)

tmp.pvals = matrix(NA, nrow = length(dat), 5)
for(i in 1:length(dat)){
  for(j in 1:5){
    tmp = c(pval.mutual1[i,j], pval.mutual2[i,j],
            pval.mutual3[i,j], pval.mutual4[i,j],
            pval.mutual5[i,j])
    tmp = tmp[order(tmp)][2:5]
    tmp.pvals[i,j] = truncatedP(tmp, 1)
  }
}
five_beta00_lambda01_combine = colMeans(tmp.pvals <= 0.05)



## power plot ## 
## lambda = 0.0, 0.1 and delta = 0.3, 0.5, 0.7
pdf("Figure/five_lambda00_delta01.pdf",  width = 7, height = 6)
par(mfrow = c(1,1),   mar = c(5,5,5,3),  cex.lab = 1.5, 
    cex.main = 2.0, cex.axis = 1.5, tcl = 0.5, lwd = 2)
beta.points = seq(0, 0.4, 0.1)
plot(beta.points, rep(NA, length(beta.points)),
     xlim = c(0, 0.44), ylim = c(0, 1),
     xlab = expression(beta^"*"), ylab = "Type I error/Power",
     main = expression(paste(delta, " = 0.1, ", bold(lambda), " = (0,0,1,0,0)" )))
points(beta.points, c(five_beta00_lambda01_combine[1], five_beta01_lambda01_combine[1], 
                      five_beta02_lambda01_combine[1], five_beta03_lambda01_combine[1],
                      five_beta04_lambda01_combine[1]),
       type = "b", cex = 1, col = "black", pch = 20, lty = 1, lwd = 2)
points(beta.points, c(five_beta00_lambda01_kruskal[1], five_beta01_lambda01_kruskal[1], 
                      five_beta02_lambda01_kruskal[1], five_beta03_lambda01_kruskal[1],
                      five_beta04_lambda01_kruskal[1]),
       type = "b", cex = 1, col = "coral2", pch = 20, lty = 2, lwd = 2)

abline(h = 0.05, lty = 2, col = "grey")
abline(v = 0.0, lty = 2, col = "grey")
legend(x = 0.15, y = 0.3, c(expression(paste("combine Stratification (", nu, "=5)" )), 
                            "Kruskal-Wallis Test"), 
       col = c("black", "coral2"),
       cex = 1.2, lty = c(1, 2), bty = "n", lwd = 2)
dev.off()

pdf("Figure/five_lambda01_delta01.pdf",  width = 7, height = 6)
par(mfrow = c(1,1),   mar = c(5,5,5,3),  cex.lab = 1.5, 
    cex.main = 2.0, cex.axis = 1.5, tcl = 0.5, lwd = 2)
beta.points = seq(0, 0.4, 0.1)
plot(beta.points, rep(NA, length(beta.points)),
     xlim = c(0, 0.44), ylim = c(0, 1),
     xlab = expression(beta^"*"), ylab = "Type I error/Power",
     main = expression(paste(delta, " = 0.1, ", bold(lambda), " = (0,0,1,0,0)" )))
points(beta.points, c(five_beta00_lambda01_combine[1], five_beta01_lambda01_combine[1], 
                      five_beta02_lambda01_combine[1], five_beta03_lambda01_combine[1],
                      five_beta04_lambda01_combine[1]),
       type = "b", cex = 1, col = "black", pch = 20, lty = 1, lwd = 2)
points(beta.points, c(five_beta00_lambda01_kruskal[1], five_beta01_lambda01_kruskal[1], 
                      five_beta02_lambda01_kruskal[1], five_beta03_lambda01_kruskal[1],
                      five_beta04_lambda01_kruskal[1]),
       type = "b", cex = 1, col = "coral2", pch = 20, lty = 2, lwd = 2)

abline(h = 0.05, lty = 2, col = "grey")
abline(v = 0.0, lty = 2, col = "grey")
legend(x = 0.15, y = 0.3, c(expression(paste("combine Stratification (", nu, "=4)" )), 
                            "Kruskal-Wallis Test"), 
       col = c("black", "coral2"),
       cex = 1.2, lty = c(1, 2), bty = "n", lwd = 2)
dev.off()