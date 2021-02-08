# Evidence factors from multiple, possibly invalid, instrumental variables

### Information

- Author: Anqi Zhao,Youjin Lee, Dylan Small, and Bikram Karmakar
- Required packages: MASS, dplyr, ivreg, senstrat, MatchIt, lmtest, sandwich, xtable, sensitivitymv


### Code

- `code/design.R`: This code investigates the effect of different design ratios.

- `code/BB_RD_2SLS.R`: This code investigates the validity of balanced block design in comparison with the reinforced analysis and two stage least squares regression.

- `code/power_analysis.R`: 


- `code/two_IVs.R`: This code generates two nested instruments and implements the mutual stratification and the Kruskal-Wallis test. This replicates Figure 1.

- `code/five_IVs.R`: This code generates five nested instruments and implements the mutual stratification and the Kruskal-Wallis test. This replicates Figure 2.

- `code/read_malaria.R`: This code analyzes the effect of malaria on stunting among children using the Asembo Bay Cohort study. This can replicate Figure 3 and Table 7 given the real dataset.

### Instructions for the use of sample data

In `Data/sample.csv`, we provide a hypothetical data with the same data structure as the Asembo Bay Cohort data. The data is provided for illustrative purpose only.


```{r}
## read the sample data
subdat = read.csv("Data/sample.csv", header = TRUE, sep = ",")

## associations betwen each instrument and the exposure variable
fitD = lm(D ~ Z1 + Z2 + Z3 + age + sickle + as.factor(village), data = subdat)

## implement the mutual stratification ##
## with IV 1
match = tryCatch(matchit(Z1 ~ Z2 + age + sickle + village, data = subdat, method = "exact"), 
                 error=function(err) NA); 
Ymat = Treat = c()
nj = length(table(match$subclass))
index = 0
for(j in 1:nj){
  index = index + 1
  treated = subdat$Y[which(match$subclass == j & subdat$Z1 == 1)]
  control = subdat$Y[which(match$subclass == j & subdat$Z1 == 0)]
  Ymat = rbind(Ymat, cbind( c(treated, control), rep(index, sum(match$subclass == j, na.rm = TRUE))))
  Treat = c(Treat, c(rep(1,length(treated)), rep(0, length(control)) ) )
}

sampleSize.mutual1.nRep = nrow(Ymat)
sc.mutual1 = hodgeslehmann(y= Ymat[,1], z= Treat,
                           st= Ymat[,2], align="hl")
res.mutual1 = senstrat(sc=sc.mutual1, z= Treat,
                       st= Ymat[,2], gamma=1.0)
res.mutual1.10 = senstrat(sc=sc.mutual1, z= Treat,
                          st= Ymat[,2], gamma=1.10)
res.mutual1.20 = senstrat(sc=sc.mutual1, z= Treat,
                       st= Ymat[,2], gamma=1.20)
res.mutual1.30 = senstrat(sc=sc.mutual1, z= Treat,
                       st= Ymat[,2], gamma=1.30)

## with IV 2
match = tryCatch(matchit(Z2 ~ Z1 + Z3 + age + sickle + village, data = subdat, method = "exact"), 
                 error=function(err) NA); 
Ymat = Treat = c()
nj = length(table(match$subclass))
index = 0
for(j in 1:nj){
  index = index + 1
  treated = subdat$Y[which(match$subclass == j & subdat$Z2 == 1)]
  control = subdat$Y[which(match$subclass == j & subdat$Z2 == 0)]
  Ymat = rbind(Ymat, cbind( c(treated, control), rep(index, sum(match$subclass == j, na.rm = TRUE))))
  Treat = c(Treat, c(rep(1,length(treated)), rep(0, length(control)) ) )
}

sampleSize.mutual2.nRep = nrow(Ymat)
sc.mutual2 = hodgeslehmann(y= Ymat[,1], z= Treat,
                           st= Ymat[,2], align="hl")
res.mutual2 = senstrat(sc=sc.mutual2, z= Treat,
                       st= Ymat[,2], gamma=1.0)
res.mutual2.10 = senstrat(sc=sc.mutual2, z= Treat,
                          st= Ymat[,2], gamma=1.10)
res.mutual2.20 = senstrat(sc=sc.mutual2, z= Treat,
                       st= Ymat[,2], gamma=1.20)
res.mutual2.30 = senstrat(sc=sc.mutual2, z= Treat,
                          st= Ymat[,2], gamma=1.30)

## with IV 3
match = tryCatch(matchit(Z3 ~ Z2 + age + sickle + village, data = subdat, method = "exact"), 
                 error=function(err) NA); 
Ymat = Treat = c()
nj = length(table(match$subclass))
index = 0
for(j in 1:nj){
  index = index + 1
  treated = subdat$Y[which(match$subclass == j & subdat$Z3 == 1)]
  control = subdat$Y[which(match$subclass == j & subdat$Z3 == 0)]
  Ymat = rbind(Ymat, cbind( c(treated, control), rep(index, sum(match$subclass == j, na.rm = TRUE))))
  Treat = c(Treat, c(rep(1,length(treated)), rep(0, length(control)) ) )
}

sampleSize.mutual3.nRep = nrow(Ymat)
sc.mutual3 = hodgeslehmann(y= Ymat[,1], z= Treat,
                           st= Ymat[,2], align="hl")
res.mutual3 = senstrat(sc=sc.mutual3, z= Treat,
                       st= Ymat[,2], gamma=1.0)
res.mutual3.10 = senstrat(sc=sc.mutual3, z= Treat,
                       st= Ymat[,2], gamma=1.10)
res.mutual3.20 = senstrat(sc=sc.mutual3, z= Treat,
                          st= Ymat[,2], gamma=1.20)
res.mutual3.30 = senstrat(sc=sc.mutual3, z= Treat,
                          st= Ymat[,2], gamma=1.30)
```                          
