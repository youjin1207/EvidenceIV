library(MatchIt)
library(senstrat)
library(xtable)

# dat: original data
attach(dat)
table(dat$HBASCH)
year=substr(as.character(dat$DATEVISIT),4,7)
year=as.numeric(year)
# Was the village cluster randomized to receive bednets
villcluster.bednet=dat$villnewflat%in%c(1,2,3,5,6,92,93,96,98)
# Mean HAZ by village assigned to receive bed nets/time period (pre-intervention, post intervention)
mean(dat$HAZ[villcluster.bednet==1 & year<=1996],na.rm=TRUE)
mean(dat$HAZ[villcluster.bednet==1 & year>=1997],na.rm=TRUE)
mean(dat$HAZ[villcluster.bednet==0 & year<=1996],na.rm=TRUE)
mean(dat$HAZ[villcluster.bednet==0 & year>=1997],na.rm=TRUE)

### Allow for repeated child observation
dat$birthmonth=as.numeric(substr(as.character(dat$newdob_main),1,2))
dat$birthyear=as.numeric(substr(as.character(dat$newdob_main),7,11))
dat$age.in.months=(dat$YEAR-dat$birthyear)*12+(dat$MONTH-dat$birthmonth)
dat$villcluster.bednet=dat$villnewflat%in%c(1,2,3,5,6,92,93,96,98)
dat$time.under.bednet.randomization=dat$villcluster.bednet*(dat$YEAR>=1997)*pmin(((dat$YEAR-1997)*12+dat$MONTH)/dat$age.in.months,1)
dat$symptomatic.malaria=(dat$TEMP>=37.5)*(dat$posneg==1)
dat$num.symptomatic.malaria.episodes=rep(0,nrow(dat))
for(i in 1:nrow(dat)){
  dat$num.symptomatic.malaria.episodes[i]=sum((dat$OBNCODE[1:i]==dat$OBNCODE[i])*dat$symptomatic.malaria[1:i],na.rm=TRUE)
}

dat$symptomatic.malaria.per.age=dat$num.symptomatic.malaria.episodes/dat$age.in.months
dat$sickle.cell.trait=rep(NA,nrow(dat))
dat$sickle.cell.trait[dat$HBASCH==1]="HbAS"
dat$sickle.cell.trait[dat$HBASCH==2]="HbAA"

# outcome Y : height for age Z-score
# exposure D : the child's clinical malaria incidence rate (number of clinical malaria episodes divided by the child's age)
# instrumental variable Z : time.under.bednet.randomization + as.factor(sickle.cell.trait)
# covairates: age (four categories) + sickle cell + as.factor(villnewflat)

dat$time.under.bednet.randomization=dat$villcluster.bednet*(dat$YEAR>=1997)*pmin(((dat$YEAR-1997)*12+dat$MONTH)/dat$age.in.months,1)
summary((dat$YEAR>=1997)*((dat$YEAR-1997)*12+dat$MONTH)/dat$age.in.months)

dat$Z1 = ifelse(dat$villcluster.bednet*(dat$YEAR>=1997) == 1, 1, 0)
dat$Z2 = ifelse(dat$time.under.bednet.randomization <= 0.2 | dat$villcluster.bednet*(dat$YEAR>=1997) == 0, 0, 1)
dat$Z3 = ifelse(dat$time.under.bednet.randomization <= 0.5 | dat$villcluster.bednet*(dat$YEAR>=1997) == 0, 0, 1)

subdata = data.frame(Z1 = dat$Z1, Z2 = dat$Z2, Z3 = dat$Z3, Y = dat$HAZ, D = dat$symptomatic.malaria.per.age,
                     age.in.months = dat$age.in.months, YEAR = dat$YEAR, 
                     sickle = as.factor(dat$sickle.cell.trait), time.under.bednet.randomization = dat$time.under.bednet.randomization,
                     village = dat$villnewflat)
subdata = na.omit(subdata)
subdata = subdata[subdata$age.in.months <= 12*5,] # 0-5 years old
subdata = subdata[subdata$YEAR >= 1995,] # outcomes measured before or after 1995

subdata$age = ifelse(subdata$age <= 6, 1, ifelse(subdata$age <= 14, 2, 
                                                  ifelse(subdata$age <= 27, 3, 4))) # categorize age by quartiles

########## sample data ########
dat = read.csv("Data/sample.csv", header = TRUE, sep = ",")

## associations betwen each instrument and the exposure variable
fitD = lm(D ~ Z1 + Z2 + Z3 + age + sickle + as.factor(village), data = subdata)

## exploratory figure ###
summary(subdata$time.under.bednet.randomization)
pdf("Figure/malaria.pdf",  width = 9, height = 6)
par(mfrow = c(1,1),   mar = c(5,5,5,3),  cex.lab = 1.5, 
    cex.main = 1.5, cex.axis = 1.5, tcl = 0.5, lwd = 2)
plot(subdata$time.under.bednet.randomization, subdata$D, cex = 0.2, ylim = c(0, 0.7),
     ylab = "Symptomatic malaria per age",
     xlab = "% of lifetime exposed to the bednet",
     main =  "Bednet randomization as an instrument of malaria")
lines(smooth.spline(subdata$time.under.bednet.randomization, subdata$D), col = "red",
      lwd = 2)
dev.off()


## implement the mutual stratification ##
## with IV 1
match = tryCatch(matchit(Z1 ~ Z2 + age + sickle + village, data = subdata, method = "exact"), 
                 error=function(err) NA); 
Ymat = Treat = c()
nj = length(table(match$subclass))
index = 0
for(j in 1:nj){
  index = index + 1
  treated = subdata$Y[which(match$subclass == j & subdata$Z1 == 1)]
  control = subdata$Y[which(match$subclass == j & subdata$Z1 == 0)]
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
match = tryCatch(matchit(Z2 ~ Z1 + Z3 + age + sickle + village, data = subdata, method = "exact"), 
                 error=function(err) NA); 
Ymat = Treat = c()
nj = length(table(match$subclass))
index = 0
for(j in 1:nj){
  index = index + 1
  treated = subdata$Y[which(match$subclass == j & subdata$Z2 == 1)]
  control = subdata$Y[which(match$subclass == j & subdata$Z2 == 0)]
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
match = tryCatch(matchit(Z3 ~ Z2 + age + sickle + village, data = subdata, method = "exact"), 
                 error=function(err) NA); 
Ymat = Treat = c()
nj = length(table(match$subclass))
index = 0
for(j in 1:nj){
  index = index + 1
  treated = subdata$Y[which(match$subclass == j & subdata$Z3 == 1)]
  control = subdata$Y[which(match$subclass == j & subdata$Z3 == 0)]
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

## table ##
tab = matrix(NA, 6, 3)
rownames(tab) = c("t_stat (p-value) psi_k",
  "The number of strata", "Gamma = 1.0", "Gamma=1.1", "Gamma=1.2", "Gamma=1.3")
colnames(tab) = c("IV1", "IV2", "IV3")
tab[1,] = c(paste0(formatC(summary(fitD)$coefficients[2,3], format = "f", digits = 2), ", (",
                  paste(formatC(summary(fitD)$coefficients[2,4], format = "f", digits = 3)), ")"),
            paste0(formatC(summary(fitD)$coefficients[3,3], format = "f", digits = 2), ", (",
                   paste(formatC(summary(fitD)$coefficients[3,4], format = "f", digits = 3)), ")"),
            paste0(formatC(summary(fitD)$coefficients[4,3], format = "f", digits = 2), ", (",
                   paste(formatC(summary(fitD)$coefficients[4,4], format = "f", digits = 3)), ")"))
tab[2,] = c(paste0(res.mutual1$Description[1], " (",  res.mutual1$Description[2], " vs ", res.mutual1$Description[3], ")"), 
            paste0(res.mutual2$Description[1], " (",  res.mutual2$Description[2], " vs ", res.mutual2$Description[3], ")"),
            paste0(res.mutual3$Description[1], " (",  res.mutual3$Description[2], " vs ", res.mutual3$Description[3], ")"))
tab[3,] = formatC(c(res.mutual1$Result[1], res.mutual2$Result[1], res.mutual3$Result[1]), format = "f", digits = 3)
tab[4,] = formatC(c(res.mutual1.10$Result[1], res.mutual2.10$Result[1], res.mutual3.10$Result[1]), format = "f", digits = 3)
tab[5,] = formatC(c(res.mutual1.20$Result[1], res.mutual2.20$Result[1], res.mutual3.20$Result[1]), format = "f", digits = 3)
tab[6,] = formatC(c(res.mutual1.30$Result[1], res.mutual2.30$Result[1], res.mutual3.30$Result[1]), format = "f", digits = 3)

print(tab)
