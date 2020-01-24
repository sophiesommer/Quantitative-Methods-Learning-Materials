require(MASS)
require(eulerr)
require(VennDiagram)
simdatafunc = function(sampsize=100, beta1=1, beta2=1, sigma2=1, corr=-.5){
  data = mvrnorm(n=sampsize, mu=c(0, 0),
                 Sigma=matrix(c(1, corr, corr, 1), nrow=2),
                 empirical=TRUE)
  X1 = data[, 1]
  X2 = data[, 2]
  Y = 1 + beta1 * X1 + beta2 * X2 + rnorm(sampsize, 0, sigma2^.5)
  return(data.frame(X1 = X1,X2 = X2,Y = Y))
}

calcSS = function(lmobj){
  r2 <- summary(lmobj)$r.squared
  res <- lmobj$residuals
  RSS <- sum(res^2)
  TSS <- (RSS)/(-(r2-1))
  MSS <- TSS - RSS
  return(data.frame(MSS=MSS, TSS=TSS, RSS=RSS))
}

simdata = simdatafunc()
lmx1x2 = lm(Y~X1+X2, data=simdata)
lmx1 = lm(Y~X1, data=simdata)
lmx2 = lm(Y~X2, data=simdata)
MSSx1 = calcSS(lmx1)$MSS
MSSx2 = calcSS(lmx2)$MSS
MSSx1x2 = calcSS(lmx1x2)$MSS
TSS = calcSS(lmx1x2)$TSS
RSS = calcSS(lmx1x2)$RSS
middle = MSSx1 + MSSx2 - MSSx1x2

fit <- euler(c('MSS X1' = 0, 'MSS X2' = 0, 'TSS' =round(TSS-(MSSx1+MSSx2-middle)), 
               "MSS X1&MSS X2" = 0, 
               "MSS X2&TSS" = round(MSSx2-middle),
               "MSS X1&TSS" = round(MSSx1-middle),
               "MSS X1&MSS X2&TSS" = round(middle)),
                names=c(1:7))

plot(fit, quantities=TRUE, legend=TRUE)


fit <- euler(c('MSS X1' = 0, 'TSS' =round(TSS-MSSx1), 
               "MSS X1&TSS" = round(MSSx1)))
plot(fit, quantities=TRUE, legend=TRUE)







fit <- euler(c(A = 0, B = 0, C=70, 
               "A&B" = 70, 
               "B&C" = 70,
               "A&C" = 70,
               "A&B&C" = 90))
plot(fit)




venn.plot <- draw.triple.venn(
  area1 = round(MSSx1),
  area2 = round(MSSx2),
  area3 = round(TSS),
  n12 = round(middle),
  n23 = round(MSSx2),
  n13 = round(MSSx1),
  n123 = round(middle),
  category = c("C1", "C2", "C3"),
  fill = c("blue", "red", "green"),
  scaled=TRUE)



su <- venneuler(c(A=162, B=104, C=86, "A&B"=206, "A&C"=112, "B&C"=90 ,"A&B&C"=2433))
plot(su)


draw.triple.venn(area1 = round(MSSx1), 
                 area2 = round(MSSx2), 
                 area3 = round(TSS), 
                 n12 = round(middle), 
                 n23 = round(MSSx2), 
                 n13 = round(MSSx1), 
                 n123 = round(MSSx1+MSSx2-middle), 
                 category = c(paste("MSS X1 =",round(MSSx1)),
                              paste("MSS X2 =",round(MSSx2)),
                              paste("MSS X1 & X2 =",round(MSSx1x2))), 
                 lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))


VennDiag <- euler(c("A" = round(MSSx1), "B" = round(MSSx2), "C" = round(TSS), 
                    "A&B" = round(middle), "B&C" = round(MSSx2), 
                    "A&C" = round(MSSx1), "A&B&C" = round(MSSx1+MSSx2-middle)))
plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("grey", "lightgrey", "darkgrey"))


slices <- c((MSSx1 - middle), RSS, (MSSx2 - middle), middle)
lbls <- c("MSS X1","RSS", "MSS X2", "MSS X1 & X2")
pie(slices, labels = lbls, main="TSS Breakdown", col=c("blue", "red", "white", "light blue"))

