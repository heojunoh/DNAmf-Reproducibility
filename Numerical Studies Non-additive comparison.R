### Non-Additive Function ###
fl <- function(x, t){
  term1 <- sin(10 * pi * x / (5+t))
  term2 <- 0.2 * sin(8 * pi * x)
  term1 + term2
}

### training data ###
n1 <- 13; n2 <- 10; n3 <- 7; n4 <- 4; n5 <- 1;
m1 <- 2.5; m2 <- 2.0; m3 <- 1.5; m4 <- 1.0; m5 <- 0.5;
d <- 1
eps <- sqrt(.Machine$double.eps)
x <- seq(0,1,0.01)

rep <- 100
result.tuononlinear.rmse <- matrix(NA, rep, 5)
result.tuononlinear.crps <- matrix(NA, rep, 5)
colnames(result.tuononlinear.rmse) <- c("DNAmf", "RNAmf", "MuFiCokriging", "Fractional BM", "BM")
colnames(result.tuononlinear.crps) <- c("DNAmf", "RNAmf", "MuFiCokriging", "Fractional BM", "BM")
beta <- delta <- c()

for(i in 1:rep) {
  set.seed(i)
  print(i)
  
  NestDesign <- NestedX(c(n1,n2,n3,n4,n5),d)
  
  X1 <- NestDesign[[1]]
  X2 <- NestDesign[[2]]
  X3 <- NestDesign[[3]]
  X4 <- NestDesign[[4]]
  X5 <- NestDesign[[5]]
  
  y1 <- fl(X1, t=m1)
  y2 <- fl(X2, t=m2)
  y3 <- fl(X3, t=m3)
  y4 <- fl(X4, t=m4)
  y5 <- fl(X5, t=m5)
  
  
  ### DNAmf-sqex ###
  fit.DNAmf.sqex <- DNAmf.sqex(X1, y1, X=rbind(X2, X3, X4, X5), y=rbind(y2, y3, y4, y5),
                               nn=c(length(y1),length(y2),length(y3),length(y4),length(y5)),
                               t=c(m1,m2,m3,m4,m5), multi.start=10, constant=TRUE)
  pred.DNAmf.sqex <- predict.DNAmf(fit.DNAmf.sqex, x, targett=0)
  predydiffu.sqex <- pred.DNAmf.sqex$mu
  predsig2diffu.sqex <- pred.DNAmf.sqex$sig2
  beta <- c(beta, fit.DNAmf.sqex$fit2$beta)
  delta <- c(delta, fit.DNAmf.sqex$fit2$delta)
  
  ### RNAmf ###
  fit.RNAmf <- RNAmf(list(X2, X3, X4), list(y2, y3, y4), kernel="sqex", constant=TRUE)
  pred.RNAmf <- predict(fit.RNAmf, x)
  predy.RNAmf <- pred.RNAmf$mu[[3]]
  predsig2.RNAmf <- pred.RNAmf$sig2[[3]]
  
  ### Fractional Brownian motion ###
  fit.fbm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5), 
                        rbind(matrix(rep(m1, nrow(X1))),matrix(rep(m2, nrow(X2))),matrix(rep(m3, nrow(X3))),
                              matrix(rep(m4, nrow(X4))),matrix(rep(m5, nrow(X5))) ), 
                        rbind(y1, y2, y3, y4, y5) )
  pred.fbm <- predict.MuFiMeshGP(fit.fbm, matrix(x), rep(0,length(x)))
  predyfbm <- pred.fbm$mean
  predsig2fbm <- pred.fbm$sd
  
  ### Brownian motion - tuononlinear Wu Yu ###
  fit.bm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5), 
                       rbind(matrix(rep(m1, nrow(X1))),matrix(rep(m2, nrow(X2))),matrix(rep(m3, nrow(X3))),
                             matrix(rep(m4, nrow(X4))),matrix(rep(m5, nrow(X5))) ), 
                       rbind(y1, y2, y3, y4, y5), 
                       H.known = 0.5)
  pred.bm <- predict.MuFiMeshGP(fit.bm, matrix(x), rep(0,length(x)))
  predybm <- pred.bm$mean
  predsig2bm <- pred.bm$sd
  
  ### Cokriging ###
  NestDesignMuFiCokriging <- NestedDesignBuild(design = list(X2,X3,X4))
  X2 <- NestDesignMuFiCokriging$PX
  X3 <- ExtractNestDesign(NestDesignMuFiCokriging,2)
  X4 <- ExtractNestDesign(NestDesignMuFiCokriging,3)
  y2 <- fl(X2, t=m2)
  y3 <- fl(X3, t=m3)
  y4 <- fl(X4, t=m4)
  
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesignMuFiCokriging, covtype="gauss",
                           lower=eps, upper=0.4,
                           # coef.trend = list(0,c(0,0),c(0,0)),
                           response = list(y2,y3,y4), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  
  ### RMSE ###
  result.tuononlinear.rmse[i,1] <- sqrt(mean((predydiffu.sqex-fl(x, t=0))^2)) #
  result.tuononlinear.rmse[i,2] <- sqrt(mean((predy.RNAmf-fl(x, t=0))^2)) #
  result.tuononlinear.rmse[i,3] <- sqrt(mean((pred.muficokm$mean-fl(x, t=0))^2)) #
  result.tuononlinear.rmse[i,4] <- sqrt(mean((predyfbm-fl(x, t=0))^2)) #
  result.tuononlinear.rmse[i,5] <- sqrt(mean((predybm-fl(x, t=0))^2)) #
  
  result.tuononlinear.crps[i,1] <- mean(crps(fl(x, t=0), predydiffu.sqex, predsig2diffu.sqex)) #
  result.tuononlinear.crps[i,2] <- mean(crps(fl(x, t=0), predy.RNAmf, predsig2.RNAmf)) #
  result.tuononlinear.crps[i,3] <- mean(crps(fl(x, t=0), pred.muficokm$mean, pred.muficokm$sig2)) #
  result.tuononlinear.crps[i,4] <- mean(crps(fl(x, t=0), predyfbm, predsig2fbm)) #
  result.tuononlinear.crps[i,5] <- mean(crps(fl(x, t=0), predybm, predsig2bm)) #
  
  print(result.tuononlinear.rmse[1:i,, drop=FALSE])
  print(apply(result.tuononlinear.rmse[1:i,, drop=FALSE], 2, mean))
  print(apply(result.tuononlinear.crps[1:i,, drop=FALSE], 2, mean))
}

apply(result.tuononlinear.rmse, 2, mean) 
apply(result.tuononlinear.crps, 2, mean) 
print(mean(beta))
print(mean(delta))

