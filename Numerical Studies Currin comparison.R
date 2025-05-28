### Currin Function ###
currin <- function(xx, t){
  x1 <- xx[1]
  x2 <- xx[2]
  fact1 <- exp(-4*t) - exp(-1/(2*x2))
  fact2 <- 2300*x1^3 + 1900*x1^2 + 2092*x1 + 60
  fact3 <- 100*x1^3 + 500*x1^2 + 4*x1 + 20
  
  fact1 * fact2/fact3
}

### training data ###
n6 <- 4; n5 <- 3+n6; n4 <- 3+n5; n3 <- 3+n4; n2 <- 3+n3; n1 <- 3+n2; 
m1 <- 0.6; m2 <- 0.5; m3 <- 0.4; m4 <- 0.3; m5 <- 0.2; m6 <- 0.1
d <- 2
eps <- sqrt(.Machine$double.eps)
set.seed(0)
x <- maximinLHS(100, d)

rep <- 100
result.currin.rmse <- matrix(NA, rep, 5)
result.currin.crps <- matrix(NA, rep, 5)
colnames(result.currin.rmse) <- c("DNAmf", "RNAmf", "MuFiCokriging", "Fractional BM", "BM")
colnames(result.currin.crps) <- c("DNAmf", "RNAmf", "MuFiCokriging", "Fractional BM", "BM")
beta <- delta <- c()

for(i in 1:rep) {
  set.seed(i)
  print(i)
  
  NestDesign <- NestedX(c(n1,n2,n3,n4,n5,n6),d)
  
  X1 <- NestDesign[[1]]
  X2 <- NestDesign[[2]]
  X3 <- NestDesign[[3]]
  X4 <- NestDesign[[4]]
  X5 <- NestDesign[[5]]
  X6 <- NestDesign[[6]]
  
  y1 <- matrix(apply(X1, 1, currin, t=m1))
  y2 <- matrix(apply(X2, 1, currin, t=m2))
  y3 <- matrix(apply(X3, 1, currin, t=m3))
  y4 <- matrix(apply(X4, 1, currin, t=m4))
  y5 <- matrix(apply(X5, 1, currin, t=m5))
  y6 <- matrix(apply(X6, 1, currin, t=m6))
  
  
  ### DNAmf-sqex ###
  fit.DNAmf.sqex <- DNAmf.sqex(X1, y1, X=rbind(X2, X3, X4, X5, X6), y=rbind(y2, y3, y4, y5, y6),
                               nn=c(length(y1),length(y2),length(y3),length(y4),length(y5),length(y6)),
                               t=c(m1,m2,m3,m4,m5,m6), multi.start=20, constant=TRUE)
  pred.DNAmf.sqex <- predict.DNAmf(fit.DNAmf.sqex, x, targett=0)
  predydiffu.sqex <- pred.DNAmf.sqex$mu
  predsig2diffu.sqex <- pred.DNAmf.sqex$sig2
  beta <- c(beta, fit.DNAmf.sqex$fit2$beta)
  delta <- c(delta, fit.DNAmf.sqex$fit2$delta)
  
  ### RNAmf ###
  fit.RNAmf <- RNAmf_three_level(X4, y4, X5, y5, X6, y6, kernel="sqex", constant=TRUE)
  pred.RNAmf <- predict(fit.RNAmf, x)
  predy.RNAmf <- pred.RNAmf$mu
  predsig2.RNAmf <- pred.RNAmf$sig2
  
  ### Fractional Brownian motion ###
  fit.fbm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5, X6), 
                               rbind(matrix(rep(m1, nrow(X1))),matrix(rep(m2, nrow(X2))),matrix(rep(m3, nrow(X3))),
                                     matrix(rep(m4, nrow(X4))),matrix(rep(m5, nrow(X5))),matrix(rep(m6, nrow(X6))) ), 
                               rbind(y1, y2, y3, y4, y5, y6) )
  pred.fbm <- predict.MuFiMeshGP(fit.fbm, x, matrix(rep(0,nrow(x))) )
  predyfbm <- pred.fbm$mean
  predsig2fbm <- pred.fbm$sd
  
  ### Brownian motion - Tuo Wu Yu ###
  fit.bm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5, X6), 
                               rbind(matrix(rep(m1, nrow(X1))),matrix(rep(m2, nrow(X2))),matrix(rep(m3, nrow(X3))),
                                     matrix(rep(m4, nrow(X4))),matrix(rep(m5, nrow(X5))),matrix(rep(m6, nrow(X6))) ), 
                               rbind(y1, y2, y3, y4, y5, y6), 
                               H.known = 0.5)
  pred.bm <- predict.MuFiMeshGP(fit.bm, x, matrix(rep(0,nrow(x))) )
  predybm <- pred.bm$mean
  predsig2bm <- pred.bm$sd
  
  ### Cokriging ###
  NestDesignMuFiCokriging <- NestedDesignBuild(design = list(X4,X5,X6))
  X4 <- NestDesignMuFiCokriging$PX
  X5 <- ExtractNestDesign(NestDesignMuFiCokriging,2)
  X6 <- ExtractNestDesign(NestDesignMuFiCokriging,3)
  y4 <- matrix(apply(X4, 1, currin, t=m4))
  y5 <- matrix(apply(X5, 1, currin, t=m5))
  y6 <- matrix(apply(X6, 1, currin, t=m6))
  
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesignMuFiCokriging, covtype="gauss",
                           lower=eps, upper=0.4,
                           # coef.trend = list(0,c(0,0),c(0,0)),
                           response = list(y4,y5,y6), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  
  ### RMSE ###
  result.currin.rmse[i,1] <- sqrt(mean((predydiffu.sqex-matrix(apply(x, 1, currin, t=0)))^2)) #
  result.currin.rmse[i,2] <- sqrt(mean((predy.RNAmf-matrix(apply(x, 1, currin, t=0)))^2)) #
  result.currin.rmse[i,3] <- sqrt(mean((pred.muficokm$mean-matrix(apply(x, 1, currin, t=0)))^2)) #
  result.currin.rmse[i,4] <- sqrt(mean((predyfbm-matrix(apply(x, 1, currin, t=0)))^2)) #
  result.currin.rmse[i,5] <- sqrt(mean((predybm-matrix(apply(x, 1, currin, t=0)))^2)) #
  
  result.currin.crps[i,1] <- mean(crps(matrix(apply(x, 1, currin, t=0)), predydiffu.sqex, predsig2diffu.sqex)) #
  result.currin.crps[i,2] <- mean(crps(matrix(apply(x, 1, currin, t=0)), predy.RNAmf, predsig2.RNAmf)) #
  result.currin.crps[i,3] <- mean(crps(matrix(apply(x, 1, currin, t=0)), pred.muficokm$mean, pred.muficokm$sig2)) #
  result.currin.crps[i,4] <- mean(crps(matrix(apply(x, 1, currin, t=0)), predyfbm, predsig2fbm)) #
  result.currin.crps[i,5] <- mean(crps(matrix(apply(x, 1, currin, t=0)), predybm, predsig2bm)) #
  
  print(result.currin.rmse[1:i,, drop=FALSE])
  print(result.currin.crps[1:i,, drop=FALSE])
  print(apply(result.currin.rmse[1:i,, drop=FALSE], 2, mean))
  print(apply(result.currin.crps[1:i,, drop=FALSE], 2, mean))
}

apply(result.currin.rmse, 2, mean) 
apply(result.currin.crps, 2, mean) 
print(mean(beta))
print(mean(delta))
