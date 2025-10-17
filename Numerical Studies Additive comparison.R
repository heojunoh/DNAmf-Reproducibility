### Additive Function ###
fl <- function(x, t){
  term1 <- exp(-1.4 * x) * cos(3.5 * pi * x)
  term2 <- t^2 * (sin(40 * x))/10
  term1 + term2
}

### training data ###
n5 <- 1; n4 <- round(c^gam*n5); n3 <- round(c^(2*gam)*n5); n2 <- round(c^(3*gam)*n5); n1 <- round(c^(4*gam)*n5);
m1 <- 2.5; m2 <- m1*c; m3 <- m2*c; m4 <- m3*c; m5 <- m4*c;
d <- 1
eps <- sqrt(.Machine$double.eps)
x <- seq(0,1,0.01)
y.test <- matrix(fl(x, t=0))

rep <- 100
result.tuolinear.rmse <- result.tuolinear.crps <- matrix(NA, rep, 6)
colnames(result.tuolinear.rmse) <- colnames(result.tuolinear.crps) <- c("DNAmf", "RNAmf", "Cokriging", "NARGP", "Fractional BM", "BM")
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
  
  if(NARGP) saveRDS(list(X1=X1, X2=X2, X3=X3, Y1=y1, Y2=y2, Y3=y3, Xtest=matrix(x,ncol=d), Ytest=y.test), file = "tmp_data.rds")
  
  
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
  # fit.RNAmf <- RNAmf(list(X2, X3, X4), list(y2, y3, y4), kernel="sqex", constant=TRUE)
  fit.RNAmf <- RNAmf(list(X1, X2, X3), list(y1, y2, y3), kernel="sqex", constant=TRUE)
  pred.RNAmf <- predict(fit.RNAmf, x)
  predy.RNAmf <- pred.RNAmf$mu[[3]]
  predsig2.RNAmf <- pred.RNAmf$sig2[[3]]
  
  ### Fractional Brownian motion ###
  fit.fbm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5),
                        rbind(matrix(rep(m1, nrow(X1))),matrix(rep(m2, nrow(X2))),matrix(rep(m3, nrow(X3))),
                              matrix(rep(m4, nrow(X4))),matrix(rep(m5, nrow(X5))) ),
                        rbind(y1, y2, y3, y4, y5) )
  pred.fbm <- predict(fit.fbm, matrix(x), rep(0,length(x)))
  predyfbm <- pred.fbm$mean
  predsig2fbm <- pred.fbm$sd
  
  ### Brownian motion - tuolinear Wu Yu ###
  fit.bm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5),
                       rbind(matrix(rep(m1, nrow(X1))),matrix(rep(m2, nrow(X2))),matrix(rep(m3, nrow(X3))),
                             matrix(rep(m4, nrow(X4))),matrix(rep(m5, nrow(X5))) ),
                       rbind(y1, y2, y3, y4, y5),
                       H.known = 0.5)
  pred.bm <- predict(fit.bm, matrix(x), rep(0,length(x)))
  predybm <- pred.bm$mean
  predsig2bm <- pred.bm$sd
  
  ### Cokriging ###
  NestDesignMuFiCokriging <- NestedDesignBuild(design = list(X1,X2,X3))
  X1 <- NestDesignMuFiCokriging$PX
  X2 <- ExtractNestDesign(NestDesignMuFiCokriging,2)
  X3 <- ExtractNestDesign(NestDesignMuFiCokriging,3)
  y1 <- fl(X1, t=m1)
  y2 <- fl(X2, t=m2)
  y3 <- fl(X3, t=m3)
  
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesignMuFiCokriging, covtype="gauss",
                           lower=eps, upper=0.2,
                           response = list(y1,y2,y3), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  
  ### NARGP ###
  if(NARGP) py <- py_run_file("NARGP.py")
  
  ### RMSE ###
  result.tuolinear.rmse[i,1] <- sqrt(mean((predydiffu.sqex-fl(x, t=0))^2)) #
  result.tuolinear.rmse[i,2] <- sqrt(mean((predy.RNAmf-fl(x, t=0))^2)) #
  result.tuolinear.rmse[i,3] <- sqrt(mean((pred.muficokm$mean-fl(x, t=0))^2)) #
  if(NARGP) result.tuolinear.rmse[i,4] <- py$error
  result.tuolinear.rmse[i,5] <- sqrt(mean((predyfbm-fl(x, t=0))^2)) #
  result.tuolinear.rmse[i,6] <- sqrt(mean((predybm-fl(x, t=0))^2)) #
  
  result.tuolinear.crps[i,1] <- mean(crps(fl(x, t=0), predydiffu.sqex, predsig2diffu.sqex)) #
  result.tuolinear.crps[i,2] <- mean(crps(fl(x, t=0), predy.RNAmf, predsig2.RNAmf)) #
  result.tuolinear.crps[i,3] <- mean(crps(fl(x, t=0), pred.muficokm$mean, pred.muficokm$sig2)) #
  if(NARGP) result.tuolinear.crps[i,4] <- py$crps
  result.tuolinear.crps[i,5] <- mean(crps(fl(x, t=0), predyfbm, predsig2fbm)) #
  result.tuolinear.crps[i,6] <- mean(crps(fl(x, t=0), predybm, predsig2bm)) #
  
  print(result.tuolinear.rmse[1:i,, drop=FALSE])
  print(apply(result.tuolinear.rmse[1:i,, drop=FALSE], 2, mean))
  print(apply(result.tuolinear.crps[1:i,, drop=FALSE], 2, mean))
}

apply(result.tuolinear.rmse, 2, mean)
apply(result.tuolinear.crps, 2, mean)
print(mean(beta))
print(mean(delta))

