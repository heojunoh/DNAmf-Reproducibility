### Borehole Example ###
borehole <- function(xx, tt)
{
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- (2*pi - tt^2) * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+tt^2+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

output.f <- function(x,t){
  factor_range <- list("rw" = c(0.05, 0.15), "r" = c(100, 50000),
                       "Tu" = c(63070, 115600), "Hu" = c(990, 1110),
                       "Tl" = c(63.1, 116), "Hl" = c(700, 820),
                       "L" = c(1120, 1680), "Kw" = c(9855, 12045))
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  borehole(x[1:8], tt=t)
} 

### training data ###
n6 <- 10; n5 <- 5+n6;  n4 <- 5+n5; n3 <- 5+n4;  n2 <- 5+n3; n1 <- 5+n2; 
m1 <- 3; m2 <- 2.5; m3 <- 2; m4 <- 1.5; m5 <- 1; m6 <- 0.5;
d <- 8
eps <- sqrt(.Machine$double.eps)
set.seed(0)
x <- maximinLHS(100, d)

rep <- 100
result.borehole.rmse <- matrix(NA, rep, 5)
result.borehole.crps <- matrix(NA, rep, 5)
colnames(result.borehole.rmse) <- c("DNAmf", "RNAmf", "MuFiCokriging", "Fractional BM", "BM")
colnames(result.borehole.crps) <- c("DNAmf", "RNAmf", "MuFiCokriging", "Fractional BM", "BM")
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
  
  y1 <- matrix(apply(X1, 1, output.f, t=m1))
  y2 <- matrix(apply(X2, 1, output.f, t=m2))
  y3 <- matrix(apply(X3, 1, output.f, t=m3))
  y4 <- matrix(apply(X4, 1, output.f, t=m4))
  y5 <- matrix(apply(X5, 1, output.f, t=m5))
  y6 <- matrix(apply(X6, 1, output.f, t=m6))
  
  
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
  y4 <- matrix(apply(X4, 1, output.f, t=m4))
  y5 <- matrix(apply(X5, 1, output.f, t=m5))
  y6 <- matrix(apply(X6, 1, output.f, t=m6))
  
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesignMuFiCokriging, covtype="gauss",
                           lower=eps, upper=0.4,
                           # coef.trend = list(0,c(0,0),c(0,0)),
                           response = list(y4,y5,y6), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  
  ### RMSE ###
  result.borehole.rmse[i,1] <- sqrt(mean((predydiffu.sqex-matrix(apply(x, 1, output.f, t=0)))^2)) #
  result.borehole.rmse[i,2] <- sqrt(mean((predy.RNAmf-matrix(apply(x, 1, output.f, t=0)))^2)) #
  result.borehole.rmse[i,3] <- sqrt(mean((pred.muficokm$mean-matrix(apply(x, 1, output.f, t=0)))^2)) #
  result.borehole.rmse[i,4] <- sqrt(mean((predyfbm-matrix(apply(x, 1, output.f, t=0)))^2)) #
  result.borehole.rmse[i,5] <- sqrt(mean((predybm-matrix(apply(x, 1, output.f, t=0)))^2)) #
  
  result.borehole.crps[i,1] <- mean(crps(matrix(apply(x, 1, output.f, t=0)), predydiffu.sqex, predsig2diffu.sqex)) #
  result.borehole.crps[i,2] <- mean(crps(matrix(apply(x, 1, output.f, t=0)), predy.RNAmf, predsig2.RNAmf)) #
  result.borehole.crps[i,3] <- mean(crps(matrix(apply(x, 1, output.f, t=0)), pred.muficokm$mean, pred.muficokm$sig2)) #
  result.borehole.crps[i,4] <- mean(crps(matrix(apply(x, 1, output.f, t=0)), predyfbm, predsig2fbm)) #
  result.borehole.crps[i,5] <- mean(crps(matrix(apply(x, 1, output.f, t=0)), predybm, predsig2bm)) #
  
  print(result.borehole.rmse[1:i,, drop=FALSE])
  print(result.borehole.crps[1:i,, drop=FALSE])
  print(apply(result.borehole.rmse[1:i,, drop=FALSE], 2, mean))
  print(apply(result.borehole.crps[1:i,, drop=FALSE], 2, mean))
}

apply(result.borehole.rmse, 2, mean) 
apply(result.borehole.crps, 2, mean) 
print(mean(beta))
print(mean(delta))
