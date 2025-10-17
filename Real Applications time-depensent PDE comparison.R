### Time dependent PDE; Heat equation; smaller t corresponds to longer time point (0.5-t)  ###
heat <- function(X, t, L = 1.0, N = 1000, D = 0.1) {
  if (!is.matrix(X)) {
    if (!is.numeric(X) || any(X < 0) || any(X > L))
      stop("'X'must be numeric (vector or matrix) in [0, ", L, "]")
    Xmat <- matrix(X, ncol = 1)
  } else {
    if (!is.numeric(X) || any(X < 0) || any(X > L))
      stop("All entries of 'X' must lie in [0, ", L, "]")
    Xmat <- X
  }
  
  if (length(t) != 1 || !is.numeric(t) || t < 0)
    stop("'t' must be a single numeric >= 0")
  
  ### grid & initial condition ###
  grid <- setup.grid.1D(x.up = 0, x.down = L, N = N)
  u_initial <- exp(-100 * (grid$x.mid - L/2)^2)
  
  ### define model ###
  heat_model <- function(time, u, parms) {
    flux <- tran.1D(C = u, D = D, dx = grid)
    list(dC = flux$dC)
  }
  
  ### solve equation at times 0 and t ###
  sol <- ode.1D(y = u_initial, times = c(0, 0.5-t), func = heat_model,
                parms = NULL, nspec = 1, dimens = N)
  u_t <- sol[2, -1] # row 2; time = t
  
  ### interpolation (clamp at endpoints) ###
  x_vec <- as.vector(Xmat)
  u_vec <- approx(x = grid$x.mid, y = u_t, xout = x_vec, rule = 2)$y
  
  ### return ###
  matrix(u_vec, nrow = nrow(Xmat), ncol = ncol(Xmat), dimnames = dimnames(Xmat))
}

### training data ###
c <- 0.7; gam <- -1.5
n5 <- 3; n4 <- round(c^gam*n5); n3 <- round(c^(2*gam)*n5); n2 <- round(c^(3*gam)*n5); n1 <- round(c^(4*gam)*n5);

m1 <- 0.5; m2 <- m1*c; m3 <- m2*c; m4 <- m3*c; m5 <- m4*c; 
d <- 1; targett=0
eps <- sqrt(.Machine$double.eps)

### test data ###
x.test <- matrix(seq(0,1,0.01))
y.test <- heat(x.test,0)

rep <- 100
result.heat.rmse <- result.heat.crps <- matrix(NA, rep, 6)
colnames(result.heat.rmse) <- colnames(result.heat.crps) <- c("DNAmf", "RNAmf", "Cokriging", "NARGP", "Fractional BM", "BM")
beta <- delta <- c()

par(mfrow=c(1,2))
for(i in 1:rep) {
  set.seed(i)
  print(i)
  
  NestDesign <- NestedX(c(n1,n2,n3,n4,n5),d)
  
  X1 <- NestDesign[[1]]
  X2 <- NestDesign[[2]]
  X3 <- NestDesign[[3]]
  X4 <- NestDesign[[4]]
  X5 <- NestDesign[[5]]
  y1 <- heat(X1,m1)
  y2 <- heat(X2,m2)
  y3 <- heat(X3,m3)
  y4 <- heat(X4,m4)
  y5 <- heat(X5,m5)
  
  if(NARGP) saveRDS(list(X1=X3, X2=X4, X3=X5, Y1=y3, Y2=y4, Y3=y5, Xtest=x.test, Ytest=y.test), file = "tmp_data.rds")
  
  ### DNAmf-sqex ###
  tic.DNAmf <- proc.time()[3]
  fit.DNAmf <- DNAmf.sqex(X1, y1,
                               X=rbind(X2, X3, X4, X5),
                               y=rbind(y2, y3, y4, y5),
                               nn=c(length(y1),length(y2),length(y3),length(y4),length(y5)),
                               t=c(m1,m2,m3,m4,m5), multi.start=10, constant=TRUE)
  pred.DNAmf <- predict.DNAmf(fit.DNAmf, x.test, targett=targett)
  predydiffu <- pred.DNAmf$mu
  predsig2diffu <- pred.DNAmf$sig2
  beta <- c(beta, fit.DNAmf$fit2$beta)
  delta <- c(delta, fit.DNAmf$fit2$delta)
  toc.DNAmf <- proc.time()[3]
  
  ### RNAmf ###
  fit.RNAmf <- RNAmf(list(X3, X4, X5), list(y3, y4, y5), kernel="sqex", constant=TRUE)
  pred.RNAmf <- predict(fit.RNAmf, x.test)
  predy.RNAmf <- pred.RNAmf$mu[[3]]
  predsig2.RNAmf <- pred.RNAmf$sig2[[3]]
  
  ### Fractional Brownian motion ###
  tic.fbm <- proc.time()[3]
  fit.fbm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5), 
                        rbind(matrix(rep(m1, nrow(X1))),matrix(rep(m2, nrow(X2))),matrix(rep(m3, nrow(X3))),
                              matrix(rep(m4, nrow(X4))),matrix(rep(m5, nrow(X5))) ), 
                        rbind(y1, y2, y3, y4, y5) )
  pred.fbm <- predict(fit.fbm, matrix(x.test), rep(0,length(x.test)))
  predyfbm <- pred.fbm$mean
  predsig2fbm <- pred.fbm$sd
  toc.fbm <- proc.time()[3]
  
  ### Brownian motion - heat Wu Yu ###
  tic.bm <- proc.time()[3]
  fit.bm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5), 
                       rbind(matrix(rep(m1, nrow(X1))),matrix(rep(m2, nrow(X2))),matrix(rep(m3, nrow(X3))),
                             matrix(rep(m4, nrow(X4))),matrix(rep(m5, nrow(X5))) ), 
                       rbind(y1, y2, y3, y4, y5), 
                       H.known = 0.5)
  pred.bm <- predict(fit.bm, matrix(x.test), rep(0,length(x.test)))
  predybm <- pred.bm$mean
  predsig2bm <- pred.bm$sd
  toc.bm <- proc.time()[3]
  
  ### Cokriging ###
  NestDesignMuFiCokriging <- NestedDesignBuild(design = list(X3,X4,X5))
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesignMuFiCokriging, covtype="gauss",
                           lower=eps, upper=0.3,
                           # coef.trend = list(0,c(0,0),c(0,0)),
                           response = list(y3,y4,y5), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, x.test, "SK")
  
  ### NARGP ###
  if(NARGP) py <- py_run_file("NARGP.py")
  
  
  ### RMSE ###
  result.heat.rmse[i,1] <- sqrt(mean((predydiffu-y.test)^2)) #
  result.heat.rmse[i,2] <- sqrt(mean((predy.RNAmf-y.test)^2)) #
  result.heat.rmse[i,3] <- sqrt(mean((pred.muficokm$mean-y.test)^2)) #
  if(NARGP) result.heat.rmse[i,4] <- py$error
  result.heat.rmse[i,5] <- sqrt(mean((predyfbm-y.test)^2)) #
  result.heat.rmse[i,6] <- sqrt(mean((predybm-y.test)^2)) #
  
  result.heat.crps[i,1] <- mean(crps(y.test, predydiffu, predsig2diffu)) #
  result.heat.crps[i,2] <- mean(crps(y.test, predy.RNAmf, predsig2.RNAmf)) #
  result.heat.crps[i,3] <- mean(crps(y.test, pred.muficokm$mean, pred.muficokm$sig2)) #
  if(NARGP) result.heat.crps[i,4] <- py$crps
  result.heat.crps[i,5] <- mean(crps(y.test, predyfbm, predsig2fbm)) #
  result.heat.crps[i,6] <- mean(crps(y.test, predybm, predsig2bm)) #
  
  # boxplot(result.heat.rmse[1:i,,drop=FALSE])
  # boxplot(result.heat.crps[1:i,,drop=FALSE])
  print(result.heat.rmse[1:i,,drop=FALSE])
  print(result.heat.crps[1:i,,drop=FALSE])
}

print(mean(beta))
print(mean(delta))

