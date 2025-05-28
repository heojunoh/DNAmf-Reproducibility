### Time dependent PDE; Heat equation; smaller t corresponds to longer time point (0.5-t)  ###
library(ReacTran)
library(deSolve)
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
n1 <- 15; n2 <- 12; n3 <- 9; n4 <- 6; n5 <- 3;
m1 <- 0.5; m2 <- 0.4; m3 <- 0.3; m4 <- 0.2; m5 <- 0.1; targett=0
d <- 1
eps <- sqrt(.Machine$double.eps)

### test data ###
x.test <- seq(0,1,0.01)
y.test <- heat(x.test,0)

rep <- 100
result.heat.rmse <- matrix(NA, rep, 3)
result.heat.crps <- matrix(NA, rep, 3)
result.heat.comptime <- matrix(NA, rep, 3)
colnames(result.heat.rmse) <- c("DNAmf", "Fractional BM", "BM")
colnames(result.heat.crps) <- c("DNAmf", "Fractional BM", "BM")
colnames(result.heat.comptime) <- c("DNAmf", "Fractional BM", "BM")
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
  y1 <- heat(X1,m1)
  y2 <- heat(X2,m2)
  y3 <- heat(X3,m3)
  y4 <- heat(X4,m4)
  y5 <- heat(X5,m5)
  
  
  ### DNAmf-sqex ###
  tic.DNAmf.sqex <- proc.time()[3]
  fit.DNAmf.sqex <- DNAmf.sqex(X1, y1,
                               X=rbind(X2, X3, X4, X5),
                               y=rbind(y2, y3, y4, y5),
                               nn=c(length(y1),length(y2),length(y3),length(y4),length(y5)),
                               t=c(m1,m2,m3,m4,m5), multi.start=10, constant=TRUE)
  pred.DNAmf.sqex <- predict.DNAmf(fit.DNAmf.sqex, x.test, targett=targett)
  predydiffu.sqex <- pred.DNAmf.sqex$mu
  predsig2diffu.sqex <- pred.DNAmf.sqex$sig2
  beta.sqex <- c(beta.sqex, fit.DNAmf.sqex$fit2$beta)
  delta.sqex <- c(delta.sqex, fit.DNAmf.sqex$fit2$delta)
  toc.DNAmf.sqex <- proc.time()[3]
  
  ### Fractional Brownian motion ###
  tic.fbm <- proc.time()[3]
  fit.fbm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5), 
                        rbind(matrix(rep(m1, nrow(X1))),matrix(rep(m2, nrow(X2))),matrix(rep(m3, nrow(X3))),
                              matrix(rep(m4, nrow(X4))),matrix(rep(m5, nrow(X5))) ), 
                        rbind(y1, y2, y3, y4, y5) )
  pred.fbm <- predict.MuFiMeshGP(fit.fbm, matrix(x.test), rep(0,length(x)))
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
  pred.bm <- predict.MuFiMeshGP(fit.bm, matrix(x.test), rep(0,length(x)))
  predybm <- pred.bm$mean
  predsig2bm <- pred.bm$sd
  toc.bm <- proc.time()[3]
  
  
  ### RMSE ###
  result.heat.rmse[i,1] <- sqrt(mean((predydiffu.sqex-y.test)^2)) #
  result.heat.rmse[i,4] <- sqrt(mean((predyfbm-y.test)^2)) #
  result.heat.rmse[i,5] <- sqrt(mean((predybm-y.test)^2)) #
  
  result.heat.crps[i,1] <- mean(crps(y.test, predydiffu.sqex, predsig2diffu.sqex)) #
  result.heat.crps[i,4] <- mean(crps(y.test, predyfbm, predsig2fbm)) #
  result.heat.crps[i,5] <- mean(crps(y.test, predybm, predsig2bm)) #
  
  result.heat.comptime[i,1] <- toc.DNAmf.sqex-tic.DNAmf.sqex
  result.heat.comptime[i,4] <- toc.fbm-tic.fbm
  result.heat.comptime[i,5] <- toc.bm-tic.bm
  
  print(result.heat.rmse[1:i,, drop=FALSE])
  print(apply(result.heat.rmse[1:i,, drop=FALSE], 2, mean))
  print(apply(result.heat.crps[1:i,, drop=FALSE], 2, mean))
}

print(mean(beta))
print(mean(delta))

