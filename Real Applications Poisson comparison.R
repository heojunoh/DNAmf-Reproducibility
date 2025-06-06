### Poisson's equation ###
f_true_0 <- function(x) {
  c <- 2 * x - 1  # Transform x
  h <- function(b) exp(c * b) * sin(pi * b)  # Define h(b)
  result <- optimize(h, interval = c(0, 1), maximum = TRUE)  # Find maximum
  return(result$objective)  # Return the maximum value
}

poisson <- function(xx, mesh){
  if(is.null(dim(xx))) n <- length(xx) else n <- nrow(xx)
  d1 <- data.frame(2*xx-1, rep(mesh, nrow(xx))) # scale X to [-1,1]
  write.csv(xx, "/Poisson-equation/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "/Poisson-equation/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("/Poisson-equation/SolvePoissonEqn.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("/Poisson-equation/generate_text/temp_to_r.txt", sep = ",")
  out <- matrix(as.numeric(unlist(d2)), nrow = n)
  return(matrix(out[,3], ncol=1)) # row-wise output
}

### training data ###
n1 <- 15; n2 <- 12; n3 <- 9; n4 <- 6; n5 <- 3
m1 <- 0.09; m2 <- 0.07; m3 <- 0.05; m4 <- 0.03; m5 <- 0.01; targett=0
d <-1

### test data ###
x.test <- matrix(seq(0,1,length=100),ncol=1)
y.test <- matrix(sapply(x.test,f_true_0),ncol=1)

rep <- 100
result.poisson.rmse <- matrix(NA, rep, 3)
result.poisson.crps <- matrix(NA, rep, 3)
result.poisson.comptime <- matrix(NA, rep, 3)
colnames(result.poisson.rmse) <- c("DNAmf", "Fractional BM", "BM")
colnames(result.poisson.crps) <- c("DNAmf", "Fractional BM", "BM")
colnames(result.poisson.comptime) <- c("DNAmf", "Fractional BM", "BM")
beta <- delta <- c()

par(mfrow=c(1,2))
### comparison using 3 PCs ###
for(i in 1:rep) {
  tryCatch({
    set.seed(i)
    print(i)
    
    NestDesign <- NestedX(c(n1,n2,n3,n4,n5),d)
    X1 <- NestDesign[[1]]
    X2 <- NestDesign[[2]]
    X3 <- NestDesign[[3]]
    X4 <- NestDesign[[4]]
    X5 <- NestDesign[[5]]
    y1 <- poisson(X1, mesh=m1)
    y2 <- poisson(X2, mesh=m2)
    y3 <- poisson(X3, mesh=m3)
    y4 <- poisson(X4, mesh=m4)
    y5 <- poisson(X5, mesh=m5)
    
    
    ### DNAmf-sqex ###
    tic.DNAmf.sqex <- proc.time()[3]
    fit.DNAmf.sqex <- DNAmf.sqex(X1, y1,
                                 X=rbind(X2, X3, X4, X5),
                                 y=rbind(y2, y3, y4, y5),
                                 nn=c(length(y1),length(y2),length(y3),length(y4),length(y5)),
                                 t=c(m1,m2,m3,m4,m5), multi.start=50, 
                                 constant=TRUE)
    pred.DNAmf.sqex <- predict.DNAmf(fit.DNAmf.sqex, x.test, targett=0)
    predydiffu.sqex <- pred.DNAmf.sqex$mu
    predsig2diffu.sqex <- pred.DNAmf.sqex$sig2
    beta.sqex <- c(beta.sqex, fit.DNAmf.sqex$fit2$beta)
    delta.sqex <- c(delta.sqex, fit.DNAmf.sqex$fit2$delta)
    toc.DNAmf.sqex <- proc.time()[3]
    
    ### Fractional Brownian motion ###
    tic.fbm <- proc.time()[3]
    fit.fbm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5), 
                          rbind(matrix(rep(m1, nrow(X1))),matrix(rep(m2, nrow(X2))),matrix(rep(m3, nrow(X3))),
                                matrix(rep(m4, nrow(X4))), matrix(rep(m5, nrow(X5))) ), 
                          rbind(y1, y2, y3, y4, y5) )
    pred.fbm <- predict.MuFiMeshGP(fit.fbm, x.test, matrix(rep(0,nrow(x.test))) )
    predyfbm <- pred.fbm$mean
    predsig2fbm <- pred.fbm$sd
    toc.fbm <- proc.time()[3]
    
    ### Brownian motion - Tuo Wu Yu ###
    tic.bm <- proc.time()[3]
    fit.bm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5), 
                         rbind(matrix(rep(m1, nrow(X1))),matrix(rep(m2, nrow(X2))),matrix(rep(m3, nrow(X3))),
                               matrix(rep(m4, nrow(X4))), matrix(rep(m5, nrow(X5))) ), 
                         rbind(y1, y2, y3, y4, y5), 
                         H.known = 0.5)
    pred.bm <- predict.MuFiMeshGP(fit.bm, x.test, matrix(rep(0,nrow(x.test))) )
    predybm <- pred.bm$mean
    predsig2bm <- pred.bm$sd
    toc.bm <- proc.time()[3]
    
    
    ### RMSE ###
    result.poisson.rmse[i,1] <- sqrt(mean((predydiffu.sqex-y.test)^2)) #
    result.poisson.rmse[i,2] <- sqrt(mean((predyfbm-y.test)^2)) #
    result.poisson.rmse[i,3] <- sqrt(mean((predybm-y.test)^2)) #
    
    result.poisson.crps[i,1] <- mean(crps(y.test, predydiffu.sqex, predsig2diffu.sqex)) #
    result.poisson.crps[i,2] <- mean(crps(y.test, predyfbm, predsig2fbm)) #
    result.poisson.crps[i,3] <- mean(crps(y.test, predybm, predsig2bm)) #
    
    result.poisson.comptime[i,1] <- toc.DNAmf.sqex-tic.DNAmf.sqex
    result.poisson.comptime[i,2] <- toc.fbm-tic.fbm
    result.poisson.comptime[i,3] <- toc.bm-tic.bm
  }, error=function(e){})
  
  boxplot(result.poisson.rmse[1:i,,drop=FALSE])
  boxplot(result.poisson.crps[1:i,,drop=FALSE])
  print(result.poisson.rmse[1:i,,drop=FALSE])
  print(result.poisson.crps[1:i,,drop=FALSE])
  
  # save.image("Poisson comparison.RData")
}

print(mean(beta))
print(mean(delta))

