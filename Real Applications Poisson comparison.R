### Poisson's equation ###
f_true_0 <- function(x) {
  c <- 2 * x - 1  # Transform x
  h <- function(b) exp(c * b) * sin(pi * b)  # Define h(b)
  result <- optimize(h, interval = c(0, 1), maximum = TRUE)  # Find maximum
  return(result$objective)  # Return the maximum value
}

# poisson <- function(xx, mesh){
#   if(is.null(dim(xx))) n <- length(xx) else n <- nrow(xx)
#   d1 <- data.frame(2*xx-1, rep(mesh, nrow(xx))) # scale X to [-1,1]
#   write.csv(xx, "Poisson-equation/generate_text/temp_to_X.txt", row.names=F)
#   write.csv(d1, "Poisson-equation/generate_text/temp_to_matlab.txt", row.names=F)
#   run_matlab_script("Poisson-equation/SolvePoissonEqn.m", verbose = FALSE, desktop = FALSE, 
#                     splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
#                     intern = TRUE)
#   d2 <- read.table("Poisson-equation/generate_text/temp_to_r.txt", sep = ",")
#   out <- matrix(as.numeric(unlist(d2)), nrow = n)
#   return(matrix(out[,3], ncol=1)) # row-wise output
# }

### training data ###
c <- 0.65; gam <- -1
n5 <- 3; n4 <- round(c^gam*n5); n3 <- round(c^(2*gam)*n5); n2 <- round(c^(3*gam)*n5); n1 <- round(c^(4*gam)*n5);

m1 <- 0.1; m2 <- m1*c; m3 <- m2*c; m4 <- m3*c; m5 <- m4*c; 
d <-1; targett=0

### test data ###
x.test <- matrix(seq(0,1,length=100),ncol=1)
y.test <- matrix(sapply(x.test,f_true_0),ncol=1)

rep <- 100
result.poisson.rmse <- result.poisson.crps <- matrix(NA, rep, 6)
colnames(result.poisson.rmse) <- colnames(result.poisson.crps) <- c("DNAmf", "RNAmf", "Cokriging", "NARGP", "Fractional BM", "BM")
beta <- delta <- c()

par(mfrow=c(1,2))
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
    y1 <- matrix(read.table(file.path("Poisson-equation/generate_text", sprintf("temp_to_r_X1_%d.txt", i)), sep = ",")$V3)
    y2 <- matrix(read.table(file.path("Poisson-equation/generate_text", sprintf("temp_to_r_X2_%d.txt", i)), sep = ",")$V3)
    y3 <- matrix(read.table(file.path("Poisson-equation/generate_text", sprintf("temp_to_r_X3_%d.txt", i)), sep = ",")$V3)
    y4 <- matrix(read.table(file.path("Poisson-equation/generate_text", sprintf("temp_to_r_X4_%d.txt", i)), sep = ",")$V3)
    y5 <- matrix(read.table(file.path("Poisson-equation/generate_text", sprintf("temp_to_r_X5_%d.txt", i)), sep = ",")$V3)
    
    if(NARGP) saveRDS(list(X1=X3, X2=X4, X3=X5, Y1=y3, Y2=y4, Y3=y5, Xtest=x.test, Ytest=y.test), file = "tmp_data.rds")
    
    ### DNAmf-sqex ###
    tic.DNAmf <- proc.time()[3]
    fit.DNAmf <- DNAmf.sqex(X1, y1,
                            X=rbind(X2, X3, X4, X5),
                            y=rbind(y2, y3, y4, y5),
                            nn=c(length(y1),length(y2),length(y3),length(y4),length(y5)),
                            t=c(m1,m2,m3,m4,m5), multi.start=50, 
                            constant=TRUE)
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
                                matrix(rep(m4, nrow(X4))), matrix(rep(m5, nrow(X5))) ), 
                          rbind(y1, y2, y3, y4, y5) )
    pred.fbm <- predict(fit.fbm, x.test, matrix(rep(0,nrow(x.test))) )
    predyfbm <- pred.fbm$mean
    predsig2fbm <- pred.fbm$sd
    toc.fbm <- proc.time()[3]
    
    ### Brownian motion - Tuo Wu Yu ###
    tic.bm <- proc.time()[3]
    fit.bm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5), 
                         rbind(matrix(rep(m1, nrow(X1))),matrix(rep(m2, nrow(X2))),matrix(rep(m3, nrow(X3))),
                               matrix(rep(m4, nrow(X4))), matrix(rep(m5, nrow(X5))) ), 
                         rbind(y1, y2, y3, y4, y5), H.known = 0.5)
    pred.bm <- predict(fit.bm, x.test, matrix(rep(0,nrow(x.test))) )
    predybm <- pred.bm$mean
    predsig2bm <- pred.bm$sd
    toc.bm <- proc.time()[3]
    
    ### Cokriging ###
    NestDesignMuFiCokriging <- NestedDesignBuild(design = list(X3,X4,X5))
    fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesignMuFiCokriging, covtype="gauss",
                             lower=eps, upper=0.3,
                             response = list(y3,y4,y5), nlevel = 3)
    pred.muficokm <- predict(fit.muficokm, x.test, "SK")
    
    ### NARGP ###
    if(NARGP) py <- py_run_file("NARGP.py")
    
    
    ### RMSE ###
    result.poisson.rmse[i,1] <- sqrt(mean((predydiffu-y.test)^2)) #
    result.poisson.rmse[i,2] <- sqrt(mean((predy.RNAmf-y.test)^2)) #
    result.poisson.rmse[i,3] <- sqrt(mean((pred.muficokm$mean-y.test)^2)) #
    if(NARGP) result.poisson.rmse[i,4] <- py$error
    result.poisson.rmse[i,5] <- sqrt(mean((predyfbm-y.test)^2)) #
    result.poisson.rmse[i,6] <- sqrt(mean((predybm-y.test)^2)) #
    
    result.poisson.crps[i,1] <- mean(crps(y.test, predydiffu, predsig2diffu)) #
    result.poisson.crps[i,2] <- mean(crps(y.test, predy.RNAmf, predsig2.RNAmf)) #
    result.poisson.crps[i,3] <- mean(crps(y.test, pred.muficokm$mean, pred.muficokm$sig2)) #
    if(NARGP) result.poisson.crps[i,4] <- py$crps
    result.poisson.crps[i,5] <- mean(crps(y.test, predyfbm, predsig2fbm)) #
    result.poisson.crps[i,6] <- mean(crps(y.test, predybm, predsig2bm)) #
  }, error=function(e){})
  
  # boxplot(result.poisson.rmse[1:i,,drop=FALSE])
  # boxplot(result.poisson.crps[1:i,,drop=FALSE])
  print(result.poisson.rmse[1:i,,drop=FALSE])
  print(result.poisson.crps[1:i,,drop=FALSE])
}

print(mean(beta))
print(mean(delta))

