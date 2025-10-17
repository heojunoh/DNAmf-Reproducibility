### Plate simulation ###
# plate <- function(xx, t){
#   d1 <- data.frame(xx, rep(t, nrow(xx))) # 
#   write.csv(d1, "Plate/generate_text/temp_to_matlab.txt", row.names=F)
#   run_matlab_script("Plate/Solve_plate.m", verbose = FALSE, desktop = FALSE, 
#                     splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE, intern = TRUE)
#   d2 <- read.table("Plate/generate_text/temp_to_r.txt", sep = ",")
#   y <- d2$V5
#   matrix(y)
# }

### training data ###
c <- 0.9; bet <- -3
n5 <- 8; n4 <- round(c^bet*n5); n3 <- round(c^(2*bet)*n5); n2 <- round(c^(3*bet)*n5); n1 <- round(c^(4*bet)*n5);

m1 <- 0.55; m2 <- m1*c; m3 <- m2*c; m4 <- m3*c; m5 <- m4*c; 
d <- 3; targett=0.3

### test data ###
set.seed(0)
x.test <- maximinLHS(100, d)
# y.test <- plate(x.test, targett)
y.test <- matrix(read.table("Plate/generate_text/ytest.txt", sep = ",")$V5)

rep <- 100
result.plate.rmse <- result.plate.crps <- matrix(NA, rep, 6)
colnames(result.plate.rmse) <- colnames(result.plate.crps) <- c("DNAmf", "RNAmf", "Cokriging", "NARGP", "Fractional BM", "BM")
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
    y1 <- matrix(read.table(file.path("Plate/generate_text", sprintf("temp_to_r_X1_%d.txt", i)), sep = ",")$V5)
    y2 <- matrix(read.table(file.path("Plate/generate_text", sprintf("temp_to_r_X2_%d.txt", i)), sep = ",")$V5)
    y3 <- matrix(read.table(file.path("Plate/generate_text", sprintf("temp_to_r_X3_%d.txt", i)), sep = ",")$V5)
    y4 <- matrix(read.table(file.path("Plate/generate_text", sprintf("temp_to_r_X4_%d.txt", i)), sep = ",")$V5)
    y5 <- matrix(read.table(file.path("Plate/generate_text", sprintf("temp_to_r_X5_%d.txt", i)), sep = ",")$V5)
    
    if(NARGP) saveRDS(list(X1=X3, X2=X4, X3=X5, Y1=y3, Y2=y4, Y3=y5, Xtest=x.test, Ytest=y.test), file = "tmp_data.rds")
    
    ### DNAmf - sqex ###
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
                          rbind(matrix(rep(m1, nrow(y1))),matrix(rep(m2, nrow(y2))),
                                matrix(rep(m3, nrow(y3))),matrix(rep(m4, nrow(y4))),matrix(rep(m5, nrow(y5)))),
                          rbind(y1, y2, y3, y4, y5))
    pred.fbm <- predict(fit.fbm, x.test, rep(0, nrow(x.test)))
    predyfbm <- pred.fbm$mean
    predsig2fbm <- pred.fbm$sd
    toc.fbm <- proc.time()[3]
    
    ### Brownian motion - Tuo Wu Yu ###
    tic.bm <- proc.time()[3]
    fit.bm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5), 
                         rbind(matrix(rep(m1, nrow(y1))),matrix(rep(m2, nrow(y2))),
                               matrix(rep(m3, nrow(y3))),matrix(rep(m4, nrow(y4))),matrix(rep(m5, nrow(y5)))),
                         rbind(y1, y2, y3, y4, y5), H.known = 0.5)
    pred.bm <- predict(fit.bm, x.test, rep(0, nrow(x.test)))
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
    result.plate.rmse[i,1] <- sqrt(mean((predydiffu-y.test)^2)) #
    result.plate.rmse[i,2] <- sqrt(mean((predy.RNAmf-y.test)^2)) #
    result.plate.rmse[i,3] <- sqrt(mean((pred.muficokm$mean-y.test)^2)) #
    if(NARGP) result.plate.rmse[i,4] <- py$error
    result.plate.rmse[i,5] <- sqrt(mean((predyfbm-y.test)^2)) #
    result.plate.rmse[i,6] <- sqrt(mean((predybm-y.test)^2)) #
    
    result.plate.crps[i,1] <- mean(crps(y.test, predydiffu, predsig2diffu)) #
    result.plate.crps[i,2] <- mean(crps(y.test, predy.RNAmf, predsig2.RNAmf)) #
    result.plate.crps[i,3] <- mean(crps(y.test, pred.muficokm$mean, pred.muficokm$sig2)) #
    if(NARGP) result.plate.crps[i,4] <- py$crps
    result.plate.crps[i,5] <- mean(crps(y.test, predyfbm, predsig2fbm)) #
    result.plate.crps[i,6] <- mean(crps(y.test, predybm, predsig2bm)) #
  }, error=function(e){})
  
  # boxplot(result.plate.rmse[1:i,,drop=FALSE])
  # boxplot(result.plate.crps[1:i,,drop=FALSE])
  print(result.plate.rmse[1:i,,drop=FALSE])
  print(result.plate.crps[1:i,,drop=FALSE])
}

print(mean(beta))
print(mean(delta))

