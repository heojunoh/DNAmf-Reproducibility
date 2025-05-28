### Plate simulation ###
plate <- function(xx, t){
  d1 <- data.frame(xx, rep(t, nrow(xx))) # 
  write.csv(d1, "/Plate/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("/Plate/Solve_plate.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE, intern = TRUE)
  d2 <- read.table("/Plate/generate_text/temp_to_r.txt", sep = ",")
  y <- d2$V5
  matrix(y)
}

### training data ###
m1 <- 0.55; m2 <- 0.50; m3 <- 0.45; m4 <- 0.40; m5 <- 0.35; targett <- 0.3; 
n1 <- 25; n2 <- 20; n3 <- 15; n4 <- 10; n5 <- 5;
d <- 3 

### test data ###
set.seed(0)
x.test <- maximinLHS(100, d)
y.test <- plate(x.test, targett)

rep <- 100
result.plate.rmse <- matrix(NA, rep, 3)
result.plate.crps <- matrix(NA, rep, 3)
result.plate.comptime <- matrix(NA, rep, 3)
colnames(result.plate.rmse) <- c("DNAmf", "Fractional BM", "BM")
colnames(result.plate.crps) <- c("DNAmf", "Fractional BM", "BM")
colnames(result.plate.comptime) <- c("DNAmf", "Fractional BM", "BM")
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
    y1 <- plate(X1, m1)
    y2 <- plate(X2, m2)
    y3 <- plate(X3, m3)
    y4 <- plate(X4, m4)
    y5 <- plate(X5, m5)
    
    ### DNAmf - sqex ###
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
                          rbind(matrix(rep(m1, nrow(y1))),matrix(rep(m2, nrow(y2))),
                                matrix(rep(m3, nrow(y3))),matrix(rep(m4, nrow(y4))),matrix(rep(m5, nrow(y5)))),
                          rbind(y1, y2, y3, y4, y5))
    pred.fbm <- predict.MuFiMeshGP(fit.fbm, x.test, rep(0, nrow(x.test)))
    toc.fbm <- proc.time()[3]
    
    ### Brownian motion - Tuo Wu Yu ###
    tic.bm <- proc.time()[3]
    fit.bm <- MuFiMeshGP(rbind(X1, X2, X3, X4, X5), 
                         rbind(matrix(rep(m1, nrow(y1))),matrix(rep(m2, nrow(y2))),
                               matrix(rep(m3, nrow(y3))),matrix(rep(m4, nrow(y4))),matrix(rep(m5, nrow(y5)))),
                         rbind(y1, y2, y3, y4, y5), H.known = 0.5)
    pred.bm <- predict.MuFiMeshGP(fit.bm, x.test, rep(0, nrow(x.test)))
    toc.bm <- proc.time()[3]
    
    ### RMSE ###
    result.plate.rmse[i,1] <- sqrt(mean((predydiffu.sqex-y.test)^2)) #
    result.plate.rmse[i,4] <- sqrt(mean((pred.fbm$mean - y.test)^2))
    result.plate.rmse[i,5] <- sqrt(mean((pred.bm$mean - y.test)^2))
    
    result.plate.crps[i,1] <- mean(crps(y.test, predydiffu.sqex, predsig2diffu.sqex)) #
    result.plate.crps[i,4] <- mean(crps(y.test, pred.fbm$mean, pred.fbm$sd))
    result.plate.crps[i,5] <- mean(crps(y.test, pred.bm$mean, pred.bm$sd))
    
    result.plate.comptime[i,1] <- toc.DNAmf.sqex-tic.DNAmf.sqex
    result.plate.comptime[i,4] <- toc.fbm - tic.fbm
    result.plate.comptime[i,5] <- toc.bm - tic.bm
    
  }, error=function(e){})
  
  boxplot(result.plate.rmse[1:i,,drop=FALSE])
  boxplot(result.plate.crps[1:i,,drop=FALSE])
  print(apply(result.plate.rmse[1:i,,drop=FALSE], 2, mean))
  print(apply(result.plate.crps[1:i,,drop=FALSE], 2, mean))
  
  # save.image("plate comparison.RData")
}

print(mean(beta))
print(mean(delta))

