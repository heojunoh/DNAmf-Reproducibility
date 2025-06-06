
checknested <- function(XX1, XX2) {
  checknest <- c()
  for (i in 1:nrow(XX2)) {
    checknest <- c(checknest, suppressWarnings(any(apply(XX1, 1, function(xx) {
      all.equal(XX2[i, ], xx, tolerance = sqrt(.Machine$double.eps))
    }))))
  }
  checknest[is.na(checknest)] <- FALSE
  all(checknest)
}

DNAmf.sqex <- function(X1, y1, X, y, nn, t, constant = TRUE, init=NULL, n.iter=500, multi.start=1, fitGP1=TRUE, ...) {
  
  g <- sqrt(.Machine$double.eps)
  lvl <- rep(seq_along(nn[-1]), times = nn[-1])  
  idxs  <- split(seq_len(nrow(X)), lvl)
  X_list <- lapply(idxs, function(i) X[i, , drop = FALSE])
  y_list <- lapply(idxs, function(i) y[i, , drop = FALSE])
  
  # check whether the design is nested #
  nested <- all(unlist(lapply(X_list, checknested, XX1=X1)))
  
  # make designs nested
  if(!nested) XX <- makenested(c(list(X1), X_list)) else XX <- X
  model <- DNAmf_sqex(X1, y1, XX, c(list(y1), y_list), t=t, nn=nn, nested = nested, constant = constant, fitGP1 = fitGP1,
                      init=init, multi.start = multi.start, n.iter = n.iter, trace = TRUE, g = g, burn.ratio = 0.75, ...)
  return(model)
}

DNAmf_sqex <- function(X1, y1, X_list, y_list, t, nn, nested = TRUE, constant = TRUE, fitGP1 = TRUE, init=NULL, 
                       multi.start = 1, n.iter = 500, trace = TRUE, g = sqrt(.Machine$double.eps), burn.ratio = 0.75, ...) {
  time0 <- proc.time()[3]
  
  if (nested) { # Nested
    X <- X_list
    y <- matrix(unlist(y_list[-1]))
    
    if (!checknested(X1, X)) stop("X is not nested by X1")
    
    tt <- rep(t[-1], times = nn[-1])
    
    if(fitGP1) fit1 <- GP(X1, y1, constant = constant, ...)
    if (2 < length(nn)){
      nnn <- c(1, nn[-1])
      yy <- y1[-(1:(nn[1]-nn[2])),, drop=FALSE]
      for(i in seq_len(length(t)-2)){
        yy <- rbind(yy, y[sum(nnn[1:i]):(sum(nnn[1:(i+1)])-1),, drop=FALSE][-(1:(nnn[i+1]-nnn[i+2])),, drop=FALSE])
      }    
      fit2 <- GP.nonsep.sqex(cbind(X, yy), y, tt, constant = constant, multi.start=multi.start, init=init, ...)
    } else {
      stop("Level should be larger than 2")
    }
    
    model <- list()
    if(fitGP1) model$fit1 <- fit1
    model$fit2 <- fit2
    model$constant <- constant
    model$t <- t
    model$TT <- tt
    model$nn <- nn
    model$level <- length(t)
    model$nested <- TRUE
  } else { # Non-nested
    XX    <- X_list
    y_all <- y_list
    L     <- length(XX$X_star)
    d     <- ncol(XX$X_star[[1]])
    
    t_list <- lapply(seq_len(L), function(i) rep(t[i], nrow(XX$X_list[[i]])))
    yy <- list(y_star  = vector("list", L), y_list  = y_all, y_tilde = vector("list", L))
    
    n.burnin <- ceiling(n.iter * burn.ratio) # number of burn in
    n.param <- n.iter - n.burnin # number of collected estimates
    param_mat <- matrix(NA, nrow=n.param, ncol=d+6) # theta_x, theta_y, theta_t, beta, delta, mu.hat, tau.hat
    
    ### Draw initial y tilde ###
    ### sample initial y1_tilde from the posterior ###
    fit1 <- GP(XX$X_list[[1]], y_all[[1]], g = g, constant = TRUE, init = NULL)
    pred1 <- pred.GP(fit1, XX$X_tilde[[1]], cov.out = TRUE) 
    yy$y_tilde[[1]] <- pred1$mu
    yy$y_star[[1]] <- rbind(yy$y_list[[1]],yy$y_tilde[[1]])
    
    ### sample initial y_2_tilde using individual GPs ###
    fit2 <- GP.nonsep.sqex(X = do.call(rbind,XX$X_list[-1]), y = do.call(rbind,yy$y_list[-1]), t = do.call(c,t_list[-1]), constant = TRUE, ...)
    for (l in 2:(L-1)) {
      pred2 <- pred.GP.nonsep(fit2, XX$X_tilde[[l]], rep(t[l], nrow(XX$X_tilde[[l]])))
      yy$y_tilde[[l]] <- pred2$mu
      yy$y_star[[l]] <- rbind(yy$y_list[[l]],yy$y_tilde[[l]])
    }
    yy$y_star[[L]] <- yy$y_list[[L]]
    
    ### initial estimates
    fit.DNAmf <- DNAmf.sqex(X1=XX$X_star[[1]], y1=yy$y_star[[1]],
                            X=do.call(rbind, XX$X_star[-1]),
                            y=do.call(rbind, yy$y_star[-1]),
                            nn=unlist(lapply(yy$y_star, length)), t=t, fitGP1=FALSE,
                            constant=constant, init = NULL, multi.start=10, ...)
    fit2 <- fit.DNAmf$fit2
    print(c(fit2$theta_x, fit2$theta_t, fit2$beta, fit2$delta, fit2$mu.hat, fit2$tau2hat))
    
    for (j in 1:n.iter) { # Imputation and Maximization
      if (trace) cat(j, '\n')
      
      # Imputation step; impute y tilde using ESS
      yy <- imputer(XX, yy, t, pred1, fit2)
      
      param.init <- c(fit2$theta_x, fit2$theta_t, fit2$beta, fit2$delta)
      # Maximization step; optimize parameters n.iter times
      fit.DNAmf <- DNAmf.sqex(X1=XX$X_star[[1]], y1=yy$y_star[[1]],
                              X=do.call(rbind, XX$X_star[-1]),
                              y=do.call(rbind, yy$y_star[-1]),
                              nn=unlist(lapply(yy$y_star, length)), t=t, fitGP1=FALSE,
                              constant=constant, init = param.init, ...)
      fit2 <- fit.DNAmf$fit2
      
      if(j > n.burnin){
        param_mat[j-n.burnin,] <- c(fit2$theta_x, fit2$theta_t, fit2$beta, fit2$delta, fit2$mu.hat, fit2$tau2hat)
      }
      if (trace) print(c(fit2$theta_x, fit2$theta_t, fit2$beta, fit2$delta, fit2$mu.hat, fit2$tau2hat))
    } # end of j for loop
    
    # average with 75% burn-in
    colnames(param_mat) <- c(paste0("theta_x", seq_len(d+1)), "theta_t", "beta", "delta", "mu_hat", "tau2_hat")
    final_params <- colMeans(param_mat)
    
    model <- fit.DNAmf 
    model$fit1 <- fit1
    model$fit2 <- fit2
    model$fit2$theta_x <- final_params[1:(d+1)]
    model$fit2$theta_t <- final_params[d+2]
    model$fit2$beta <- final_params[d+3]
    model$fit2$delta <- final_params[d+4]
    model$fit2$mu.hat <- final_params[d+5]
    model$fit2$tau2hat <- final_params[d+6]
    
    model$XX <- XX
    model$yy <- yy
    model$t <- t
    model$pred1 <- pred1
    model$estim <- final_params
    model$nested <- FALSE
  }
  model$time <- proc.time()[3] - time0
  class(model) <- "DNAmf"
  return(model)
}
