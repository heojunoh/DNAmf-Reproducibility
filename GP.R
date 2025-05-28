GP <- function(X, y, g = sqrt(.Machine$double.eps), constant = FALSE, p=0.1, min_cor = 0.2, max_cor = 0.5, init=NULL, lower=NULL, upper=NULL) { # p=0.05 for hetGP
  if (constant) {
    if (is.null(dim(X))) X <- matrix(X, ncol = 1)

    Xscaled <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% 
      diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X)) 
    if(is.null(lower)) lower <- -quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], p) / log(min_cor) * 
      (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2 
    if(is.null(upper)) upper <- -quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], 1-p) / log(max_cor) * 
      (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2 
    if(is.null(init)) init <- sqrt(lower * upper)

    n <- length(y)

    nlsep <- function(par, X, Y) {
      theta <- par # lengthscale
      K <- covar.sep(X, d = theta, g = g)
      Ki <- solve(K)
      ldetK <- determinant(K, logarithm = TRUE)$modulus

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- drop((t(one.vec) %*% Ki %*% Y) / (t(one.vec) %*% Ki %*% one.vec))

      tau2hat <- drop(t(Y - mu.hat) %*% Ki %*% (Y - mu.hat) / n)
      ll <- -(n / 2) * log(tau2hat) - (1 / 2) * ldetK
      return(drop(-ll))
    }

    gradnlsep <- function(par, X, Y) {
      theta <- par
      K <- covar.sep(X, d = theta, g = g)
      Ki <- solve(K)

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- drop((t(one.vec) %*% Ki %*% Y) / (t(one.vec) %*% Ki %*% one.vec))

      KiY <- Ki %*% (Y - mu.hat)
      ## loop over theta components
      dlltheta <- rep(NA, length(theta))
      for (k in 1:length(dlltheta)) {
        dotK <- K * distance(X[, k]) / (theta[k]^2)
        dlltheta[k] <- (n / 2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dotK))
      }
      return(-c(dlltheta))
    }

    outg <- optim(init, nlsep, gradnlsep, method = "L-BFGS-B", lower = lower, upper = upper, X = X, Y = y)
    
    theta <- outg$par
    K <- covar.sep(X, d = theta, g = g)
    Ki <- solve(K)
    one.vec <- matrix(1, ncol = 1, nrow = n)
    mu.hat <- drop((t(one.vec) %*% Ki %*% y) / (t(one.vec) %*% Ki %*% one.vec))
    tau2hat <- drop(t(y - mu.hat) %*% Ki %*% (y - mu.hat) / nrow(X))
    names(theta) <- NULL

    return(list(K = K, Ki = Ki, X = X, y = y, theta = theta, g = g, mu.hat = mu.hat, tau2hat = tau2hat, constant = constant))
  } else {
    if (is.null(dim(X))) X <- matrix(X, ncol = 1)

    Xscaled <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% 
      diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X)) 
    if(is.null(lower)) lower <- -quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], p) / log(min_cor) * 
      (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2 
    if(is.null(upper)) upper <- -quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], 1-p) / log(max_cor) * 
      (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2 
    if(is.null(init)) init <- sqrt(lower * upper)

    n <- length(y)

    nlsep <- function(par, X, Y) {
      theta <- par # lengthscale
      K <- covar.sep(X, d = theta, g = g)
      Ki <- solve(K)
      ldetK <- determinant(K, logarithm = TRUE)$modulus
      ll <- -(n / 2) * log(t(Y) %*% Ki %*% Y) - (1 / 2) * ldetK
      return(drop(-ll))
    }
    
    gradnlsep <- function(par, X, Y) {
      theta <- par
      K <- covar.sep(X, d = theta, g = g)
      Ki <- solve(K)
      
      KiY <- Ki %*% Y
      ## loop over theta components
      dlltheta <- rep(NA, length(theta))
      for (k in 1:length(dlltheta)) {
        dotK <- K * distance(X[, k]) / (theta[k]^2)
        dlltheta[k] <- (n / 2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dotK))
      }
      return(-c(dlltheta))
    }

    outg <- optim(init, nlsep, gradnlsep, method = "L-BFGS-B", lower = lower, upper = upper, X = X, Y = y)
    
    theta <- outg$par
    K <- covar.sep(X, d = theta, g = g)
    Ki <- solve(K)
    mu.hat <- 0
    tau2hat <- drop(t(y) %*% Ki %*% y / n)
    names(theta) <- NULL

    return(list(K = K, Ki = Ki, X = X, y = y, theta = theta, g = g, mu.hat = mu.hat, tau2hat = tau2hat, constant = constant))
  }
}


pred.GP <- function(fit, xnew, cov.out=FALSE) {
  xnew <- as.matrix(xnew)
  
  Ki <- fit$Ki
  theta <- fit$theta
  g <- fit$g
  X <- fit$X
  y <- fit$y
  tau2hat <- fit$tau2hat
  mu.hat <- fit$mu.hat
  
  KXX <- covar.sep(xnew, d = theta, g = g)
  KX <- covar.sep(xnew, X, d = theta, g = 0)
  
  mup2 <- mu.hat + KX %*% Ki %*% (y - mu.hat)
  Sigmap2 <- tau2hat * (KXX - KX %*% Ki %*% t(KX))
  Sigmap2 <- (t(Sigmap2)+Sigmap2)/2
  if(cov.out){
    return(list(mu = mup2, cov=Sigmap2, sig2 = pmax(0,diag(Sigmap2))))
  }else{
    return(list(mu = mup2, sig2 = pmax(0,diag(Sigmap2))))
  }
}
