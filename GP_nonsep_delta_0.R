theta_bounds <- function(X, t, p = 0.05, beta = 0.5, min_cor = 0.1, max_cor = 0.9) {

  Xscaled <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*%
    diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X))
  tscaled <- t / max(t)

  repr_dist_x_low <- quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], p)
  repr_dist_x_high <- quantile(distance(Xscaled)[lower.tri(distance(Xscaled))], 1-p)
  repr_dist_t <- quantile(distance(tscaled)[lower.tri(distance(tscaled))], 0.5)

  d <- ncol(X)

  tmpfun <- function(theta, repr_dist_t, repr_dist_x, beta, d, value) {
    theta_x <- theta[1]
    theta_t <- theta[2]
    delta <- 0.5 # or 1?

    term1 <- (repr_dist_t^2 / theta_t + 1)^(- (beta * d) / 2 - delta)
    term2 <- exp(- (repr_dist_t^2 / theta_t + 1)^(-beta) * (repr_dist_x^2 / theta_x))

    return(term1 * term2 - value)
  }

  theta_min <- tryCatch(
    stats::uniroot(
      function(theta) tmpfun(c(theta, theta), repr_dist_t, repr_dist_x_low, beta, d, min_cor),
      interval = c(sqrt(.Machine$double.eps), 100)
    )$root,
    error = function(e) {
      warning("Lower bound estimation failed. Using default values.")
      return(c(1e-2, 1e-2))
    }
  )
  theta_max <- tryCatch(
    stats::uniroot(
      function(theta) tmpfun(c(theta, theta), repr_dist_t, repr_dist_x_high, beta, d, max_cor),
      interval = c(sqrt(.Machine$double.eps), 100)
    )$root,
    error = function(e) {
      warning("Upper bound estimation failed. Using default values.")
      return(c(5, 5))
    }
  )

  return(list(lower = theta_min, upper = theta_max))
}

cor.sep.sqex <- function(X, x = NULL, t, tnew=NULL, param, alpha=1) {
  d <- NCOL(X)
  n <- NROW(X)
  if (length(t) != n) {
    stop("Length of mesh size should be the number of rows of inputs")
  }

  theta_x <- param[1:d]
  theta_t <- param[d+1]
  beta <- param[d+2]
  delta <- param[d+3]

  if (is.null(x)) {
    if(length(unique(t))==1){
      Tcomponent=1
    }else{
      Tcomponent <- 1/(outer(t, t, "-")^(2*alpha)/theta_t + 1)
    }
    K <- Tcomponent^(beta*d/2+delta) * covar.sep(X, d=theta_x, g=0)^(Tcomponent^beta)
  } else {
    n.new <- NROW(x)
    if (length(tnew) != n.new) {
      stop("Length of mesh size should be the number of rows of inputs")
    }
    if(length(unique(t))==1){
      Tcomponent=1
    }else{
      Tcomponent <- 1/(outer(t, tnew, "-")^(2*alpha)/theta_t + 1)
    }
    K <- Tcomponent^(beta*d/2+delta) * covar.sep(X, x, d=theta_x, g=0)^(Tcomponent^beta)
  }
  return(K)
}

GP.nonsep.sqex <- function(X, y, t, g = sqrt(.Machine$double.eps), constant = FALSE,
                           p=0.1, min_cor = 0.2, max_cor = 0.5,
                           init=NULL, lower=NULL, upper=NULL, multi.start=1,
                           alpha=1) { # p=0.05 for hetGP
  if (constant) {
    if (is.null(dim(X))) X <- matrix(X, ncol = 1)

    theta_bounds <- theta_bounds(X, t, p = p, min_cor = min_cor, max_cor = max_cor)

    if(is.null(lower)) {
      Xlower <- theta_bounds$lower * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2
      tlower <- 1e-2 * max(t) # theta_bounds_t$lower * (max(max(t)-min(t), eps))^2
      lower <- c(Xlower, tlower, eps, 0)
    }
    if(is.null(upper)) {
      Xupper <- theta_bounds$upper * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2
      tupper <- 1e1 * max(t) # theta_bounds_t$upper * (max(max(t)-min(t), eps))^2
      upper <- c(Xupper, tupper, 1-eps, eps)
    }
    if(is.null(init)) {
      init <- c(sqrt(lower*upper)[-c(length(lower)-1, length(lower))], 0.5, 1)
    }

    # lower[init < lower] <- init[init < lower]
    # upper[init > upper] <- init[init > upper]

    n <- nrow(X)

    nlsep <- function(par, X, Y, tt, alpha) {
      K <- cor.sep.sqex(X, t=tt, param=par, alpha=alpha)
      Ki <- solve(K + diag(g, n))
      ldetK <- determinant(K + diag(g, n), logarithm = TRUE)$modulus

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- drop((t(one.vec) %*% Ki %*% Y) / (t(one.vec) %*% Ki %*% one.vec))

      tau2hat <- drop(t(Y - mu.hat) %*% Ki %*% (Y - mu.hat) / n)

      ll <- -(n / 2) * log(tau2hat) - (1 / 2) * ldetK

      return(drop(-ll))
    }

    gradnlsep <- function(par, X, Y, tt, alpha) {
      d <-ncol(X)
      dd <- d-1
      n <- nrow(X)
      theta <- par[1:d]
      theta_t <- par[d+1]
      beta <- par[d+2]
      delta <- par[d+3]

      K <- cor.sep.sqex(X, t=tt, param=par, alpha=alpha)
      Ki <- solve(K + diag(g, n))

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- drop((t(one.vec) %*% Ki %*% Y) / (t(one.vec) %*% Ki %*% one.vec))

      KiY <- Ki %*% (Y - mu.hat)

      Rt <- sqrt(distance(tt / sqrt(theta_t)))
      u <- (1 + Rt^2)^(-1)
      v <- matrix(rowSums(sapply(1:d, function(i) distance(X[, i])/theta[i])), ncol=n)


      du_dtheta_t <- (Rt^2 / theta_t) * u^2
      dk_dbeta <- K * log(u) * ((dd+1)/2 - u^beta*v)
      dk_ddelta <- K * log(u)

      dlltheta <- rep(NA, length(par))

      dk_du <- ((dd+1)/2 * beta + delta - beta*u^beta*v) * u^((dd+1)/2*beta+delta-1) * exp(-u^beta * v)
      dk_dv <- -u^((dd+3)/2*beta+delta) * exp(-u^beta * v)

      dk_dtheta_t <- dk_du * du_dtheta_t

      for (i in 1:d) {
        dk_dtheta_i <- -distance(X[,i] / theta[i]) * dk_dv
        dlltheta[i] <- (n / 2) * t(KiY) %*% dk_dtheta_i %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_i))
      }

      dlltheta[d+1] <- (n / 2) * t(KiY) %*% dk_dtheta_t %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_t))
      dlltheta[d+2] <- (n / 2) * t(KiY) %*% dk_dbeta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dbeta))
      dlltheta[d+3] <- (n / 2) * t(KiY) %*% dk_ddelta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_ddelta))

      return(-c(dlltheta))
    }

    if(multi.start > 1){
      start <- randomLHS(multi.start - 1, ncol(X)+3)
      start <- t(t(start) * (upper - lower) + lower)
      start <- rbind(init, start)
      for(i in 1:nrow(start)) {
        outi <- optim(start[i,], nlsep, gradnlsep,
                      method = "L-BFGS-B", lower = lower, upper = upper,
                      X = X, Y = y, tt = t, alpha=alpha)
        if(i == 1) {
          out <- outi
        }else if(outi$value < out$value) {
          out <- outi
        }
      }
    }else{
      out <- optim(init, nlsep, gradnlsep,
                   method = "L-BFGS-B", lower = lower, upper = upper,
                   X = X, Y = y, tt = t, alpha=alpha)
    }

    theta_x <- out$par[1:ncol(X)]
    theta_t <- out$par[ncol(X)+1]
    beta <- out$par[ncol(X)+2]
    delta <- out$par[ncol(X)+3]
    K <- cor.sep.sqex(X, t=t, param=out$par, alpha=alpha)
    Ki <- solve(K + diag(g, n))
    one.vec <- matrix(1, ncol = 1, nrow = n)
    mu.hat <- drop((t(one.vec) %*% Ki %*% y) / (t(one.vec) %*% Ki %*% one.vec))
    tau2hat <- drop(t(y - mu.hat) %*% Ki %*% (y - mu.hat) / nrow(X))

    return(list(K = K, Ki = Ki, X = X, y = y, t=t, theta_x=theta_x, theta_t=theta_t, beta=beta, delta=delta, alpha=alpha, g = g, mu.hat = mu.hat, tau2hat = tau2hat, constant = constant))
  } else {
    if (is.null(dim(X))) X <- matrix(X, ncol = 1)

    theta_bounds <- theta_bounds(X, t, p = p, min_cor = min_cor, max_cor = max_cor)

    if(is.null(lower)) {
      Xlower <- theta_bounds$lower * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2
      tlower <- 1e-2 * max(t) # theta_bounds_t$lower * (max(max(t)-min(t), eps))^2
      lower <- c(Xlower, tlower, eps, 0.5)
    }
    if(is.null(upper)) {
      Xupper <- theta_bounds$upper * (apply(X, 2, range)[2,] - apply(X, 2, range)[1,])^2
      tupper <- 1e1 * max(t) # theta_bounds_t$upper * (max(max(t)-min(t), eps))^2
      upper <- c(Xupper, tupper, 1-eps, 2-eps)
    }
    if(is.null(init)) {
      init <- c(sqrt(lower*upper)[-c(length(lower)-1, length(lower))], 0.5, 1)
    }

    lower[init < lower] <- init[init < lower]
    upper[init > upper] <- init[init > upper]

    n <- nrow(X)

    nlsep <- function(par, X, Y, tt, alpha) {
      K <- cor.sep.sqex(X, t=tt, param=par, alpha=alpha)
      Ki <- solve(K + diag(g, n))
      ldetK <- determinant(K + diag(g, n), logarithm = TRUE)$modulus

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- 0

      tau2hat <- drop(t(Y - mu.hat) %*% Ki %*% (Y - mu.hat) / n)

      ll <- -(n / 2) * log(tau2hat) - (1 / 2) * ldetK

      return(drop(-ll))
    }

    gradnlsep <- function(par, X, Y, tt, alpha) {
      d <-ncol(X)
      dd <- d-1
      n <- nrow(X)
      theta <- par[1:d]
      theta_t <- par[d+1]
      beta <- par[d+2]
      delta <- par[d+3]

      K <- cor.sep.sqex(X, t=tt, param=par, alpha=alpha)
      Ki <- solve(K + diag(g, n))

      one.vec <- matrix(1, ncol = 1, nrow = n)
      mu.hat <- 0

      KiY <- Ki %*% (Y - mu.hat)

      Rt <- sqrt(distance(tt / sqrt(theta_t)))
      u <- (1 + Rt^2)^(-1)
      v <- matrix(rowSums(sapply(1:d, function(i) distance(X[, i])/theta[i])), ncol=n)


      du_dtheta_t <- (Rt^2 / theta_t) * u^2
      dk_dbeta <- K * log(u) * ((dd+1)/2 - u^beta*v)
      dk_ddelta <- K * log(u)

      dlltheta <- rep(NA, length(par))

      dk_du <- ((dd+1)/2 * beta + delta - beta*u^beta*v) * u^((dd+1)/2*beta+delta-1) * exp(-u^beta * v)
      dk_dv <- -u^((dd+3)/2*beta+delta) * exp(-u^beta * v)

      dk_dtheta_t <- dk_du * du_dtheta_t

      for (i in 1:d) {
        dk_dtheta_i <- -distance(X[,i] / theta[i]) * dk_dv
        dlltheta[i] <- (n / 2) * t(KiY) %*% dk_dtheta_i %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_i))
      }

      dlltheta[d+1] <- (n / 2) * t(KiY) %*% dk_dtheta_t %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dtheta_t))
      dlltheta[d+2] <- (n / 2) * t(KiY) %*% dk_dbeta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_dbeta))
      dlltheta[d+3] <- (n / 2) * t(KiY) %*% dk_ddelta %*% KiY / (t(Y) %*% KiY) - (1 / 2) * sum(diag(Ki %*% dk_ddelta))

      return(-c(dlltheta))
    }

    if(multi.start > 1){
      start <- randomLHS(multi.start - 1, ncol(X)+3)
      start <- t(t(start) * (upper - lower) + lower)
      start <- rbind(init, start)
      for(i in 1:nrow(start)) {
        outi <- optim(start[i,], nlsep, gradnlsep,
                      method = "L-BFGS-B", lower = lower, upper = upper,
                      X = X, Y = y, tt = t, alpha=alpha)
        if(i == 1) {
          out <- outi
        }else if(outi$value < out$value) {
          out <- outi
        }
      }
    }else{
      out <- optim(init, nlsep, gradnlsep,
                   method = "L-BFGS-B", lower = lower, upper = upper,
                   X = X, Y = y, tt = t, alpha=alpha)
    }

    theta_x <- out$par[1:ncol(X)]
    theta_t <- out$par[ncol(X)+1]
    beta <- out$par[ncol(X)+2]
    delta <- out$par[ncol(X)+3]
    K <- cor.sep.sqex(X, t=t, param=out$par, alpha=alpha)
    Ki <- solve(K + diag(g, n))
    one.vec <- matrix(1, ncol = 1, nrow = n)
    mu.hat <- 0
    tau2hat <- drop(t(y - mu.hat) %*% Ki %*% (y - mu.hat) / nrow(X))

    return(list(K = K, Ki = Ki, X = X, y = y, t=t, theta_x=theta_x, theta_t=theta_t, beta=beta, delta=delta, alpha=alpha, g = g, mu.hat = mu.hat, tau2hat = tau2hat, constant = constant))
  }
}

pred.GP.nonsep <- function(fit, xnew, tnew) {
  xnew <- as.matrix(xnew)

  Ki <- fit$Ki
  theta_x <- fit$theta_x
  theta_t <- fit$theta_t
  beta <- fit$beta
  delta <- fit$delta
  g <- fit$g
  X <- fit$X
  y <- fit$y
  t <- fit$t
  tau2hat <- fit$tau2hat
  mu.hat <- fit$mu.hat

  KXX <- cor.sep.sqex(xnew, t=tnew, param=c(theta_x,theta_t,beta,delta))
  KX <- cor.sep.sqex(xnew, X, t=tnew, tnew=t,param=c(theta_x,theta_t,beta,delta))

  mup2 <- mu.hat + KX %*% Ki %*% (y - mu.hat)
  Sigmap2 <- pmax(0, diag(tau2hat * (KXX - KX %*% Ki %*% t(KX))))

  return(list(mu = mup2, sig2 = Sigmap2))
}
