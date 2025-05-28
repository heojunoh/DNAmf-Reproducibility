h1 <- function(x, mu.hat, theta_x, X2, Telement, a, w1.x2, previous_mu, previous_sig2, d, delta, beta) {
  
  dist_matrix <- distance(t(t(x) / sqrt(theta_x[-(d + 1)])), t(t(X2) / sqrt(theta_x[-(d + 1)])))
  
  exp_dist <- exp(-sweep(dist_matrix, 2, Telement, `*`))
  
  sqrt_term <- sqrt(theta_x[d + 1] / (theta_x[d + 1] + 2 * previous_sig2))
  
  diff_mu_w1 <- drop(outer(previous_mu, w1.x2, FUN = "-"))
  denom <- theta_x[d + 1] + 2 * previous_sig2  
  exp_mu <- exp(-sweep((diff_mu_w1^2 / denom), 2, Telement, `*`))
  
  # Combine terms and compute sum over training points
  weights <- a * Telement^((d+1)/2 + delta/beta)
  pred <- mu.hat + (exp_dist * sqrt_term * exp_mu) %*% weights
  
  return(pred)
}

h2_single <- function(x, mu.hat, theta_x, X2, Telement, a, w1.x2, previous_mu, previous_sig2, 
                      d, delta, beta, tau2hat, Ci, predy, c_sum, diff_w1_sq, Telement_outer, power_term) {
  
  dist <- distance(t(x / sqrt(theta_x[-(d + 1)])), t(t(X2) / sqrt(theta_x[-(d + 1)])))
  v <- exp(-dist * t(Telement))
  mat1 <- drop(v %o% v)
  
  # Zeta components
  denom <- theta_x[d + 1] + 2 * c_sum * previous_sig2
  sqrt_term <- sqrt(theta_x[d + 1] / denom)
  
  term1 <- drop(outer(Telement * (w1.x2 - previous_mu)^2, Telement * (w1.x2 - previous_mu)^2, "+"))
  term2 <- (2 / theta_x[d + 1]) * Telement_outer * previous_sig2 * diff_w1_sq
  exp_term <- exp(-(term1 + term2) / denom)
  
  mat <- mat1 * sqrt_term * exp_term * power_term
  
  # Compute h2
  quad <- drop(t(a) %*% mat %*% a)
  trace <- sum(Ci * mat)
  h2_val <- tau2hat - (predy - mu.hat)^2 + quad - tau2hat * trace
  return(pmax(0, h2_val))
}

h2 <- function(x, mu.hat, theta_x, X2, Telement, a, w1.x2, previous_mu, previous_sig2, d, delta, beta, tau2hat, Ci, predy) {
  # Precompute terms outside the loop
  c_sum <- outer(c(Telement), c(Telement), "+")
  diff_w1_sq <- drop(outer(w1.x2, w1.x2, "-"))^2
  Telement_outer <- Telement %*% t(Telement)
  power_term <- (Telement_outer)^((d+1)/2 + delta/beta)
  
  # Apply h2_single to each test point
  predsig2 <- numeric(nrow(x))
  for (i in 1:nrow(x)) {
    predsig2[i] <- h2_single(x[i, ], mu.hat, theta_x, X2, Telement, a, w1.x2, previous_mu[i], previous_sig2[i], 
                             d, delta, beta, tau2hat, Ci, predy[i], c_sum, diff_w1_sq, Telement_outer, power_term)
  }
  return(predsig2)
}

closed_form <- function(fit1, fit2, targett, nn, tt, TT, level, x, XX=NULL, pseudo_yy=NULL, ...){
  
  theta_x <- fit2$theta_x
  theta_t <- fit2$theta_t
  beta <- fit2$beta
  alpha <- fit2$alpha
  delta <- fit2$delta
  g <- fit2$g
  
  if(is.null(pseudo_yy)){
    ### prediction ###
    y <- fit2$y
    n <- length(y)
    Ci <- fit2$Ki
    tau2hat <- fit2$tau2hat
    mu.hat <- fit2$mu.hat
    X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
    w1.x2 <- fit2$X[, (d + 1)]
  }else{
    y <- do.call(rbind, pseudo_yy$y_star[-1])
    n <- length(y)
    
    # update K based on the new pseudo_yy
    L <- length(XX$X_star)
    X_m1 <- do.call(rbind, XX$X_star[-1])
    Y_mL <- do.call(rbind, lapply(2:L, function(l) {
      pseudo_yy$y_star[[l-1]][checkindices(XX$X_star[[l-1]], XX$X_star[[l]]), , drop = FALSE]
    }))
    t_m1 <- unlist(lapply(2:L, function(l) rep(tt[l], nrow(XX$X_star[[l]]))))
    
    K <- cor.sep.sqex(cbind(X_m1,Y_mL), t=t_m1, param=c(theta_x,theta_t,beta,delta), alpha=alpha)
    Ci <- solve(K + diag(g, n))
    one.vec <- matrix(1, ncol = 1, nrow = n)
    mu.hat <- drop((t(one.vec) %*% Ci %*% y) / (t(one.vec) %*% Ci %*% one.vec))
    tau2hat <- drop(t(y - mu.hat) %*% Ci %*% (y - mu.hat) / n)
    X2 <- X_m1
    w1.x2 <- c(Y_mL)
  }
  a <- Ci %*% (y - mu.hat)
  d <- ncol(fit1$X)
  x <- matrix(x, ncol = d)
  pred.fit <- pred.GP(fit1, x)
  
  ### calculate the closed form ###
  
  if (level > 2){
    pred_result <- list()
    pred_result[["mu_1"]] <- pred.fit$mu
    pred_result[["sig2_1"]] <- pred.fit$sig2
    
    for (j in 2:level) { ### from level 2 to level of current design
      Telement <- 1/(outer(TT, tt[j], "-")^(2*alpha)/theta_t + 1)^beta # c_i = c(t,t_i)
      # mean & var
      predy_temp <- h1(x, mu.hat, theta_x, X2, Telement, a, w1.x2, pred_result[[paste0("mu_", j-1)]], pred_result[[paste0("sig2_", j-1)]], d, delta, beta)
      predsig2_temp <- h2(x, mu.hat, theta_x, X2, Telement, a, w1.x2, pred_result[[paste0("mu_", j-1)]], pred_result[[paste0("sig2_", j-1)]], d, delta, beta, tau2hat, Ci, predy_temp)
      
      pred_result[[paste0("mu_", j)]] <- predy_temp
      pred_result[[paste0("sig2_", j)]] <- predsig2_temp
    }
    
    ### for target mesh size ###
    if(any(tt == targett)){ # if target mesh belongs to initial design, use previous outputs with that specific mesh size (for interpolation)
      pred_result[["mu"]] <- pred_result[[paste0("mu_", which(tt == targett))]]
      pred_result[["sig2"]] <- pred_result[[paste0("sig2_", which(tt == targett))]]
    }else{
      # find the smallest mesh which is larger than target t
      tt_place <- findInterval(-targett, -tt) 
      Telement <- 1/(outer(TT, targett, "-")^(2*alpha)/theta_t + 1)^beta # c_i = c(t,t_i)
      # mean & var
      predy <- h1(x, mu.hat, theta_x, X2, Telement, a, w1.x2, pred_result[[paste0("mu_", tt_place)]], pred_result[[paste0("sig2_", tt_place)]], d, delta, beta)
      predsig2 <- h2(x, mu.hat, theta_x, X2, Telement, a, w1.x2, pred_result[[paste0("mu_", tt_place)]], pred_result[[paste0("sig2_", tt_place)]], d, delta, beta, tau2hat, Ci, predy)
      
      pred_result[["mu"]] <- predy
      pred_result[["sig2"]] <- predsig2
    }
  } else {
    stop("Level should be larger than 2")
  }
  return(pred_result)
}

predict.DNAmf <- function(object, x, targett=0, nimpute=50,  ...) {
  t1 <- proc.time()
  
  ### check the object ###
  if (!inherits(object, "DNAmf")) {
    stop("The object is not of class \"DNAmf\" \n")
  }
  
  ### prediction ###
  nn <- object$nn
  tt <- object$t 
  TT <- object$TT 
  XX <- object$XX 
  yy <- object$yy 
  level  <- object$level
  
  fit1 <- object$fit1
  fit2 <- object$fit2
  
  if(object$nested){
    fit2 <- object$fit2
    pred_result <- closed_form(fit1, fit2, targett, nn, tt, TT, level, x)
  }else{
    pred1 <- object$pred1
    
    sum_mu <- sum_mu2_plus_sig2 <- sig2_star <- vector("list", level+1)
    names(sum_mu) <- mu_names <- c(paste0("mu_", 1:level), "mu")
    names(sum_mu2_plus_sig2) <- names(sig2_star) <- sig2_names <- c(paste0("sig2_", 1:level), "sig2")
    
    n <- nrow(matrix(x, ncol = ncol(fit1$X)))
    sum_mu <- setNames(rep(list(rep(0, n)), level + 1), mu_names)
    sum_mu2_plus_sig2 <- setNames(rep(list(rep(0, n)), level + 1), sig2_names)
    
    for (m in 1:nimpute) {
      # Generate imputed dataset for the m-th imputation
      yy <- imputer(XX, yy, tt, pred1, fit2)
      pred <- closed_form(fit1, fit2, targett, nn, tt, TT, level, x, XX, yy)
      
      sum_mu <- mapply(`+`, sum_mu, pred[mu_names], SIMPLIFY = FALSE)
      sum_mu2_plus_sig2 <- mapply(function(old, mu_vec, sig2_vec) old + mu_vec^2 + sig2_vec,
                                  sum_mu2_plus_sig2, pred[mu_names], pred[sig2_names], SIMPLIFY = FALSE)
    }
    
    # Compute final mu* and sig2* for each component
    mu_star <- lapply(sum_mu, function(s) s / nimpute)
    sig2_star <- mapply(function(mu2_plus_sig2, mu_star_vec) (mu2_plus_sig2 / nimpute) - (mu_star_vec)^2, 
                        sum_mu2_plus_sig2, mu_star, SIMPLIFY = FALSE, USE.NAMES = TRUE)
    
    # Combine results into pred_result
    pred_result <- c(mu_star, sig2_star)
  }
  
  pred_result[["time"]] <- (proc.time() - t1)[3]
  
  return(pred_result)
}

