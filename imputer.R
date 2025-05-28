
### Construct nested X^* ###
makenested <- function(X_list){
  # Initialize the list to store nested matrices
  L <- length(X_list)
  X_star <- vector("list",L)
  
  # Set X^*_L = XL
  X_star[[L]] <- X_list[[L]]
  
  # Recursively build X^*_l for l = L down to 1
  for (l in (L-1):1) {
    X_star[[l]] <- unique(rbind(X_list[[l]], X_star[[l+1]])) # intersections are on X_list
  }
  
  # Compute X_tilde where X_tilde[[l]] = X_star[[l]] \ X_list[[l]]
  X_tilde <- vector("list", L)
  for (l in 1:(L-1)) {
    combined <- rbind(X_list[[l]], X_star[[l]])
    dup <- duplicated(combined)
    indices <- (nrow(X_list[[l]]) + 1) : nrow(combined)
    X_tilde[[l]] <- X_star[[l]][!dup[indices], , drop = FALSE]
  }
  
  return(list(X_star=X_star, X_list=X_list, X_tilde = X_tilde))
}

### Find the indices in nested design ###
checkindices <- function(XX1, XX2) {
  dist_matrix <- fields::rdist(XX2, XX1)
  matches <- which(dist_matrix < sqrt(.Machine$double.eps), arr.ind = TRUE)
  return(matches[, 2])
}

### Draw y_2 tilde - y_{L-1} tilde by ESS ###
imputer <- function(XX, yy, # list; XX(Xstar, X, Xtilde), yy(ystar, y, ytilde)
                  t, pred1, fit2, nstep=1){ 
  
  L <- length(XX$X_star)
  X <- fit2$X
  y <- fit2$y
  K <- fit2$K
  g <- fit2$g
  alpha <- fit2$alpha
  K <- K + diag(g, length(y))
  
  n_list <- sapply(yy$y_list, length) 
  n_star <- sapply(yy$y_star, length) 
  ncum_star <- c(0, cumsum(n_star[-1]))   # cum sum of n_star
  n_tilde <- sapply(yy$y_tilde, length) 
  ncum_tilde <- c(0,cumsum(n_tilde[2:(L-1)])) # cum sum of n_tilde
  
  X_m1 <- do.call(rbind, XX$X_star[-1])
  Y_m1 <- do.call(rbind, yy$y_star[-1])
  t_m1 <- unlist(lapply(2:L, function(l) rep(t[l], nrow(XX$X_star[[l]]))))
  params <- c(fit2$theta_x, fit2$theta_t, fit2$beta, fit2$delta)
  mu.hat <- fit2$mu.hat
  tau2hat <- fit2$tau2hat
  
  # ESS loop
  for (i in 1:nstep) {
    
    # Draw from prior distribution
    yy$y_tilde[[1]] <- t(mvtnorm::rmvnorm(1, mean = pred1$mu, sigma = pred1$cov))
    yy$y_star[[1]] <- rbind(yy$y_list[[1]],yy$y_tilde[[1]])
    
    ### sampling from Y_tilde give Y_list
    list_idx <- unlist(mapply(seq, ncum_star[1:(L-1)]+1, ncum_star[1:(L-1)]+n_list[2:L]))
    K_list <- K[list_idx, list_idx]     # K corresponding to X_list
    Ki_list <- solve(K_list)
    K_tilde <- K[-list_idx, -list_idx]  # K corresponding to X_tilde
    K_list_tilde <- K[list_idx, -list_idx] # K corresponding to X_list given X_tilde
    y_list <- y[list_idx]
    
    cond_mean <- mu.hat + t(K_list_tilde) %*% Ki_list %*% (y_list - mu.hat)
    cond_var <- tau2hat * (K_tilde - t(K_list_tilde) %*% Ki_list %*% K_list_tilde)
    cond_var <- (cond_var+t(cond_var))/2
    
    y_prior <- t(mvtnorm::rmvnorm(1, mean = cond_mean, sigma = cond_var))
    
    for(l in 2:(L-1)){
      yy$y_tilde[[l]] <- matrix(y_prior[(ncum_tilde[l-1]+1):ncum_tilde[l]],ncol=1)
      yy$y_star[[l]] <- rbind(yy$y_list[[l]], yy$y_tilde[[l]])
    }
    
    Y_m1 <- do.call(rbind, yy$y_star[-1])
    Y_mL <- do.call(rbind, lapply(2:L, function(l) {
      yy$y_star[[l-1]][checkindices(XX$X_star[[l-1]], XX$X_star[[l]]), , drop = FALSE]
    }))
    y <- Y_m1
    X[,ncol(X)] <- Y_mL
    K <- cor.sep.sqex(X, t=t_m1, param=params, alpha=alpha) + diag(g, length(y))
  } # End of ESS loop
  return(yy)
}

