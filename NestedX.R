subsetX <- function(X1 = NULL, X2 = NULL) {
  d <- dim(as.matrix(X2))[2] # d2
  n2 <- dim(as.matrix(X2))[1] # n2
  n1 <- dim(as.matrix(X1))[1] # n1

  dist <- 0
  for (i in 1:d) {
    grid <- expand.grid(X2[, i], X1[, i])
    dist <- dist + (grid[, 1] - grid[, 2])^2
  }
  dist.mat <- matrix(dist, n2, n1)
  indice <- max.col(-(dist.mat)) # find the minimum distance of column at each row

  X1 <- matrix(X1[-indice, ], ncol=d)
  X1 <- rbind(X1, matrix(X2, ncol=d))

  return(list(X = X1, le = length(indice)))
}


NestedX <- function(n, d) { # n; vector, d; dim
  
  # Input validation
  if (is.unsorted(rev(n), strictly = TRUE)) {
    stop("The number of design at each level must be descending order \n")
  }
  if (!all(n > 0)) {
    stop("The number of design at each level must be positive \n")
  }
  if (length(d) != 1) {
    stop("The dimension of design at each level must be same \n")
  }
  
  level <- length(n)
  if (level < 2) {
    stop("The level of design should be larger than 1 \n")
  }
  
  # Generating initial designs
  X <- list()
  for (i in 1:level) {
    X[[i]] <- maximinLHS(n[i], d)
  }
  
  # Subsetting designs
  indices <- list()
  for (i in (level - 1):1) {
    SB <- subsetX(matrix(X[[i]], ncol=d), matrix(X[[i + 1]], ncol=d))
    X[[i]] <- SB$X
    n_rows <- nrow(SB$X)
    indices[[i]] <- seq(n_rows - SB$le + 1, n_rows, by = 1)
  }
  
  # Constructing nested designs
  X_list <- list()
  current_idx <- 1:nrow(X[[1]])
  for (k in 1:level) {
    if (k > 1) {
      current_idx <- current_idx[indices[[k-1]]]
    }
    X_list[[k]] <- matrix(X[[1]][current_idx, ], ncol = d)
  }
  
  return(X = X_list)
}
