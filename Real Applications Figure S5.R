### Figure S5 ###
# Define a reusable plotting function
create_plot <- function(i, mesh_size, x, pred_mu, pred_sig2, X_points = NULL, y_points = NULL, add_points = TRUE) {
  data_long <- data.frame(
    x = x.test,
    mu = pred_mu,
    lower = pred_mu - qnorm(0.995) * sqrt(pred_sig2),
    upper = pred_mu + qnorm(0.995) * sqrt(pred_sig2),
    group = factor(i)
  )
  
  facet_labels <- setNames(paste0("Mesh size = ", mesh_size), as.character(i))
  
  p <- ggplot() +
    geom_line(data = data_long, aes(x = x, y = mu), size = 1, color = scales::hue_pal()(6)[5]) +
    geom_ribbon(data = data_long, aes(x = x, ymin = lower, ymax = upper), alpha = 0.2, fill = scales::hue_pal()(6)[5]) +
    facet_wrap(~group, scales = "free_y", nrow = 1, labeller = as_labeller(facet_labels), strip.position = "top") +
    ylim(0.6, 1.75) + theme_minimal() + theme_bw() +
    theme(legend.position = "none",
          strip.text = element_text(size = 12),
          strip.text.x = element_text(size = 12, margin = margin(t = 5, b = 5)),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank()) +
    labs(y = NULL)
  
  if (mesh_size == 0) {
    true_data <- data.frame(x = x.test, y = y.test, group = factor(i))
    p <- p + geom_line(data = true_data, aes(x = x, y = y),
                       inherit.aes = FALSE,
                       color = "black", linetype = "dashed", size = 1)
  }
  
  if (add_points && !is.null(X_points) && !is.null(y_points)) {
    data_points <- data.frame(x = X_points, y = y_points, group = factor(i))
    p <- p + geom_point(data = data_points, aes(x = x, y = y), color = "red", shape = 16, size = 2.1)
  }
  
  return(p)
}

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
  write.csv(xx, "Poisson-equation/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "Poisson-equation/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("Poisson-equation/SolvePoissonEqn.m", verbose = FALSE, desktop = FALSE,
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("Poisson-equation/generate_text/temp_to_r.txt", sep = ",")
  out <- matrix(as.numeric(unlist(d2)), nrow = n)
  return(matrix(out[,3], ncol=1)) # row-wise output
}

### training data ###
n1 <- 15; n2 <- 12; n3 <- 9; n4 <- 6; n5 <- 3
m1 <- 0.09; m2 <- 0.07; m3 <- 0.05; m4 <- 0.03; m5 <- 0.01; targett=0
d <-1
eps <- sqrt(.Machine$double.eps)
### test data ###
x.test <- matrix(seq(0,1,length=100),ncol=1)
y.test <- matrix(sapply(x.test,f_true_0),ncol=1)

set.seed(2)
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


### DNAmf ###
fit.DNAmf.sqex <- DNAmf.sqex(X1, y1,
                             X=rbind(X2, X3, X4, X5),
                             y=rbind(y2, y3, y4, y5),
                             nn=c(length(y1),length(y2),length(y3),length(y4),length(y5)),
                             t=c(m1,m2,m3,m4,m5), multi.start=50, constant=TRUE)
pred.DNAmf <- predict.DNAmf(fit.DNAmf.sqex, x.test, targett=targett)
predydiffu.sqex <- pred.DNAmf$mu
predsig2diffu.sqex <- pred.DNAmf$sig2

### Plotting ###
mesh_sizes <- c(m1, m2, m3, m4, m5, 0)
mu_list <- list(pred.DNAmf$mu_1, pred.DNAmf$mu_2, pred.DNAmf$mu_3,
                pred.DNAmf$mu_4, pred.DNAmf$mu_5, pred.DNAmf$mu)
sig2_list <- list(pred.DNAmf$sig2_1, pred.DNAmf$sig2_2, pred.DNAmf$sig2_3,
                  pred.DNAmf$sig2_4, pred.DNAmf$sig2_5, pred.DNAmf$sig2)
X_list <- list(X1, X2, X3, X4, X5, NULL)
y_list <- list(y1, y2, y3, y4, y5, NULL)

plots <- mapply(function(i, m, mu, sig2, X, y) {
  create_plot(i, m, x, mu, sig2, X, y, add_points = !is.null(X))
}, i = 1:6, m = mesh_sizes, mu = mu_list, sig2 = sig2_list, X = X_list, y = y_list, SIMPLIFY = FALSE)

# Arrange plots
poissonrow <- annotate_figure(ggarrange(plotlist = plots, ncol = 6, nrow = 1) + ylab("y"),
                            left = text_grob("y", rot = 90, size = 10))

# Define a reusable plotting function
create_plot <- function(i, mesh_size, x, pred_mu, pred_sig2, X_points = NULL, y_points = NULL, add_points = TRUE) {
  data_long <- data.frame(
    x = x,
    mu = pred_mu,
    lower = pred_mu - qnorm(0.995) * sqrt(pred_sig2),
    upper = pred_mu + qnorm(0.995) * sqrt(pred_sig2),
    group = factor(i)
  )
  
  facet_labels <- setNames(paste0("Mesh size = ", mesh_size), as.character(i))
  
  p <- ggplot() +
    geom_line(data = data_long, aes(x = x, y = mu), size = 1, color = scales::hue_pal()(6)[5]) +
    geom_ribbon(data = data_long, aes(x = x, ymin = lower, ymax = upper), alpha = 0.2, fill = scales::hue_pal()(6)[5]) +
    facet_wrap(~group, scales = "free_y", nrow = 1, labeller = as_labeller(facet_labels), strip.position = "top") +
    ylim(-0, 1) + theme_minimal() + theme_bw() +
    theme(legend.position = "none",
          strip.text = element_text(size = 12),
          strip.text.x = element_text(size = 12, margin = margin(t = 5, b = 5)),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank()) +
    labs(y = NULL) 
  
  if (mesh_size == 0) {
    p <- p + geom_function(fun = function(x) heat(x, mesh_size),
                           inherit.aes = FALSE,
                           color = "black", linetype = "dashed", size = 1,
                           data = data.frame(group = factor(i)))
  }
  
  if (add_points && !is.null(X_points) && !is.null(y_points)) {
    data_points <- data.frame(x = X_points, y = y_points, group = factor(i))
    p <- p + geom_point(data = data_points, aes(x = x, y = y), color = "red", shape = 16, size = 2.1)
  }
  
  return(p)
}

### Time dependent PDE; Heat equation; smaller t corresponds to longer time point (0.5-t)  ###
library(ReacTran)
library(deSolve)
heat <- function(X, t, L = 1.0, N = 1000, D = 0.1) {
  if (!is.matrix(X)) {
    if (!is.numeric(X) || any(X < 0) || any(X > L))
      stop("'X'must be numeric (vector or matrix) in [0, ", L, "]")
    Xmat <- matrix(X, ncol = 1)
  } else {
    if (!is.numeric(X) || any(X < 0) || any(X > L))
      stop("All entries of 'X' must lie in [0, ", L, "]")
    Xmat <- X
  }
  
  if (length(t) != 1 || !is.numeric(t) || t < 0)
    stop("'t' must be a single numeric >= 0")
  
  ### grid & initial condition ###
  grid <- setup.grid.1D(x.up = 0, x.down = L, N = N)
  u_initial <- exp(-100 * (grid$x.mid - L/2)^2)
  
  ### define model ###
  heat_model <- function(time, u, parms) {
    flux <- tran.1D(C = u, D = D, dx = grid)
    list(dC = flux$dC)
  }
  
  ### solve equation at times 0 and t ###
  sol <- ode.1D(y = u_initial, times = c(0, 0.5-t), func = heat_model,
                parms = NULL, nspec = 1, dimens = N)
  u_t <- sol[2, -1] # row 2; time = t
  
  ### interpolation (clamp at endpoints) ###
  x_vec <- as.vector(Xmat)
  u_vec <- approx(x = grid$x.mid, y = u_t, xout = x_vec, rule = 2)$y
  
  ### return ###
  matrix(u_vec, nrow = nrow(Xmat), ncol = ncol(Xmat), dimnames = dimnames(Xmat))
}

### training data ###
n1 <- 15; n2 <- 12; n3 <- 9; n4 <- 6; n5 <- 3;
m1 <- 0.5; m2 <- 0.4; m3 <- 0.3; m4 <- 0.2; m5 <- 0.1; targett=0
d <- 1
eps <- sqrt(.Machine$double.eps)
x <- seq(0,1,0.01)

set.seed(2)
NestDesign <- NestedX(c(n1,n2,n3,n4,n5),d)
X1 <- NestDesign[[1]]
X2 <- NestDesign[[2]]
X3 <- NestDesign[[3]]
X4 <- NestDesign[[4]]
X5 <- NestDesign[[5]]
y1 <- heat(X1,m1)
y2 <- heat(X2,m2)
y3 <- heat(X3,m3)
y4 <- heat(X4,m4)
y5 <- heat(X5,m5)


### DNAmf ###
fit.DNAmf.sqex <- DNAmf.sqex(X1, y1,
                             X=rbind(X2, X3, X4, X5),
                             y=rbind(y2, y3, y4, y5),
                             nn=c(length(y1),length(y2),length(y3),length(y4),length(y5)),
                             t=c(m1,m2,m3,m4,m5), multi.start=10, constant=TRUE)
pred.DNAmf <- predict.DNAmf(fit.DNAmf.sqex, x, targett=targett)
predydiffu.sqex <- pred.DNAmf$mu
predsig2diffu.sqex <- pred.DNAmf$sig2

### Plotting ###
mesh_sizes <- c(m1, m2, m3, m4, m5, 0)
mu_list <- list(pred.DNAmf$mu_1, pred.DNAmf$mu_2, pred.DNAmf$mu_3,
                pred.DNAmf$mu_4, pred.DNAmf$mu_5, pred.DNAmf$mu)
sig2_list <- list(pred.DNAmf$sig2_1, pred.DNAmf$sig2_2, pred.DNAmf$sig2_3,
                  pred.DNAmf$sig2_4, pred.DNAmf$sig2_5, pred.DNAmf$sig2)
X_list <- list(X1, X2, X3, X4, X5, NULL)
y_list <- list(y1, y2, y3, y4, y5, NULL)

plots <- mapply(function(i, m, mu, sig2, X, y) {
  create_plot(i, m, x, mu, sig2, X, y, add_points = !is.null(X))
}, i = 1:6, m = mesh_sizes, mu = mu_list, sig2 = sig2_list, X = X_list, y = y_list, SIMPLIFY = FALSE)

# Arrange plots
heatrow <- annotate_figure(ggarrange(plotlist = plots, ncol = 6, nrow = 1) + ylab("y"),
                           left = text_grob("y", rot = 90, size = 10))

# Combine two plots
figureS5 <- ggarrange(poissonrow, heatrow, nrow = 2)

