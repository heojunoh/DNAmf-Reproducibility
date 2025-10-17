### Figure 3 ###
# Define a reusable plotting function
create_plot <- function(i, mesh_size, x, pred_mu, pred_sig2, X_points = NULL, y_points = NULL, add_points = TRUE, yylim) {
  data_long <- data.frame(
    x = x,
    mu = pred_mu,
    lower = pred_mu - qnorm(0.995) * sqrt(pred_sig2),
    upper = pred_mu + qnorm(0.995) * sqrt(pred_sig2),
    group = factor(i)
  )
  
  facet_labels <- setNames(paste0("Mesh size = ", round(mesh_size, digits=3)), as.character(i))
  
  p <- ggplot() +
    geom_line(data = data_long, aes(x = x, y = mu), linewidth = 1, color = scales::hue_pal()(6)[5]) +
    geom_ribbon(data = data_long, aes(x = x, ymin = lower, ymax = upper), alpha = 0.2, fill = scales::hue_pal()(6)[5]) +
    facet_wrap(~group, scales = "free_y", nrow = 1, labeller = as_labeller(facet_labels), strip.position = "top") +
    ylim(-yylim, yylim) + theme_minimal() + theme_bw() +
    theme(legend.position = "none",
          strip.text = element_text(size = 12),
          strip.text.x = element_text(size = 12, margin = margin(t = 5, b = 5)),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank()) +
    labs(y = NULL) +
    geom_function(fun = function(x) fl(x, mesh_size),
                  inherit.aes = FALSE,
                  color = "black", linetype = "dashed", linewidth = 1,
                  data = data.frame(group = factor(i)))
  
  if (add_points && !is.null(X_points) && !is.null(y_points)) {
    data_points <- data.frame(x = X_points, y = y_points, group = factor(i))
    p <- p + geom_point(data = data_points, aes(x = x, y = y), color = "red", shape = 16, size = 2.1)
  }
  
  return(p)
}

### Additive function ###
fl <- function(x, t){
  term1 <- exp(-1.4 * x) * cos(3.5 * pi * x)
  term2 <- t^2 * (sin(40 * x))/10
  term1 + term2
}
### training data ###
# n1 <- 13; n2 <- 10; n3 <- 7; n4 <- 4; n5 <- 1; 
# m1 <- 3.0; m2 <- 2.5; m3 <- 2.0; m4 <- 1.5; m5 <- 1.0;

c <- 0.7; bet <- -2
n5 <- 1; n4 <- round(c^bet*n5); n3 <- round(c^(2*bet)*n5); n2 <- round(c^(3*bet)*n5); n1 <- round(c^(4*bet)*n5);

m1 <- 2.5; m2 <- m1*c; m3 <- m2*c; m4 <- m3*c; m5 <- m4*c;
d <- 1
eps <- sqrt(.Machine$double.eps)
x <- seq(0,1,0.01)

set.seed(3)
NestDesign <- NestedX(c(n1,n2,n3,n4,n5),d)
X1 <- NestDesign[[1]]
X2 <- NestDesign[[2]]
X3 <- NestDesign[[3]]
X4 <- NestDesign[[4]]
X5 <- NestDesign[[5]]
y1 <- fl(X1, t=m1)
y2 <- fl(X2, t=m2)
y3 <- fl(X3, t=m3)
y4 <- fl(X4, t=m4)
y5 <- fl(X5, t=m5)


### DNAmf ###
fit.DNAmf.sqex <- DNAmf.sqex(X1, y1,
                             X=rbind(X2, X3, X4, X5),
                             y=rbind(y2, y3, y4, y5),
                             nn=c(length(y1),length(y2),length(y3),length(y4),length(y5)),
                             t=c(m1,m2,m3,m4,m5), multi.start=10, 
                             constant=TRUE)
pred.DNAmf <- predict.DNAmf(fit.DNAmf.sqex, x, targett=0)
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
  create_plot(i, m, x, mu, sig2, X, y, add_points = !is.null(X), yylim=3)
}, i = 1:6, m = mesh_sizes, mu = mu_list, sig2 = sig2_list, X = X_list, y = y_list, SIMPLIFY = FALSE)

# Arrange plots
firstrow <- annotate_figure(ggarrange(plotlist = plots, ncol = 6, nrow = 1) + ylab("y"),
                            left = text_grob("y", rot = 90, size = 10))


### Non-additive function ###
fl <- function(x, t){
  term1 <- sin(10 * pi * x / (5+t))
  term2 <- 0.2 * sin(8 * pi * x)
  term1 + term2
}
### training data ###
# m1 <- 2.5; m2 <- 2.0; m3 <- 1.5; m4 <- 1.0; m5 <- 0.5;

c <- 0.7; bet <- -2
n5 <- 1; n4 <- round(c^bet*n5); n3 <- round(c^(2*bet)*n5); n2 <- round(c^(3*bet)*n5); n1 <- round(c^(4*bet)*n5);

m1 <- 2.5; m2 <- m1*c; m3 <- m2*c; m4 <- m3*c; m5 <- m4*c;

set.seed(81)
NestDesign <- NestedX(c(n1,n2,n3,n4,n5),d)
X1 <- NestDesign[[1]]
X2 <- NestDesign[[2]]
X3 <- NestDesign[[3]]
X4 <- NestDesign[[4]]
X5 <- NestDesign[[5]]
y1 <- fl(X1, t=m1)
y2 <- fl(X2, t=m2)
y3 <- fl(X3, t=m3)
y4 <- fl(X4, t=m4)
y5 <- fl(X5, t=m5)


### DNAmf ###
fit.DNAmf.sqex <- DNAmf.sqex(X1, y1,
                             X=rbind(X2, X3, X4, X5),
                             y=rbind(y2, y3, y4, y5),
                             nn=c(length(y1),length(y2),length(y3),length(y4),length(y5)),
                             t=c(m1,m2,m3,m4,m5), multi.start=10, 
                             constant=TRUE)
pred.DNAmf <- predict.DNAmf(fit.DNAmf.sqex, x, targett=0)
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
  create_plot(i, m, x, mu, sig2, X, y, add_points = !is.null(X), yylim=1.5)
}, i = 1:6, m = mesh_sizes, mu = mu_list, sig2 = sig2_list, X = X_list, y = y_list, SIMPLIFY = FALSE)

# Arrange plots
secondrow <- annotate_figure(ggarrange(plotlist = plots, ncol = 6, nrow = 1) + ylab("y"),
                             left = text_grob("y", rot = 90, size = 10))

# Combine with firstrow
figure3 <- ggarrange(firstrow, secondrow, nrow = 2)
