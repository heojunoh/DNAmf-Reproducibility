### Non-Additive function ###
fl <- function(x, t){
  term1 <- sin(10 * pi * x / (5+t))
  term2 <- 0.2 * sin(8 * pi * x)
  term1 + term2
}
### training data ###
# n1 <- 18; n2 <- 14; n3 <- 11; n4 <- 8; n5 <- 5; 
n1 <- 13; n2 <- 10; n3 <- 7; n4 <- 4; n5 <- 1;
m1 <- 2.5; m2 <- 2; m3 <- 1.5; m4 <- 1; m5 <- 0.5; 

set.seed(4) 
N <- n1+n2+n3+n4+n5
X_star <- matrix(seq(0,1,length=N),ncol=1)
idx <- sample(1:N, n1, replace = FALSE)
X1 <- X_star[idx,,drop=FALSE]
X_star <- X_star[-idx,,drop=FALSE]; N <- N - n1
idx <- sample(1:N, n2, replace = FALSE)
X2 <- X_star[idx,,drop=FALSE] 
X_star <- X_star[-idx,,drop=FALSE]; N <- N - n2
idx <- sample(1:N, n3, replace = FALSE)
X3 <- X_star[idx,,drop=FALSE] 
X_star <- X_star[-idx,,drop=FALSE]; N <- N - n3
idx <- sample(1:N, n4, replace = FALSE)
X4 <- X_star[idx,,drop=FALSE] 
X_star <- X_star[-idx,,drop=FALSE]; N <- N - n4
idx <- sample(1:N, n5, replace = FALSE)
X5 <- X_star[idx,,drop=FALSE] 
y1 <- fl(X1, t=m1)
y2 <- fl(X2, t=m2)
y3 <- fl(X3, t=m3)
y4 <- fl(X4, t=m4)
y5 <- fl(X5, t=m5)

fit.DNAmf.sqex <- DNAmf.sqex(X1, y1, X=rbind(X2, X3, X4, X5), y=rbind(y2, y3, y4, y5),
                             nn=c(length(y1),length(y2),length(y3),length(y4),length(y5)),
                             t=c(m1,m2,m3,m4,m5), multi.start=10, n.iter=100, constant=TRUE)
X_list_2 <- fit.DNAmf.sqex$XX$X_tilde
y_list_2 <- fit.DNAmf.sqex$yy$y_tilde


fit.DNAmf.sqex <- DNAmf.sqex(X1, y1, X=rbind(X2, X3, X4, X5), y=rbind(y2, y3, y4, y5),
                             nn=c(length(y1),length(y2),length(y3),length(y4),length(y5)),
                             t=c(m1,m2,m3,m4,m5), multi.start=10, n.iter=100, constant=TRUE)
X_list_3 <- fit.DNAmf.sqex$XX$X_tilde
y_list_3 <- fit.DNAmf.sqex$yy$y_tilde


fit.DNAmf.sqex <- DNAmf.sqex(X1, y1, X=rbind(X2, X3, X4, X5), y=rbind(y2, y3, y4, y5),
                             nn=c(length(y1),length(y2),length(y3),length(y4),length(y5)),
                             t=c(m1,m2,m3,m4,m5), multi.start=10, n.iter=100, constant=TRUE)
X_list_5 <- fit.DNAmf.sqex$XX$X_tilde
y_list_5 <- fit.DNAmf.sqex$yy$y_tilde

### --- Build df_points (Original + Pseudo outputs) --- ###
make_pseudo_df <- function(X_list, y_list, label) {
  df_list <- lapply(1:4, function(i) {
    xi <- X_list[[i]]
    yi <- y_list[[i]]
    if (length(xi) > 0 && length(yi) > 0 && length(xi) == length(yi)) {
      data.frame(
        x = xi,
        y = yi,
        level = factor(i),
        type = label
      )
    } else {
      NULL
    }
  })
  do.call(rbind, df_list)
}

df_points_Original <- data.frame(
  x = c(X1, X2, X3, X4),
  y = c(y1, y2, y3, y4),
  level = factor(rep(1:4, times = c(length(X1), length(X2), length(X3), length(X4)))),
  type = "Original"
)

df_points_pseudo1 <- make_pseudo_df(X_list_2, y_list_2, "Iteration 1")
df_points_pseudo2 <- make_pseudo_df(X_list_3, y_list_3, "Iteration 2")
df_points_pseudo3 <- make_pseudo_df(X_list_5, y_list_5, "Iteration 3")

df_points <- rbind(df_points_Original, df_points_pseudo1, df_points_pseudo2, df_points_pseudo3)
df_points$type <- factor(df_points$type, levels = c("Original", "Iteration 1", "Iteration 2", "Iteration 3"))

### --- Build curve_df --- ###
curve_df <- data.frame(
  x = rep(x, 4),
  y = c(fl(x, t=m1),fl(x, t=m2),fl(x, t=m3),fl(x, t=m4)),
  level = factor(rep(1:4, each = length(x)))
)

### --- Create mesh_labels --- ###
mesh_labels <- c(
  "1" = paste0("Mesh size = ", m1),
  "2" = paste0("Mesh size = ", m2),
  "3" = paste0("Mesh size = ", m3),
  "4" = paste0("Mesh size = ", m4)
)

### --- Plotting --- ###
levels <- levels(df_points$level)

make_panel <- function(i, show_legend = FALSE) {
  pts   <- df_points[df_points$level == i, ]
  curve <- curve_df[curve_df$level == i, ]
  
  ggplot() +
    geom_point(data = pts, aes(x = x, y = y, color = type, shape = type), size = 3, stroke = 0.8) +
    geom_line(data = curve, aes(x = x, y = y), size = 1, linetype = "dashed", color = "black") +
    scale_color_manual(values = c(
      "Original" = "red",
      "Iteration 1" = "green",
      "Iteration 2" = "blue",
      "Iteration 3" = "purple"
    )) +
    scale_shape_manual(values = c(
      "Original" = 16,
      "Iteration 1" = 0,
      "Iteration 2" = 2,
      "Iteration 3" = 5
    )) +
    facet_wrap(
      ~level, nrow = 1, scales = "free_y",
      strip.position = "top",
      labeller = as_labeller(mesh_labels)
    ) +
    coord_cartesian(ylim = c(-1.86, 1.86)) +
    labs(x = NULL, y = NULL, color = NULL, shape = NULL) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_text(size = 12),
      legend.position = if (show_legend) "bottom" else "none"
    )
}

# Extract one panel with the legend
plot_with_legend <- make_panel(levels[1], show_legend = TRUE)
legend <- get_legend(plot_with_legend)

# Create all panels
panels <- lapply(levels, function(lvl) make_panel(lvl, show_legend = FALSE))

# Combine panels
combined_plot <- ggarrange(
  plotlist = panels,
  nrow = 1,
  ncol = length(panels),
  widths = rep(1, length(panels))
)

# Final assembled plot with y-axis label and legend
figure7 <- ggarrange(
  annotate_figure(combined_plot,
                  left = text_grob("y", size = 10, rot = 90),
                  bottom = text_grob("x", size = 10)),
  legend,
  ncol = 1,
  heights = c(1, 0.1)
)

