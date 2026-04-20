library(ggplot2)

# Main function to calculate the budget-constrained design
calculate_new_design <- function(C_total = 500, t_1 = 1, t_L = 0.01, n_min = 2, 
                                 lambda_max = 1, eta = 4, c1 = 1, gamma = 1, d = 1) {
  
  # Step 1 Parameters
  a <- d / (eta + d)
  b <- 1 - a
  
  valid_designs <- list() # Storage for valid designs as we test L
  
  # Step 2 & 3: Iterate over candidate L >= 2 to find the maximum valid L
  L <- 2
  while (TRUE) {
    l_seq <- 1:L
    k_seq <- 1:L
    
    # Compute T for the current candidate L
    T_val <- (t_1 / t_L)^(1 / (L - 1))
    
    # Pre-calculate the denominator sum (over k=1 to L)
    denom_terms <- (lambda_max^(a * (L - k_seq + 1))) * ((t_1 / t_L)^((gamma * b * (k_seq - 1)) / (L - 1)))
    denom_sum <- sum(denom_terms)
    
    # Calculate the numerator (for each l=1 to L)
    num_terms <- (lambda_max^(a * (L - l_seq + 1))) * ((t_1 / t_L)^(-(gamma * a * (l_seq - 1)) / (L - 1)))
    
    # Calculate continuous sample sizes using Eq (10)
    n_l_exact <- ((C_total * (t_1^gamma)) / c1) * (num_terms / denom_sum)
    
    # Floor to integers to strictly ensure we don't exceed the total budget
    n_l <- floor(n_l_exact)
    
    # Step 3: Check minimum sample size requirement
    if (min(n_l) >= n_min) {
      # This L is valid! Calculate exact tuning parameters and costs for the record
      t_l <- t_1 * T_val^(-(l_seq - 1))
      C_l <- c1 * t_l^(-gamma)
      
      # Save this valid design
      valid_designs[[as.character(L)]] <- data.frame(
        l = l_seq,
        n_l = n_l,
        t_l = t_l,
        C_l = C_l,
        T_val = T_val
      )
      
      # Try the next deeper hierarchy
      L <- L + 1
    } else {
      # The sample sizes dropped below n_min. The ceiling is hit.
      break
    }
  }
  
  # Retrieve the largest L that passed the test
  optimal_L <- L - 1
  
  if (optimal_L < 2) {
    stop("Budget C_total is too small to even support L=2 with the given n_min.")
  }
  # Get the winning design
  final_design <- valid_designs[[as.character(optimal_L)]]
  
  # 1. Calculate the initial leftover budget
  C_spent <- sum(final_design$n_l * final_design$C_l)
  C_leftover <- C_total - C_spent
  
  # 2. Loop backward from highest fidelity (L) down to lowest (1)
  for (lvl in optimal_L:1) {
    # Check if we can afford at least one sample at this level
    if (C_leftover >= final_design$C_l[lvl]) {
      # Buy as many as possible
      extra_samples <- floor(C_leftover / final_design$C_l[lvl])
      
      # Add them to the design
      final_design$n_l[lvl] <- final_design$n_l[lvl] + extra_samples
      
      # Deduct the spent money from the leftover budget
      C_leftover <- C_leftover - (extra_samples * final_design$C_l[lvl])
    }
  }
  
  return(final_design)
}

# Generates 1D nested coordinates for plotting
get_nested_design_1d <- function(n_vector, lower_bound = 0, upper_bound = 1) {
  L <- length(n_vector)
  X_nested <- list()
  
  X_nested[[1]] <- seq(lower_bound, upper_bound, length.out = n_vector[1])
  if (L > 1) {
    for (l in 2:L) {
      prev_n <- length(X_nested[[l - 1]])
      current_n <- n_vector[l]
      
      if (current_n == 1) {
        idx <- round((1 + prev_n) / 2)
      } else {
        idx <- round(seq(1, prev_n, length.out = current_n))
      }
      X_nested[[l]] <- X_nested[[l - 1]][idx]
    }
  }
  names(X_nested) <- paste0("Level_", 1:L)
  return(X_nested)
}


C_values <- c(500, 700, 900) 
all_designs_C <- list()

for (i in seq_along(C_values)) {
  curr_C <- C_values[i]
  
  # Calculate optimal allocation using the new inputs
  des_df <- calculate_new_design(C_total = curr_C, t_1 = 2, t_L = 0.1, 
                                 n_min = 3, lambda_max = 1, eta = 4, 
                                 c1 = 3.5, gamma = 1.6, d = 1)
  
  nested_pts <- get_nested_design_1d(des_df$n_l)
  curr_L <- nrow(des_df)
  curr_T <- des_df$T_val[1] # T is identical across the dataframe
  
  # Flatten to dataframe for ggplot
  temp_df <- data.frame(
    x = unlist(nested_pts, use.names = FALSE),
    t_l = rep(des_df$t_l, times = sapply(nested_pts, length)),
    Level = as.factor(rep(des_df$l, times = sapply(nested_pts, length)))
  )
  
  # Create a clean math label indicating the resulting L and T for this Budget
  label_str <- sprintf(
    "italic(C)[total] == %d * ',' ~ italic(L) == %d * ',' ~ italic(T) == '%.2f'",
    curr_C, curr_L, curr_T
  )
  temp_df$FacetLabel <- label_str
  
  all_designs_C[[i]] <- temp_df
}

plot_df_C <- do.call(rbind, all_designs_C)
plot_df_C$FacetLabel <- factor(plot_df_C$FacetLabel, levels = unique(plot_df_C$FacetLabel))

# Plotting
figure3 <- ggplot(plot_df_C, aes(x = x, y = t_l)) +
  geom_point(size = 3, show.legend = FALSE) +
  facet_wrap(~ FacetLabel, labeller = label_parsed, ncol = 3) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 2)) + 
  labs(x = "x", y = "t") +
  theme_bw(base_family = "serif", base_size = 14) +
  theme(strip.background = element_rect(fill = "lightgray"),
        panel.spacing = unit(1, "lines"))

print(figure3)
