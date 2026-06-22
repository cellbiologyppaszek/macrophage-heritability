# ==============================================================================
# Stochastic Colony Growth Simulation & Experimental Matching
# Author: Apurv Srivastav
# Description: Gillespie simulation of cell colony growth with stochastic 
#              ON/OFF switching phenotypic transitions, matched against 
#              experimental data distributions with bootstrap LOESS trends.
# ==============================================================================

# --- 1. DEPENDENCIES ---
library(readxl)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(patchwork)

# --- 2. GLOBAL PARAMETERS ---
set.seed(42)

# Growth & Transition Rates
STATE_ON_GROWTH  <- log(2) / 14.4
STATE_OFF_GROWTH <- log(2) / 14.4
SWITCH_ON_TO_OFF <- 0.05493723  # k_off
SWITCH_OFF_TO_ON <- 0.009457892 # k_on

SIM_PARAMS   <- c(STATE_ON_GROWTH, STATE_OFF_GROWTH, SWITCH_ON_TO_OFF, SWITCH_OFF_TO_ON)
NUM_COLONIES <- 10000 
T_END_SIM    <- 144   

# File Paths & Settings
DATA_PATH  <- "..\data\external\excel_inputs\IL1 and CD36 fractions repository.xlsx"
DATA_SHEET <- "Cd36 fractions" # Change to "IL1B fraction" as needed
EPSILON    <- 1e-4

# Theoretical steady state fraction
THEORETICAL_FRAC <- SWITCH_OFF_TO_ON / (SWITCH_ON_TO_OFF + SWITCH_OFF_TO_ON)

# --- 3. CORE SIMULATION FUNCTIONS ---

#' Gillespie Algorithm for Single Colony Growth
#' @param init Vector of c(n_on, n_off)
#' @param params Vector of c(g_on, g_off, k_off, k_on)
#' @param tend End time of simulation
grow_colony <- function(init, params, tend) {
  n_on  <- init[1]
  n_off <- init[2]
  t     <- 0
  
  g_on  <- params[1]
  g_off <- params[2]
  k_off <- params[3]
  k_on  <- params[4]
  
  while (t <= tend) {
    props <- c(
      g_on * n_on,    # Birth ON
      g_off * n_off,  # Birth OFF
      k_off * n_on,   # Switch ON -> OFF
      k_on * n_off    # Switch OFF -> ON
    )
    
    prop_sum <- sum(props)
    if (prop_sum <= 0) break
    
    t <- t + rexp(1, rate = prop_sum)
    if (t > tend) break
    
    rand <- runif(1) * prop_sum
    if (rand <= props[1]) {
      n_on <- n_on + 1
    } else if (rand <= sum(props[1:2])) {
      n_off <- n_off + 1
    } else if (rand <= sum(props[1:3])) {
      n_on <- n_on - 1
      n_off <- n_off + 1
    } else {
      n_on <- n_on + 1
      n_off <- n_off - 1
    }
  }
  return(c(n_on = n_on, n_off = n_off))
}

#' Generate Large Pool of Simulated Colonies
run_stochastic_simulation <- function(num_colonies, params, tend, target_frac) {
  message(sprintf("Simulating %s colonies via Gillespie Algorithm...", format(num_colonies, big.mark=",")))
  
  results <- map_df(1:num_colonies, function(i) {
    init <- if(runif(1) <= target_frac) c(1, 0) else c(0, 1)
    out  <- grow_colony(init, params, tend)
    tot  <- sum(out)
    
    data.frame(
      Total_N     = tot,
      log2_N      = log2(tot),
      Fraction_ON = if(tot > 0) out[1]/tot else 0
    )
  })
  
  return(results %>% filter(Total_N > 1))
}

# --- 4. BOOTSTRAP LOESS FUNCTIONS ---

#' Bootstrap LOESS Regressions over a standard grid
boot_loess <- function(data, x_col, y_col, x_grid, iterations = 1000) {
  map_df(1:iterations, function(i) {
    sample_df <- data[sample(1:nrow(data), replace = TRUE), ]
    
    # Dynamic formula construction
    model_formula <- as.formula(paste(y_col, "~", x_col))
    m <- loess(model_formula, data = sample_df, span = 0.75)
    
    # Predict over standard grid
    grid_df <- data.frame(new_x = x_grid) %>% setNames(x_col)
    preds   <- predict(m, newdata = grid_df)
    
    data.frame(x = x_grid, y = preds, iter = i)
  })
}

#' Calculate Confidence Intervals from Bootstrapped Data
get_bootstrap_ci <- function(boot_df) {
  boot_df %>%
    group_by(x) %>%
    summarise(
      mean_y = mean(y, na.rm = TRUE),
      low    = quantile(y, 0.025, na.rm = TRUE),
      high   = quantile(y, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

# --- 5. PIPELINE EXECUTION ---

# 5.1 Run Simulation Pool
colony_pool <- run_stochastic_simulation(NUM_COLONIES, SIM_PARAMS, T_END_SIM, THEORETICAL_FRAC)

# 5.2 Load and Clean Experimental Data
message("Loading and cleaning experimental data...")
raw_excel_data <- read_excel(DATA_PATH, sheet = DATA_SHEET, col_names = TRUE)

# Remove completely empty columns, handle data types and bound 0/1 issues
experimental_data <- raw_excel_data[, colSums(is.na(raw_excel_data)) != nrow(raw_excel_data)] %>%
  select(log2Number = 1, FractionON = 2) %>%
  mutate(
    log2Number = as.numeric(log2Number),
    FractionON = as.numeric(FractionON)
  ) %>%
  na.omit() %>%
  mutate(FractionON = pmin(pmax(FractionON, EPSILON), 1 - EPSILON))

# 5.3 Nearest-Neighbor Matching
message("Matching simulation pool to experimental growth profile...")
matched_indices <- map_int(experimental_data$log2Number, function(exp_val) {
  which.min(abs(colony_pool$log2_N - exp_val))
})

results_matched <- colony_pool[matched_indices, ] %>%
  mutate(log2_group = floor(log2_N))

# 5.4 Grouped Bootstrapping (Mean & CV Analysis)
message("Bootstrapping matched subset for Mean/CV statistics...")
boot_results <- map_df(1:1000, ~ results_matched %>% 
                         group_by(log2_group) %>% 
                         sample_frac(replace = TRUE) %>% 
                         summarise(b_mean = mean(Fraction_ON), 
                                   b_cv  = sd(Fraction_ON)/mean(Fraction_ON), 
                                   .groups="drop"))

final_stats <- boot_results %>% 
  group_by(log2_group) %>% 
  summarise(mean_val  = mean(b_mean, na.rm=T), 
            mean_low  = quantile(b_mean, 0.025, na.rm=T), 
            mean_high = quantile(b_mean, 0.975, na.rm=T),
            cv_val    = mean(b_cv, na.rm=T), 
            cv_low    = quantile(b_cv, 0.025, na.rm=T), 
            cv_high   = quantile(b_cv, 0.975, na.rm=T),
            .groups   = "drop")

# 5.5 Bootstrap LOESS Trends
message("Executing Bootstrap LOESS Trends...")
x_grid   <- seq(2, 10.5, length.out = 100)
exp_boot <- boot_loess(experimental_data, "log2Number", "FractionON", x_grid)
sim_boot <- boot_loess(results_matched, "log2_N", "Fraction_ON", x_grid)

exp_ci   <- get_bootstrap_ci(exp_boot)
sim_ci   <- get_bootstrap_ci(sim_boot)

# --- 6. PLOTTING & DATA EXPORT ---

# P1: Main Scatter & Local Fits
p1 <- ggplot() + 
  geom_jitter(data = results_matched, aes(x = log2_N, y = Fraction_ON, color = "Simulation (Matched)"), alpha = 0.4, width = 0.05) +
  geom_point(data = experimental_data, aes(x = log2Number, y = FractionON, color = "Clone"), size = 3, shape = 17) +
  geom_smooth(data = results_matched, aes(x = log2_N, y = Fraction_ON), color = "blue", method = "loess", se = FALSE, linetype = "dashed", size = 1) + 
  geom_smooth(data = experimental_data, aes(x = log2Number, y = FractionON), color = "red", method = "loess", se = FALSE, linetype = "dashed", size = 1) +
  scale_color_manual(values = c("Simulation (Matched)" = "blue", "Clone" = "red")) +
  labs(title = "Matched Population: Constant Rates vs Experimental", x = "log2(N)", y = "Fraction ON") +
  theme_minimal() + theme(legend.position = "bottom")

# P2: Distribution Alignment Check
p_dist_N <- ggplot(bind_rows(results_matched %>% mutate(Type = "Simulation (Matched)") %>% select(N = log2_N, Type),
                             experimental_data %>% mutate(Type = "Experimental") %>% select(N = log2Number, Type)), 
                   aes(x = N, fill = Type)) +
  geom_density(alpha = 0.4) +
  labs(title = "Validation: log2(N) Overlap", x = "log2(N)", y = "Density") +
  theme_minimal() + theme(legend.position = "bottom")

# P3: Combined Bootstrap LOESS Trend Comparison Plot
p_loess_ci <- ggplot() +
  geom_ribbon(data = exp_ci, aes(x = x, ymin = low, ymax = high), fill = "red", alpha = 0.2) +
  geom_line(data = exp_ci, aes(x = x, y = mean_y), color = "red", size = 1) +
  geom_ribbon(data = sim_ci, aes(x = x, ymin = low, ymax = high), fill = "blue", alpha = 0.2) +
  geom_line(data = sim_ci, aes(x = x, y = mean_y), color = "blue", size = 1, linetype = "dashed") +
  labs(title = "LOESS Trend Line: Experiment vs. Simulation",
       subtitle = "Ribbons represent 95% Confidence Intervals (1,000 bootstrap iterations)",
       x = "log2(N)", y = "Fraction ON") +
  theme_minimal()

# Display plots sequentially or assemble via patchwork if desired
print(p1)
print(p_dist_N)
print(p_loess_ci)

# Optional Export Data Prep
# sim_export <- sim_ci %>% select(log2N = x, mean = mean_y, `2.5%` = low, `97.5%` = high)
# exp_export <- exp_ci %>% select(log2N = x, mean = mean_y, `2.5%` = low, `97.5%` = high)
# write.csv(sim_export, "simulated_loess_trend_data_constant_rates_IL1B_wCI.csv", row.names = FALSE)