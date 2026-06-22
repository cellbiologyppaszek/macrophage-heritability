# ==============================================================================
# Stochastic Colony Growth Simulation with Density-Dependent (Hill) k_on Rate
# Author: Apurv Srivastav
# Description: Gillespie simulation of cell colony growth with non-linear,
#              size-dependent phenotypic state transitions, matched against
#              experimental data with bootstrap-validated LOESS trend line profiles.
# ==============================================================================

# --- 1. DEPENDENCIES ---
library(readxl)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(patchwork)

# --- 2. GLOBAL CONSTANTS & HYPERPARAMETERS ---
set.seed(42)

# Growth & Time Params
GROWTH_RATE   <- log(2) / 14.4  
DOUBLING_TIME <- 14.4
T_END_SIM     <- 144          
NUM_COLONIES  <- 10000        

# Non-linear Transition Parameters (Hill Activation Function)
KMIN_ON       <- 0.001622386
K_PARAM       <- 2^6.4
N_HILL        <- 3   
KMAX_CONSTANT <- 0.01707356  
KOFF_CONSTANT <- 0.02277356

# Input/Output Configurations
DATA_PATH     <- "..\data\external\excel_inputs\IL1 and CD36 fractions repository.xlsx"
DATA_SHEET    <- "Cd36 fractions"
EPSILON       <- 1e-4
X_GRID        <- seq(2, 10.5, length.out = 100)

# --- 3. CORE CORE SIMULATION ENGINE ---

#' Gillespie Engine with Size-Dependent Phenotypic Transitions
#' @return Numeric vector containing final c(n_on, n_off)
grow_colony <- function(init, tend, kmin, kmax_val, K, n, koff_val) {
  n_on <- init[1]; n_off <- init[2]; t <- 0
  
  while (t < tend) {
    N <- n_on + n_off
    if (N <= 0) break
    
    # Calculate density-dependent activation rate (Hill Equation)
    curr_kon <- kmin + (kmax_val - kmin) * (N^n / (K^n + N^n))
    
    props <- c(
      GROWTH_RATE * n_on,  # Birth ON
      GROWTH_RATE * n_off, # Birth OFF
      koff_val * n_on,     # Switch ON -> OFF
      curr_kon * n_off     # Switch OFF -> ON
    )
    p_sum <- sum(props)
    if (is.na(p_sum) || p_sum <= 0) break
    
    t <- t + rexp(1, p_sum)
    if (t > tend) break
    
    rand <- runif(1) * p_sum
    if (rand <= props[1]) {
      n_on <- n_on + 1
    } else if (rand <= sum(props[1:2])) {
      n_off <- n_off + 1
    } else if (rand <= sum(props[1:3])) {
      n_on <- n_on - 1; n_off <- n_off + 1
    } else {
      n_on <- n_on + 1; n_off <- n_off - 1
    }
  }
  return(c(n_on = n_on, n_off = n_off))
}

#' Executes simulation loops over the defined colony universe
run_simulation_universe <- function(n_colonies) {
  message(sprintf("Running %s simulations with constant kmax...", format(n_colonies, big.mark=",")))
  
  results <- map_df(1:n_colonies, function(id) {
    # Phenotype initialization step
    init <- if(runif(1) <= 0.075) c(1, 0) else c(0, 1)
    
    out <- grow_colony(init, T_END_SIM, KMIN_ON, KMAX_CONSTANT, K_PARAM, N_HILL, KOFF_CONSTANT)
    tot <- sum(out)
    
    data.frame(
      Total_N       = tot, 
      Fraction_ON   = if(tot > 0) out[1]/tot else 0, 
      log2_N        = log2(tot),
      kmax_assigned = KMAX_CONSTANT
    )
  })
  return(results %>% filter(Total_N > 1))
}

# --- 4. BOOTSTRAP RESAMPLING FUNCTIONS ---

#' Computes LOESS models on resampled profiles
boot_loess <- function(data, x_col, y_col, iterations = 1000) {
  map_df(1:iterations, function(i) {
    sample_df <- data[sample(1:nrow(data), replace = TRUE), ]
    model_formula <- as.formula(paste(y_col, "~", x_col))
    m <- loess(model_formula, data = sample_df, span = 0.75)
    
    grid_df <- data.frame(new_x = X_GRID) %>% setNames(x_col)
    preds   <- predict(m, newdata = grid_df)
    
    data.frame(x = X_GRID, y = preds, iter = i)
  })
}

#' Calculates empirical confidence bounds from resampled iterations
get_bootstrap_ci <- function(boot_df) {
  boot_df %>%
    group_by(x) %>%
    summarise(
      mean_y  = mean(y, na.rm = TRUE),
      low     = quantile(y, 0.025, na.rm = TRUE),
      high    = quantile(y, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

# --- 5. PIPELINE PROCESSING WORKFLOW ---

# 5.1 Run Core Simulations
results_all <- run_simulation_universe(NUM_COLONIES)

# 5.2 Load and Parse Empirical Data Metrics
message("Loading and formatting raw experimental benchmarks...")
df_raw <- read_excel(DATA_PATH, sheet = DATA_SHEET, col_names = TRUE)
df_filtered <- df_raw[, colSums(is.na(df_raw)) != nrow(df_raw)]

CD36_Mock_fraction <- df_filtered %>%
  select(log2Number = 1, FractionON = 2) %>%
  mutate(
    log2Number = as.numeric(log2Number),
    FractionON = as.numeric(FractionON)
  ) %>%
  na.omit() %>%
  mutate(FractionON = pmin(pmax(FractionON, EPSILON), 1 - EPSILON))

# 5.3 Distribution-Matching Pipeline via Nearest Neighbors
message("Aligning simulated profiles via Nearest Neighbor matching...")
matched_indices <- map_int(CD36_Mock_fraction$log2Number, function(exp_val) {
  which.min(abs(results_all$log2_N - exp_val))
})

results_matched <- results_all[matched_indices, ] %>%
  mutate(log2_group = floor(log2_N))

# 5.4 Compute Stratified Statistical Metrics (Mean, CV Bounds)
message("Bootstrapping localized subset stats...")
boot_results <- map_df(1:1000, ~ results_matched %>% 
                         group_by(log2_group) %>% 
                         sample_frac(replace = TRUE) %>% 
                         summarise(b_mean = mean(Fraction_ON), 
                                   b_cv   = sd(Fraction_ON)/mean(Fraction_ON), 
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

# 5.5 Analytical Calculations for Analytical Steady-States
theo_curve_data <- data.frame(log2_N = seq(min(results_matched$log2_N), max(results_matched$log2_N), length.out = 300)) %>%
  mutate(
    N             = 2^log2_N,
    curr_kon      = KMIN_ON + (KMAX_CONSTANT - KMIN_ON) * (N^N_HILL / (K_PARAM^N_HILL + N^N_HILL)),
    theo_fraction = curr_kon / (curr_kon + KOFF_CONSTANT)
  )

# 5.6 Bootstrap Regression Lines with CIs
message("Executing Nonparametric Bootstrap LOESS loops...")
exp_boot <- boot_loess(CD36_Mock_fraction, "log2Number", "FractionON")
sim_boot <- boot_loess(results_matched, "log2_N", "Fraction_ON")

exp_ci <- get_bootstrap_ci(exp_boot)
sim_ci <- get_bootstrap_ci(sim_boot)

# --- 6. PLOT VISUALIZATION AND ANALYSIS PARSERS ---

# P1: Master Scatter Overlay + Local Trend lines + Theoretical Baseline
p_master <- ggplot() + 
  geom_jitter(data = results_matched, aes(x = log2_N, y = Fraction_ON, color = "Simulation (Matched)"), alpha = 0.3, width = 0.05) + 
  geom_point(data = CD36_Mock_fraction, aes(x = log2Number, y = FractionON, color = "Experimental"), size = 3, shape = 17) + 
  geom_smooth(data = results_matched, aes(x = log2_N, y = Fraction_ON), color = "blue", method = "loess", se = FALSE, linetype = "dashed", size = 1) + 
  geom_smooth(data = CD36_Mock_fraction, aes(x = log2Number, y = FractionON), color = "red", method = "loess", se = FALSE, size = 1.2) + 
  geom_line(data = theo_curve_data, aes(x = log2_N, y = theo_fraction, color = "Theoretical (Fixed kmax)"), size = 1.2) + 
  scale_color_manual(values = c("Simulation (Matched)" = "blue", "Experimental" = "red", "Theoretical (Fixed kmax)" = "black")) +
  labs(title = "Constant kmax Model: Simulation vs Experimental Comparison", subtitle = paste0("kmax fixed at ", round(KMAX_CONSTANT, 5)),
       x = "log2(N)", y = "Fraction ON", color = "Legend") +
  theme_minimal() + theme(legend.position = "bottom")

# P2: Empirical Density Validation Overlay
dist_data <- bind_rows(results_matched %>% mutate(Type = "Simulation (Matched)") %>% select(N = log2_N, Type),
                       CD36_Mock_fraction %>% mutate(Type = "Experimental") %>% select(N = log2Number, Type))

p_dist_validation <- ggplot(dist_data, aes(x = N, fill = Type)) +
  geom_density(alpha = 0.4) +
  labs(title = "Validation Plot: Size Profile Overlap Density", x = "log2(N)", y = "Density") +
  theme_minimal() + theme(legend.position = "bottom")

# P3: Localized Mean & CV Summary Profiles
p_mean <- ggplot(final_stats, aes(x = log2_group, y = mean_val)) + 
  geom_ribbon(aes(ymin = mean_low, ymax = mean_high), fill = "steelblue", alpha = 0.3) + 
  geom_line(color = "steelblue", linewidth = 1) + geom_point() +
  labs(title = "Mean Fraction ON Profile", x = "log2(N) Group", y = "Mean") + theme_minimal()

p_cv   <- ggplot(final_stats, aes(x = log2_group, y = cv_val)) + 
  geom_ribbon(aes(ymin = cv_low, ymax = cv_high), fill = "firebrick", alpha = 0.3) + 
  geom_line(color = "firebrick", linewidth = 1) + geom_point() +
  labs(title = "CV Coefficient Variation Profile", x = "log2(N) Group", y = "CV") + theme_minimal()

# P4: Clear Bootstrap Error Ribbon Comparative Layout
p_loess_ci <- ggplot() +
  geom_ribbon(data = exp_ci, aes(x = x, ymin = low, ymax = high), fill = "red", alpha = 0.2) +
  geom_line(data = exp_ci, aes(x = x, y = mean_y), color = "red", size = 1) +
  geom_ribbon(data = sim_ci, aes(x = x, ymin = low, ymax = high), fill = "blue", alpha = 0.2) +
  geom_line(data = sim_ci, aes(x = x, y = mean_y), color = "blue", size = 1, linetype = "dashed") +
  labs(title = "LOESS Curve Trends: Experiment vs. Simulation Bounds",
       subtitle = "Ribbons outline 95% Confidence Bounds via 1,000 Bootstraps", x = "log2(N)", y = "Fraction ON") +
  theme_minimal()

# Render Diagnostic Suite
print(p_master)
print(p_dist_validation)
print(p_mean / p_cv)
print(p_loess_ci)

# --- 7. FILE STORAGE DATA EXPORTS ---
message("Writing filtered outputs to disk...")
dataset_density <- results_matched %>% select(log2_N, Fraction_ON)
write.csv(dataset_density, "CD36_simulated_density_constant_rates.csv", row.names = FALSE)

sim_export <- sim_ci %>% select(log2N = x, mean = mean_y, `2.5%` = low, `97.5%` = high)
write.csv(sim_export, "simulated_loess_trend_data_kON_constant_wCI.csv", row.names = FALSE)