# ==============================================================================
# Stochastic Lineage Engine with Continuous-Time Agent Tracking
# Author: Apurv Srivastav
# Description: Simulates single-cell colony growth tracking exact dividing mother-
#              daughter relationships. State transitions are size-dependent, driven
#              by a Hill Equation model and Beta-distributed maximum activation rates.
# ==============================================================================

# --- 1. DEPENDENCIES & ENVIRONMENT SETUP ---
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

set.seed(42)

# --- 2. GLOBAL HYPERPARAMETERS & SYSTEM CONSTANTS ---
GROWTH_RATE   <- log(2) / 18  
T_END_SIM     <- 144   
NUM_COLONIES  <- 10
RECORD_TIMES  <- seq(0, T_END_SIM, by = 1)

# Activation Kinetics Parameters (Hill Function Input)
KMIN_ON       <- 7.041946e-04 
K_PARAM       <- 7.787997e+01
N_HILL        <- 9.5  

# Heterogeneity Bounds (Beta Model Distribution Parameters)
ALPHA_K       <- 0.261833506   
BETA_K        <- 1.146577844  
KMAX_SCALE    <- 0.177582098
KOFF_CONSTANT <- 0.02077356

# --- 3. HIGH-PERFORMANCE LINEAGE AGENT ENGINE ---

#' Simulates Continuously Logged Colony Trees via Vector Memory Allocations
#' @return List holding an event tracking dataframe log and initial phenotypic state
grow_colony_lineage <- function(init_state, tend, kmin, kmax_val, K, n, koff_val) {
  
  # Pre-allocate large flat vector stores to avoid row-binding memory penalties
  max_estimated_events <- 20000
  v_time      <- numeric(max_estimated_events)
  v_event     <- character(max_estimated_events)
  v_parent    <- integer(max_estimated_events)
  v_d1        <- integer(max_estimated_events)
  v_d2        <- integer(max_estimated_events)
  v_cell      <- integer(max_estimated_events)
  v_new_state <- character(max_estimated_events)
  
  # Tracks actively living agents in system
  active_ids    <- c(1L)
  active_states <- c(init_state)
  
  t <- 0
  next_id <- 2L
  log_idx <- 1L
  
  while (t < tend && length(active_ids) > 0) {
    N <- length(active_ids)
    curr_kon <- kmin + (kmax_val - kmin) * (N^n / (K^n + N^n))
    
    # Vectorized rate extraction assignment
    switch_rates <- ifelse(active_states == "ON", koff_val, curr_kon)
    total_rates  <- GROWTH_RATE + switch_rates
    sum_all_rate <- sum(total_rates)
    
    if (is.na(sum_all_rate) || sum_all_rate <= 0) break
    
    t <- t + rexp(1, sum_all_rate)
    if (t > tend) break
    
    # Randomly select candidate cell based on relative hazard scales
    chosen_idx  <- sample.int(N, size = 1, prob = total_rates)
    chosen_id   <- active_ids[chosen_idx]
    chosen_stat <- active_states[chosen_idx]
    
    # Decide event path: True Division vs Phenotypic State Switching
    if (runif(1) < (GROWTH_RATE / total_rates[chosen_idx])) {
      # 1. Process Cellular Division
      v_time[log_idx]   <- t
      v_event[log_idx]  <- "division"
      v_parent[log_idx] <- chosen_id
      v_d1[log_idx]     <- next_id
      v_d2[log_idx]     <- next_id + 1L
      
      # Update tracking elements
      active_ids    <- c(active_ids[-chosen_idx], next_id, next_id + 1L)
      active_states <- c(active_states[-chosen_idx], chosen_stat, chosen_stat)
      next_id       <- next_id + 2L
    } else {
      # 2. Process State Toggle Switch
      new_state <- ifelse(chosen_stat == "ON", "OFF", "ON")
      
      v_time[log_idx]      <- t
      v_event[log_idx]     <- "switch"
      v_cell[log_idx]      <- chosen_id
      v_new_state[log_idx] <- new_state
      
      active_states[chosen_idx] <- new_state
    }
    log_idx <- log_idx + 1L
    
    # Expand vector lengths dynamically if memory limits are breached
    if (log_idx > max_estimated_events) {
      expansion_size <- 10000
      v_time      <- c(v_time, numeric(expansion_size))
      v_event     <- c(v_event, character(expansion_size))
      v_parent    <- c(v_parent, integer(expansion_size))
      v_d1        <- c(v_d1, integer(expansion_size))
      v_d2        <- c(v_d2, integer(expansion_size))
      v_cell      <- c(v_cell, integer(expansion_size))
      v_new_state <- c(v_new_state, character(expansion_size))
    }
  }
  
  # Trim allocations down into a clean tracking dataframe
  valid_indices <- 1:(log_idx - 1L)
  lineage_df <- data.frame(
    time      = v_time[valid_indices],
    event     = v_event[valid_indices],
    parent    = v_parent[valid_indices],
    d1        = v_d1[valid_indices],
    d2        = v_d2[valid_indices],
    cell      = v_cell[valid_indices],
    new_state = v_new_state[valid_indices],
    stringsAsFactors = FALSE
  )
  
  return(list(lineage = lineage_df, init_state = init_state))
}

#' Slices continuous event log chronologically to build static snapshots
get_snapshot <- function(time_point, lineage_df, init_state) {
  states <- data.frame(cell_id = 1L, mother_id = NA_integer_, state = init_state, actual_birth_time = 0)
  
  if (nrow(lineage_df) == 0) {
    states$snapshot_time <- time_point
    return(states)
  }
  
  # Filter events up to target evaluation marker
  valid_events <- lineage_df %>% filter(time <= time_point)
  if (nrow(valid_events) == 0) {
    states$snapshot_time <- time_point
    return(states)
  }
  
  for (i in seq_len(nrow(valid_events))) {
    event <- valid_events[i, ]
    if (event$event == "division") {
      p_idx <- which(states$cell_id == event$parent)
      parent_state <- states$state[p_idx]
      states <- states[-p_idx, ]
      
      # Re-insert generated daughters to track lineage tree branches
      states <- bind_rows(
        states,
        data.frame(cell_id = c(event$d1, event$d2),
                   mother_id = as.integer(event$parent),
                   state = parent_state,
                   actual_birth_time = event$time)
      )
    } else if (event$event == "switch") {
      states$state[states$cell_id == event$cell] <- event$new_state
    }
  }
  states$snapshot_time <- time_point
  return(states)
}

# --- 4. EXECUTION PIPELINE ---

message(sprintf("Running lineage tracking simulation across %s colonies...", NUM_COLONIES))
all_snapshots_list <- vector("list", NUM_COLONIES)

for (col_id in 1:NUM_COLONIES) {
  init_state <- ifelse(runif(1) <= 0.075, "ON", "OFF")
  kmax_val   <- rbeta(1, ALPHA_K, BETA_K) * KMAX_SCALE
  
  out <- grow_colony_lineage(init_state, T_END_SIM, KMIN_ON, kmax_val, K_PARAM, N_HILL, KOFF_CONSTANT)
  
  # Isolate dividing parents to append actual final division markers
  div_times <- out$lineage %>%
    filter(event == "division") %>%
    select(cell_id = parent, division_time = time)
  
  # Slice arrays over evaluated time vectors
  snaps <- map_df(RECORD_TIMES, ~ get_snapshot(.x, out$lineage, init_state)) %>%
    left_join(div_times, by = "cell_id") %>%
    mutate(colony_id = col_id)
  
  all_snapshots_list[[col_id]] <- snaps
}

# Final Data Consolidation
snapshots_all <- bind_rows(all_snapshots_list)

final_table <- snapshots_all %>%
  select(cell_id, mother_id, state, snapshot_time, colony_id, actual_birth_time, division_time) %>%
  arrange(colony_id, snapshot_time, cell_id)

message("\nPreviewing clean agent output file:")
print(head(final_table, 10))

# --- 5. AUTOMATED EXPORT & STORAGE ---
# write.csv(final_table, "simulated_10colonies_time144h_CD36.csv", row.names = FALSE)

# --- 6. PLOT VISUALIZATION SUITE ---

snap_summary <- snapshots_all %>%
  group_by(colony_id, snapshot_time) %>%
  summarise(frac_on = mean(state == "ON"), total_cells = n(), .groups = "drop")

# P1: Local Fraction ON Expression Paths
p1 <- ggplot(snap_summary, aes(snapshot_time, frac_on)) +
  geom_line(color = "darkorchid4", linewidth = 1) +
  facet_wrap(~ colony_id) +
  labs(title = "Phenotypic Fraction ON Paths", subtitle = "Hill Function Agent Model Tracking", x = "Time (Hours)", y = "Fraction ON") +
  theme_minimal()

# P2: Exponential Growth Profiling (Semi-Log scale check)
p2 <- ggplot(snap_summary, aes(snapshot_time, total_cells)) +
  geom_line(color = "grey10", linewidth = 1) +
  scale_y_log10() +
  facet_wrap(~ colony_id, scales = "free_y") +
  labs(title = "Colony Kinetics Population Scaling Profile", subtitle = "Log10 Transformation Metric Check", x = "Time (Hours)", y = "Total Living Cells (Count)") +
  theme_minimal()

print(p1)
print(p2)