# ==============================================================================
# Constant-Rate Stochastic Lineage Engine (Gillespie Tracking)
# Author: Apurv Srivastav
# Description: Continuous-time simulation of single-cell colony expansions with 
#              constant phenotypic state transitions (independent of colony size). 
#              Tracks exact mother-daughter cell lineage histories.
# ==============================================================================

# --- 1. DEPENDENCIES & ENVIRONMENT SETUP ---
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

set.seed(42)

# --- 2. GLOBAL CONSTANTS & EXPERIMENTAL HYPERPARAMETERS ---
GROWTH_RATE   <- log(2) / 18  
T_END_SIM     <- 144   
NUM_COLONIES  <- 10
RECORD_TIMES  <- seq(0, T_END_SIM, by = 2) 

# Constant Phenotypic Transition Kinetic Rates
KOFF_CONSTANT <- 0.04429719 
KON_CONSTANT  <- 0.007508616 

# Analytical Steady-State Target Profile
THEORETICAL_FRAC <- KON_CONSTANT / (KON_CONSTANT + KOFF_CONSTANT)

# --- 3. HIGH-PERFORMANCE LINEAGE ENGINE ---

#' Simulates Single-Cell Trees Using Pre-Allocated Fixed Arrays
#' @return List holding an tracking logging dataframe of history mutations
grow_colony_lineage_const <- function(init_state, tend, kon_val, koff_val) {
  
  # Allocate large vector stores up front to avoid dataframe memory fragmenting
  max_estimated_events <- 25000
  v_time      <- numeric(max_estimated_events)
  v_event     <- character(max_estimated_events)
  v_parent    <- integer(max_estimated_events)
  v_d1        <- integer(max_estimated_events)
  v_d2        <- integer(max_estimated_events)
  v_cell      <- integer(max_estimated_events)
  v_new_state <- character(max_estimated_events)
  
  # Track active simulation cell arrays natively
  active_ids    <- c(1L)
  active_states <- c(init_state)
  
  t       <- 0
  next_id <- 2L
  log_idx <- 1L
  
  while (t < tend && length(active_ids) > 0) {
    N <- length(active_ids)
    
    # Vectorized rates calculations
    switch_rates <- ifelse(active_states == "ON", koff_val, kon_val)
    total_rates  <- GROWTH_RATE + switch_rates
    sum_all_rate <- sum(total_rates)
    
    if (is.na(sum_all_rate) || sum_all_rate <= 0) break
    
    t <- t + rexp(1, sum_all_rate)
    if (t > tend) break
    
    # Choose operational cell index via base probability arrays
    chosen_idx  <- sample.int(N, size = 1, prob = total_rates)
    chosen_id   <- active_ids[chosen_idx]
    chosen_stat <- active_states[chosen_idx]
    
    # Branching decision matrix logic
    if (runif(1) < (GROWTH_RATE / total_rates[chosen_idx])) {
      # Path A: Symmetric Mitotic Division
      v_time[log_idx]   <- t
      v_event[log_idx]  <- "division"
      v_parent[log_idx] <- chosen_id
      v_d1[log_idx]     <- next_id
      v_d2[log_idx]     <- next_id + 1L
      
      active_ids    <- c(active_ids[-chosen_idx], next_id, next_id + 1L)
      active_states <- c(active_states[-chosen_idx], chosen_stat, chosen_stat)
      next_id       <- next_id + 2L
    } else {
      # Path B: Stochastic State Toggle Conversion Step
      new_state <- ifelse(chosen_stat == "ON", "OFF", "ON")
      
      v_time[log_idx]      <- t
      v_event[log_idx]     <- "switch"
      v_cell[log_idx]      <- chosen_id
      v_new_state[log_idx] <- new_state
      
      active_states[chosen_idx] <- new_state
    }
    log_idx <- log_idx + 1L
    
    # Safe protection guard checking array boundary constraints
    if (log_idx > max_estimated_events) {
      expansion_blocks     <- 15000
      v_time      <- c(v_time, numeric(expansion_blocks))
      v_event     <- c(v_event, character(expansion_blocks))
      v_parent    <- c(v_parent, integer(expansion_blocks))
      v_d1        <- c(v_d1, integer(expansion_blocks))
      v_d2        <- c(v_d2, integer(expansion_blocks))
      v_cell      <- c(v_cell, integer(expansion_blocks))
      v_new_state <- c(v_new_state, character(expansion_blocks))
    }
  }
  
  valid_slices <- 1:(log_idx - 1L)
  lineage_df   <- data.frame(
    time      = v_time[valid_slices],
    event     = v_event[valid_slices],
    parent    = v_parent[valid_slices],
    d1        = v_d1[valid_slices],
    d2        = v_d2[valid_slices],
    cell      = v_cell[valid_slices],
    new_state = v_new_state[valid_slices],
    stringsAsFactors = FALSE
  )
  
  return(list(lineage = lineage_df))
}

#' Slices continuous event streams chronologically to build static snapshot tables
get_snapshot <- function(time_point, lineage_df, init_state) {
  states <- data.frame(
    cell_id    = 1L, 
    mother_id  = NA_integer_, 
    state      = init_state, 
    birth_time = 0
  )
  
  if (nrow(lineage_df) == 0) {
    states$snapshot_time <- time_point
    return(states)
  }
  
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
      
      states <- bind_rows(
        states,
        data.frame(cell_id    = c(event$d1, event$d2),
                   mother_id  = as.integer(event$parent),
                   state      = parent_state,
                   birth_time = event$time)
      )
    } else if (event$event == "switch") {
      states$state[states$cell_id == event$cell] <- event$new_state
    }
  }
  states$snapshot_time <- time_point
  return(states)
}

# --- 4. PIPELINE RUNNER PROCESSOR ---

message(sprintf("Executing constant rate simulation suite for %s colonies...", NUM_COLONIES))
all_snapshots_list <- vector("list", NUM_COLONIES)

for (col_id in 1:NUM_COLONIES) {
  init_state <- ifelse(runif(1) <= THEORETICAL_FRAC, "ON", "OFF")
  out <- grow_colony_lineage_const(init_state, T_END_SIM, KON_CONSTANT, KOFF_CONSTANT)
  
  lineage_df <- out$lineage
  
  # Gather validation metrics indexing parent-birth historical traces
  division_times <- lineage_df %>%
    filter(event == "division") %>%
    select(cell_id = parent, division_time = time)
  
  snaps <- map_df(RECORD_TIMES, ~ get_snapshot(.x, lineage_df, init_state)) %>%
    left_join(division_times, by = "cell_id") %>%
    mutate(colony_id = col_id)
  
  all_snapshots_list[[col_id]] <- snaps
}

# Consolidate complete snapshots tracking dataframes
final_table <- bind_rows(all_snapshots_list) %>%
  select(cell_id, mother_id, state, snapshot_time, colony_id, birth_time, division_time) %>%
  arrange(colony_id, snapshot_time, cell_id)

message("\nPreviewing final clean organized datatable output structure:")
print(head(final_table, 10))

# --- 5. VISUALIZATION COMPILATION PANEL ---

# Extract and shape summarizing metrics tables
snap_summary <- final_table %>%
  group_by(colony_id, snapshot_time) %>%
  summarise(frac_on = mean(state == "ON"), total_cells = n(), .groups = "drop")

snap_counts <- final_table %>%
  group_by(colony_id, snapshot_time, state) %>%
  summarise(n = n(), .groups = "drop")

# P1: Local Continuous Phenotype Ratio Expressions
p1 <- ggplot(snap_summary, aes(snapshot_time, frac_on)) +
  geom_line(color = "purple", linewidth = 1) +
  facet_wrap(~ colony_id) +
  labs(title = "Fraction ON Expressions Trajectories", subtitle = "Constant Rate System Model Configurations", x = "Time (Hours)", y = "Fraction ON") +
  theme_minimal()

# P2: Raw Spatial Abundance Tracking Comparisons (ON vs. OFF counts)
p2 <- ggplot(snap_counts, aes(snapshot_time, n, color = state)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("ON" = "dodgerblue3", "OFF" = "firebrick")) +
  facet_wrap(~ colony_id, scales = "free_y") +
  labs(title = "Absolute Cellular Abundance Transitions", subtitle = "ON vs. OFF Single-Agent Population Growth Profiles", x = "Time (Hours)", y = "Total Cells Count") +
  theme_minimal() + theme(legend.position = "bottom")

# P3: Net Growth Rates (Semi-log scale projection checking)
p3 <- ggplot(snap_summary, aes(snapshot_time, total_cells)) +
  geom_line(color = "black", linewidth = 1) +
  scale_y_log10() +
  facet_wrap(~ colony_id, scales = "free_y") +
  labs(title = "Net Colony Biomass Accumulation Profiles", subtitle = "Log10 Transformation Metric Check", x = "Time (Hours)", y = "Total Population Count (Log Scale)") +
  theme_minimal()

print(p1)
print(p2)
print(p3)