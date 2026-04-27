setwd("~/Library/CloudStorage/OneDrive-NTNU/NARM/MSNARM/Footprint/Data/R.master/leirelv")

library(tidyverse)
library(MASS)
library(mgcv)
library(sf)
library(performance)
library(readxl)
library(ggplot2)
library(unmarked)


library(future.apply)

n_workers <- max(1, parallel::detectCores() - 1)
future::plan(future::multisession, workers = n_workers)

#load data:

#baseline data moose 
moose_baseline <- st_read("moose4.shp")
#baseline data roedeer
roedeer_baseline <- st_read("roedeer4.shp")
#baseline data redfox 
redfox_baseline <- st_read("redfox4.shp")

#load species SDMs 

model_moose_sdm2 <- readRDS("baselines/model_moose_sdm_original.rds")
model_roedeer_sdm2 <- readRDS("baselines/model_roedeer_sdm_original.rds")
model_redfox_sdm2 <- readRDS("baselines/model_redfox_sdm_original.rds")

#-----BASELINE MODELS AND PREDICTIONS -------

summary(model_moose_sdm2)
summary(model_roedeer_sdm2)
summary(model_redfox_sdm2)

moose_baseline$p_true <- predict(model_moose_sdm2,
                                 newdata = moose_baseline,
                                 type = "response")

redfox_baseline$p_true <- predict(model_redfox_sdm2,
                                  newdata = redfox_baseline,
                                  type = "response")

roedeer_baseline$p_true <- predict(model_roedeer_sdm2,
                                   newdata = roedeer_baseline,
                                   type = "response")


#subset to studyarea:
moose_studyarea <- moose_baseline[moose_baseline$studyarea != 0, ]
roedeer_studyarea <- roedeer_baseline[roedeer_baseline$studyarea != 0, ]
redfox_studyarea <- redfox_baseline[redfox_baseline$studyarea != 0, ]

#-------Adjusting the mean ------------------------

logit <- function(p) qlogis(p)
inv_logit <- function(x) plogis(x)

# Solve for a constant delta (intercept shift) so the mean drops by `drop`
find_delta_for_mean_drop <- function(p, drop = 0.2, eps = 1e-6) {
  p <- as.numeric(p)
  p <- p[is.finite(p)]
  if (length(p) == 0) stop("No finite p values.")

  p_clamp <- pmin(pmax(p, eps), 1 - eps)
  
  mu0 <- mean(p_clamp)
  mu_target <- mu0 - drop
  mu_target <- max(min(mu_target, 1 - eps), eps)
  
  f <- function(delta) mean(inv_logit(logit(p_clamp) + delta)) - mu_target
  
  uniroot(f, interval = c(-50, 0))$root
}

apply_delta <- function(p, delta, eps = 1e-6) {
  p <- as.numeric(p)
  p <- pmin(pmax(p, eps), 1 - eps)
  inv_logit(logit(p) + delta)
}

# ---- for each species within the study area ----
drop <- 0.2

drop_moose <- 0.4

# Moose
delta_moose <- find_delta_for_mean_drop(moose_studyarea$p_true, drop = drop_moose)
moose_studyarea$p_true_minus02 <- apply_delta(moose_studyarea$p_true, delta_moose)

# Roe deer
delta_roedeer <- find_delta_for_mean_drop(roedeer_studyarea$p_true, drop = drop)
roedeer_studyarea$p_true_minus02 <- apply_delta(roedeer_studyarea$p_true, delta_roedeer)

# Red fox
delta_redfox <- find_delta_for_mean_drop(redfox_studyarea$p_true, drop = drop)
redfox_studyarea$p_true_minus02 <- apply_delta(redfox_studyarea$p_true, delta_redfox)

# ---- Sanity checks ----
check <- data.frame(
  species = c("Moose", "Roe deer", "Red fox"),
  mean_before = c(mean(moose_studyarea$p_true, na.rm=TRUE),
                  mean(roedeer_studyarea$p_true, na.rm=TRUE),
                  mean(redfox_studyarea$p_true, na.rm=TRUE)),
  mean_after  = c(mean(moose_studyarea$p_true_minus02, na.rm=TRUE),
                  mean(roedeer_studyarea$p_true_minus02, na.rm=TRUE),
                  mean(redfox_studyarea$p_true_minus02, na.rm=TRUE))
)
check$drop_achieved <- check$mean_before - check$mean_after
check

delta_moose
delta_roedeer
delta_redfox

summary(moose_studyarea$p_true_minus02)
summary(roedeer_studyarea$p_true_minus02)
summary(redfox_studyarea$p_true_minus02)


#--------- Add fixed corridor effect to studyarea subset p_true_minus02-----------------
logit <- function(p) qlogis(p)
inv_logit <- function(x) plogis(x)

corridor_OR <- 1.3
beta_corridor <- log(corridor_OR)

eps <- 1e-6

add_corridor_effect <- function(df, p_col = "p_true_minus02",
                                studyarea_col = "studyarea",
                                corridor_value = 2,
                                out_col = "p_true_adj") {
  p <- as.numeric(df[[p_col]])
  p <- pmin(pmax(p, eps), 1 - eps)  
  
  corridor_flag <- as.integer(df[[studyarea_col]] == corridor_value)
  
  df[[out_col]] <- inv_logit(logit(p) + beta_corridor * corridor_flag)
  df[["corridor_flag"]] <- corridor_flag  
  df
}

# Apply to study-area subsets (studyarea != 0 already)
moose_studyarea   <- add_corridor_effect(moose_studyarea)
roedeer_studyarea <- add_corridor_effect(roedeer_studyarea)
redfox_studyarea  <- add_corridor_effect(redfox_studyarea)

# Quick checks
summary(moose_studyarea$p_true_adj)
summary(roedeer_studyarea$p_true_adj)
summary(redfox_studyarea$p_true_adj)
table(moose_studyarea$corridor_flag)



#-------------SIMULATION SETUP-----------------------
# ============================================================
# Random clustered site selection scenario (clusters of 9)
# Spatial bias ONLY in site selection. Visit allocation unchanged.
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(future.apply)
  library(unmarked)
})

# ---- SIMULATION SETUP -------------------------------------------------------
p_detect_vals  <- c(0.6, 0.7, 0.8, 0.9, 0.999)

n_samples_vals <- c(200, 300, 400, 600, 1000)

n_reps <- 1000

visit_budget_map <- c(
  "200" =  400,
  "300" =  600,
  "400" =  800,
  "600" = 1200,
  "1000" = 2000
)

# true corridor effect
corridor_OR <- 1.3
beta_corridor_true <- log(corridor_OR)

# occupancy (state) covariates used in the fitted model
state_covariates <- c(
  "corridor_flag",
  "building_c",
  "road_area",
  "areatype20",
  "areatype30",
  "areatype50",
  "areatype60",
  "areatype81"
)

# ---- SPECIES DATA -----------------------------------------------------------
species_data_study_sf <- redfox_studyarea
species_name <- "Red fox"

seed_map <- c("Moose" = 100, "Roe deer" = 200, "Red fox" = 300)
set.seed(seed_map[[species_name]])

# ---- ADD CENTROID COORDINATES (for clustering) ------------------------------
if (sf::st_is_longlat(species_data_study_sf)) {
  stop("species_data_study_sf is lon/lat. Transform it to a projected CRS in meters before running.")
}

xy <- sf::st_coordinates(sf::st_centroid(species_data_study_sf))
species_data_study_sf$X <- xy[, 1]
species_data_study_sf$Y <- xy[, 2]

# ---- DROP GEOMETRY FOR unmarked --------------------------------------------
species_data_study <- species_data_study_sf %>%
  sf::st_drop_geometry() %>%
  as.data.frame()

# ---- SCALE CONTINUOUS COVARIATES ONCE --------------------------------------
cont_covs <- c("road_area", "areatype20", "areatype30", "areatype50", "areatype60", "areatype81")
missing_covs <- setdiff(c(cont_covs, "p_true_adj", "X", "Y"), names(species_data_study))
if (length(missing_covs) > 0) {
  stop("Missing columns in species_data_study: ", paste(missing_covs, collapse = ", "))
}

species_data_study <- species_data_study %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(cont_covs), ~ as.numeric(scale(.x))))

# ---- HELPERS ----------------------------------------------------------------
clamp01 <- function(p, eps = 1e-9) pmin(pmax(p, eps), 1 - eps)

allocate_visits <- function(n_sites, total_visits, cap_visits = Inf) {
  if (total_visits < n_sites) stop("Need total_visits >= n_sites (>=1 visit per site).")
  if (is.finite(cap_visits) && total_visits > n_sites * cap_visits) {
    stop("total_visits exceeds n_sites*cap_visits; adjust cap or budget.")
  }
  
  # baseline: 1 visit per site
  K <- rep.int(1L, n_sites)
  remaining <- total_visits - n_sites
  
  # distribute remaining uniformly at random across sites
  if (remaining > 0) {
    extra <- as.vector(rmultinom(1, size = remaining, prob = rep(1, n_sites)))
    K <- K + extra
  }
  
  # enforce cap_visits
  if (is.finite(cap_visits)) {
    cap_visits <- as.integer(cap_visits)
    overflow <- sum(pmax(K - cap_visits, 0L))
    K <- pmin(K, cap_visits)
    
    while (overflow > 0) {
      idx <- which(K < cap_visits)
      if (length(idx) == 0) stop("Could not redistribute overflow under cap_visits.")
      add <- as.vector(rmultinom(1, size = overflow, prob = rep(1, length(idx))))
      K[idx] <- pmin(cap_visits, K[idx] + add)
      overflow <- sum(pmax(K - cap_visits, 0L))
      K <- pmin(K, cap_visits)
    }
  }
  
  stopifnot(sum(K) == total_visits)
  K
}

# ---- CLUSTERED RANDOM SITE SELECTOR -----------------------------------------
# Returns n_samples unique row indices, clustered in groups of 'cluster_size'
# Cluster centers are drawn uniformly from ALL cells.
select_sites_clustered_random <- function(df, n_samples, cluster_size = 9) {
  if (cluster_size < 2) stop("cluster_size must be >= 2.")
  if (n_samples > nrow(df)) stop("n_samples exceeds total cells (no replacement).")
  if (!all(c("X", "Y") %in% names(df))) stop("Missing X/Y coordinates in df (needed for clustering).")
  
  picked <- integer(0)
  available <- seq_len(nrow(df))
  
  while (length(picked) < n_samples) {
    remaining <- n_samples - length(picked)
    k <- min(cluster_size, remaining)
    
    # pick a random cluster center from remaining cells
    center <- sample(available, 1)
    
    # take k nearest available cells (including center)
    cx <- df$X[center]; cy <- df$Y[center]
    ax <- df$X[available]; ay <- df$Y[available]
    d2 <- (ax - cx)^2 + (ay - cy)^2
    
    ord <- order(d2)
    cluster <- available[ord][seq_len(k)]
    
    picked <- c(picked, cluster)
    available <- setdiff(available, cluster)
    
    if (length(available) == 0 && length(picked) < n_samples) {
      stop("Ran out of cells before reaching n_samples (unexpected).")
    }
  }
  
  stopifnot(length(picked) == n_samples, length(unique(picked)) == n_samples)
  picked
}

# ---- ONE REPLICATE: random clustered site selection --------------------------
run_one_sim_clustered_random <- function(species_data_study,
                                         n_samples,
                                         total_visits,
                                         p_detect,
                                         state_covariates,
                                         cluster_size = 9,
                                         cap_visits = 20) {
  
  # 1) clustered random site selection
  sampled_idx <- select_sites_clustered_random(
    df = species_data_study,
    n_samples = n_samples,
    cluster_size = cluster_size
  )
  sampled_data <- species_data_study[sampled_idx, , drop = FALSE]
  
  # 2) simulate true occupancy
  sampled_data$z_true <- rbinom(n_samples, 1, prob = clamp01(sampled_data$p_true_adj))
  
  # 3) allocate visits across sampled sites 
  K_i <- allocate_visits(n_sites = n_samples, total_visits = total_visits, cap_visits = cap_visits)
  K_max <- max(K_i)
  
  # 4) detection histories (NA = no visit)
  y_mat <- matrix(NA_integer_, nrow = n_samples, ncol = K_max)
  ii <- rep.int(seq_len(n_samples), times = K_i)
  jj <- sequence(K_i)
  p_ij <- sampled_data$z_true[ii] * p_detect
  
  y_mat[cbind(ii, jj)] <- rbinom(length(p_ij), 1, prob = p_ij)
  colnames(y_mat) <- paste0("y_", seq_len(K_max))
  
  detected_any <- rowSums(y_mat == 1, na.rm = TRUE) > 0
  n_sampled_inside  <- sum(sampled_data$corridor_flag == 1)
  n_sampled_outside <- sum(sampled_data$corridor_flag == 0)
  n_detected_inside  <- sum(detected_any & sampled_data$corridor_flag == 1)
  n_detected_outside <- sum(detected_any & sampled_data$corridor_flag == 0)
  
  # 5) fit occupancy model
  site_covs <- sampled_data[, state_covariates, drop = FALSE]
  umf <- unmarked::unmarkedFrameOccu(y = y_mat, siteCovs = site_covs)
  
  state_rhs <- paste(state_covariates, collapse = " + ")
  full_form <- stats::as.formula(paste("~ 1 ~", state_rhs))
  
  fit <- tryCatch(unmarked::occu(full_form, data = umf), error = function(e) e)
  
  if (inherits(fit, "error")) {
    return(list(
      beta_hat = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_hat = NA_real_,
      fit_failed = TRUE,
      fit_error = conditionMessage(fit),
      
      n_sampled_inside = n_sampled_inside,
      n_sampled_outside = n_sampled_outside,
      n_detected_inside = n_detected_inside,
      n_detected_outside = n_detected_outside,
      
      K_mean = mean(K_i),
      K_prop1 = mean(K_i == 1)
    ))
  }
  
  # 6) extract estimates robustly 
  extract <- tryCatch({
    state_coef <- unmarked::coef(fit, type = "state")
    state_vcov <- unmarked::vcov(fit, type = "state")
    
    j <- grep("corridor_flag", names(state_coef))
    if (length(j) != 1) stop("Could not uniquely find corridor_flag.")
    
    est <- unname(state_coef[j])
    
    vj <- unname(diag(state_vcov)[j])
    if (!is.finite(vj) || vj <= 0) stop("Non-positive variance for corridor_flag (vcov diagonal).")
    
    se <- sqrt(vj)
    ci <- est + c(-1, 1) * 1.96 * se
    
    alpha_hat <- unmarked::coef(fit, type = "det")[1]
    p_hat <- plogis(alpha_hat)
    
    list(
      beta_hat = est,
      se = se,
      ci_low = unname(ci[1]),
      ci_high = unname(ci[2]),
      p_hat = unname(p_hat),
      fit_failed = FALSE,
      fit_error = NA_character_
    )
  }, error = function(e) {
    list(
      beta_hat = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_hat = NA_real_,
      fit_failed = TRUE, fit_error = conditionMessage(e)
    )
  })
  
  c(extract, list(
    n_sampled_inside = n_sampled_inside,
    n_sampled_outside = n_sampled_outside,
    n_detected_inside = n_detected_inside,
    n_detected_outside = n_detected_outside,
    
    K_mean = mean(K_i),
    K_prop1 = mean(K_i == 1)
  ))
}

# ---- RUN SCENARIO GRID ------------------------------------------------------
scenario_ID <- "CL3R"
scenario_name <- "Random clustered site selection (clusters of 9), occupancy (random visits; fixed budget)"
cluster_size <- 9

results <- list()
idx <- 1

for (p_detect in p_detect_vals) {
  for (n_samp in n_samples_vals) {
    
    total_visits <- unname(visit_budget_map[as.character(n_samp)])
    if (is.na(total_visits)) stop("No total visit budget found for n_samples=", n_samp)
    
    outs <- future.apply::future_lapply(
      X = seq_len(n_reps),
      FUN = function(r) {
        run_one_sim_clustered_random(
          species_data_study = species_data_study,
          n_samples = n_samp,
          total_visits = total_visits,
          p_detect = p_detect,
          state_covariates = state_covariates,
          cluster_size = cluster_size,
          cap_visits = 20
        )
      },
      future.seed = TRUE
    )
    
    # collect outputs
    betas <- vapply(outs, function(o) o$beta_hat, numeric(1))
    ses   <- vapply(outs, function(o) o$se, numeric(1))
    cil   <- vapply(outs, function(o) o$ci_low, numeric(1))
    cih   <- vapply(outs, function(o) o$ci_high, numeric(1))
    
    p_hats <- vapply(outs, function(o) o$p_hat, numeric(1))
    
    K_mean_vec  <- vapply(outs, function(o) o$K_mean, numeric(1))
    K_prop1_vec <- vapply(outs, function(o) o$K_prop1, numeric(1))
    
    n_sampled_inside_vec   <- vapply(outs, function(o) o$n_sampled_inside, numeric(1))
    n_sampled_outside_vec  <- vapply(outs, function(o) o$n_sampled_outside, numeric(1))
    n_detected_inside_vec  <- vapply(outs, function(o) o$n_detected_inside, numeric(1))
    n_detected_outside_vec <- vapply(outs, function(o) o$n_detected_outside, numeric(1))
    
    fit_failed_vec <- vapply(outs, function(o) isTRUE(o$fit_failed), logical(1))
    fit_error_vec  <- vapply(outs, function(o) if (!is.null(o$fit_error)) o$fit_error else NA_character_, character(1))
    
    # performance metrics (corridor effect)
    bias <- mean(betas - beta_corridor_true, na.rm = TRUE)
    rmse <- sqrt(mean((betas - beta_corridor_true)^2, na.rm = TRUE))
    coverage <- mean(cil <= beta_corridor_true & cih >= beta_corridor_true, na.rm = TRUE)
    power <- mean(cil > 0 | cih < 0, na.rm = TRUE)
    typeS <- mean(sign(betas) != sign(beta_corridor_true), na.rm = TRUE)
    
    # detection performance
    p_bias <- mean(p_hats - p_detect, na.rm = TRUE)
    p_rmse <- sqrt(mean((p_hats - p_detect)^2, na.rm = TRUE))
    
    n_failed  <- sum(fit_failed_vec, na.rm = TRUE)
    fail_rate <- mean(fit_failed_vec, na.rm = TRUE)
    top_error <- if (n_failed > 0) {
      names(sort(table(fit_error_vec[fit_failed_vec]), decreasing = TRUE))[1]
    } else {
      NA_character_
    }
    
    results[[idx]] <- data.frame(
      scenario_ID  = scenario_ID,
      species      = species_name,
      p_detect     = p_detect,
      n_samples    = n_samp,
      total_visits = total_visits,
      
      bias        = bias,
      rmse        = rmse,
      coverage_95 = coverage,
      power_95    = power,
      typeS_error = typeS,
      mean_beta_hat = mean(betas, na.rm = TRUE),
      sd_beta_hat   = sd(betas, na.rm = TRUE),
      
      mean_p_hat = mean(p_hats, na.rm = TRUE),
      sd_p_hat   = sd(p_hats, na.rm = TRUE),
      p_bias     = p_bias,
      p_rmse     = p_rmse,
      
      mean_K_mean  = mean(K_mean_vec, na.rm = TRUE),
      mean_K_prop1 = mean(K_prop1_vec, na.rm = TRUE),
      
      mean_sampled_inside   = mean(n_sampled_inside_vec, na.rm = TRUE),
      mean_sampled_outside  = mean(n_sampled_outside_vec, na.rm = TRUE),
      mean_detected_inside  = mean(n_detected_inside_vec, na.rm = TRUE),
      mean_detected_outside = mean(n_detected_outside_vec, na.rm = TRUE),
      
      n_failed  = n_failed,
      fail_rate = fail_rate,
      top_error = top_error,
      
  
      scenario_name = scenario_name,
      cluster_size  = cluster_size
    )
    
    idx <- idx + 1
  }
}

scenario_perf <- dplyr::bind_rows(results) %>%
  mutate(
    beta_se_mean = sd_beta_hat / sqrt(n_reps),
    beta_ci_low  = mean_beta_hat - 1.96 * beta_se_mean,
    beta_ci_high = mean_beta_hat + 1.96 * beta_se_mean,
    
    p_se_mean = sd_p_hat / sqrt(n_reps),
    p_ci_low  = mean_p_hat - 1.96 * p_se_mean,
    p_ci_high = mean_p_hat + 1.96 * p_se_mean
  )

print(head(scenario_perf))


#----------- Visualizaation .--------------


ggplot(scenario_perf,
       aes(x = p_detect, y = mean_beta_hat, color = factor(total_visits))) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = beta_ci_low, ymax = beta_ci_high), width = 0.015) +
  geom_line() +
  geom_hline(yintercept = beta_corridor_true, linetype = "dashed") +
  facet_wrap(~ total_visits, nrow = 1) +
  labs(title = paste(species_name, "- Corridor Effect Recovery by Sample Size"),
       subtitle = scenario_name,
       x = "True detection probability (p_detect)",
       y = "Mean estimated corridor effect (log-odds)",
       color = "Total visits") +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2")


#with consistent y-axis - what range to use??? -0.5 to 0.5? 
ggplot(scenario_perf,
       aes(x = p_detect, y = mean_beta_hat, color = factor(total_visits))) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = beta_ci_low, ymax = beta_ci_high), width = 0.015) +
  geom_line() +
  geom_hline(yintercept = beta_corridor_true, linetype = "dashed") +
  facet_wrap(~ total_visits, nrow = 1, scales = "fixed") +
  labs(title = paste(species_name, "- Corridor Effect Recovery by Sample Size"),
       subtitle = scenario_name,
       x = "True detection probability (p_detect)",
       y = "Mean estimated corridor effect (log-odds)",
       color = "Total visits") +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(-0.5, 0.5))


ggplot(scenario_perf,
       aes(x = p_detect, y = mean_p_hat, color = factor(total_visits), group = total_visits)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = p_ci_low, ymax = p_ci_high), width = 0.015) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  facet_wrap(~ n_samples, nrow = 1) +
  labs(
    title = paste(species_name, "- Detection Probability Recovery by Sample Size"),
    subtitle = scenario_name,
    x = "True detection probability (p_detect)",
    y = "Mean estimated detection probability (p_hat)",
    color = "Total visits"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0, 1))

ggplot(scenario_perf,
       aes(x = p_detect, y = mean_p_hat, color = factor(total_visits))) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = p_ci_low, ymax = p_ci_high), width = 0.01) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(
    title = paste(species_name, "- Detection Probability Calibration"),
    subtitle = scenario_name,
    x = "True detection probability",
    y = "Estimated detection probability (mean p_hat)",
    color = "Total visits"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))

scenario_perf %>%
  dplyr::mutate(p_error = mean_p_hat - p_detect) %>%
  ggplot(aes(x = p_detect, y = p_error, color = factor(total_visits), group = total_visits)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~ n_samples, nrow = 1) +
  labs(
    title = paste(species_name, "- Detection probability error (p_hat - p_true)"),
    subtitle = scenario_name,
    x = "True detection probability",
    y = "Mean error in estimated detection probability",
    color = "Total visits"
  ) +
  theme_minimal()

ggplot(scenario_perf,
       aes(x = p_detect, y = mean_p_hat, color = factor(total_visits))) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(0.58, 1.01), ylim = c(0.58, 1.01)) +
  theme_minimal()

#save
stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

saveRDS(
  scenario_perf,
  sprintf("saved/occupancy/clustered/%s_%s_scenario_summary_%s.rds",
          scenario_ID, species_name, stamp)
)

#---------PLOT VISITS AND DETECTIONS -----------------

library(dplyr)
library(tidyr)

scenario_perf2 <- scenario_perf %>%
  dplyr::mutate(
    visits_inside  = mean_sampled_inside,
    visits_outside = mean_sampled_outside,
    det_inside     = mean_detected_inside,
    det_outside    = mean_detected_outside,
    det_rate_inside  = det_inside / visits_inside,
    det_rate_outside = det_outside / visits_outside,
    visit_share_inside = visits_inside / (visits_inside + visits_outside),
    det_share_inside   = det_inside / (det_inside + det_outside),
    det_rate_ratio_in_out = det_rate_inside / det_rate_outside
  )


library(ggplot2)
library(tidyr)

#Visits

visits_long <- scenario_perf2 %>%
  tibble::as_tibble() %>%                                   # ensure tibble
  dplyr::select(p_detect, n_samples, visits_inside, visits_outside) %>%
  tidyr::pivot_longer(cols = c(visits_inside, visits_outside),
                      names_to = "where", values_to = "visits") %>%
  dplyr::mutate(where = ifelse(where == "visits_inside",
                               "Inside corridor", "Outside corridor"))
ggplot(visits_long,
       aes(x = factor(p_detect), y = visits, fill = where)) +
  geom_col(position = "stack") +
  facet_wrap(~ n_samples, nrow = 1) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(title = "Visits (sampled grid cells): Inside vs Outside corridor",
       x = "Detection probability", y = "Mean visits per scenario") +
  theme_minimal()


#Detections 

det_long <- scenario_perf2 %>%
  tibble::as_tibble() %>%
  dplyr::select(p_detect, n_samples, det_inside, det_outside) %>%
  tidyr::pivot_longer(cols = c(det_inside, det_outside),
                      names_to = "where", values_to = "detections") %>%
  dplyr::mutate(where = dplyr::recode(where,
                                      "det_inside" = "Inside corridor",
                                      "det_outside" = "Outside corridor"))

ggplot2::ggplot(det_long,
                ggplot2::aes(x = factor(p_detect),
                             y = detections,
                             fill = where)) +
  ggplot2::geom_col(position = "stack") +
  ggplot2::facet_wrap(~ n_samples, nrow = 1) +
  ggplot2::scale_fill_brewer(palette = "Set2", name = NULL) +
  ggplot2::labs(
    title = "Detections: Inside vs Outside the Corridor",
    x = "Detection probability",
    y = "Mean detections"
  ) +
  ggplot2::theme_minimal()

#Detection rates 
rates_long <- scenario_perf2 %>%
  tibble::as_tibble() %>%
  dplyr::select(p_detect, n_samples, det_rate_inside, det_rate_outside) %>%
  tidyr::pivot_longer(cols = c(det_rate_inside, det_rate_outside),
                      names_to = "where", values_to = "det_rate") %>%
  dplyr::mutate(where = dplyr::recode(where,
                                      "det_rate_inside" = "Inside corridor",
                                      "det_rate_outside" = "Outside corridor"))

ggplot2::ggplot(rates_long,
                ggplot2::aes(x = p_detect,
                             y = det_rate,
                             color = where)) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::facet_wrap(~ n_samples, nrow = 1) +
  ggplot2::scale_color_brewer(palette = "Dark2", name = NULL) +
  ggplot2::labs(
    title = "Detection rate (detections per visit)",
    x = "Detection probability",
    y = "Detection rate"
  ) +
  ggplot2::theme_minimal()

#Inside share (%) visits and detections 
ggplot2::ggplot(scenario_perf2,
                ggplot2::aes(x = p_detect,
                             y = visit_share_inside,
                             color = factor(n_samples))) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = "Share of visits inside the corridor",
    x = "Detection probability",
    y = "% visits inside",
    color = "Sample size"
  ) +
  ggplot2::theme_minimal()

#Detections inside share
ggplot2::ggplot(scenario_perf2,
                ggplot2::aes(x = p_detect,
                             y = det_share_inside,
                             color = factor(n_samples))) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = "Share of detections inside the corridor",
    x = "Detection probability",
    y = "% detections inside",
    color = "Sample size"
  ) +
  ggplot2::theme_minimal()

#Inside/outside detection rate ratio 
ggplot2::ggplot(scenario_perf2,
                ggplot2::aes(x = p_detect,
                             y = det_rate_ratio_in_out,
                             color = factor(n_samples))) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::labs(
    title = "Detection rate ratio (Inside / Outside)",
    x = "Detection probability",
    y = "Rate ratio",
    color = "Sample size"
  ) +
  ggplot2::theme_minimal()

#side by side grouped baars for counts (inside vs outside)

#visits
visits_long <- scenario_perf2 %>%
  tibble::as_tibble() %>%
  dplyr::select(p_detect, n_samples, visits_inside, visits_outside) %>%
  tidyr::pivot_longer(cols = c(visits_inside, visits_outside),
                      names_to = "where", values_to = "visits") %>%
  dplyr::mutate(where = dplyr::recode(where,
                                      "visits_inside" = "Inside corridor",
                                      "visits_outside" = "Outside corridor"))

ggplot2::ggplot(visits_long,
                ggplot2::aes(x = factor(p_detect),
                             y = visits,
                             fill = where)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7)) +
  ggplot2::facet_wrap(~ n_samples, nrow = 1) +
  ggplot2::scale_fill_brewer(palette = "Set2", name = NULL) +
  ggplot2::labs(
    title = "Visits: Inside vs Outside (grouped bars)",
    x = "Detection probability",
    y = "Mean visits"
  ) +
  ggplot2::theme_minimal()

#detections 
det_long <- scenario_perf2 %>%
  tibble::as_tibble() %>%
  dplyr::select(p_detect, n_samples, det_inside, det_outside) %>%
  tidyr::pivot_longer(cols = c(det_inside, det_outside),
                      names_to = "where", values_to = "detections") %>%
  dplyr::mutate(where = dplyr::recode(where,
                                      "det_inside" = "Inside corridor",
                                      "det_outside" = "Outside corridor"))

ggplot2::ggplot(det_long,
                ggplot2::aes(x = factor(p_detect),
                             y = detections,
                             fill = where)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7)) +
  ggplot2::facet_wrap(~ n_samples, nrow = 1) +
  ggplot2::scale_fill_brewer(palette = "Set2", name = NULL) +
  ggplot2::labs(
    title = "Detections: Inside vs Outside (grouped bars)",
    x = "Detection probability",
    y = "Mean detections"
  ) +
  ggplot2::theme_minimal()


# % of visits that result in a detection - inside vs outside corridor

rates_pct_long <- scenario_perf2 %>%
  dplyr::select(p_detect, n_samples, det_rate_inside, det_rate_outside) %>%
  tidyr::pivot_longer(
    cols = c(det_rate_inside, det_rate_outside),
    names_to = "where",
    values_to = "rate"
  ) %>%
  dplyr::mutate(
    where = dplyr::recode(where,
                          "det_rate_inside" = "Inside corridor",
                          "det_rate_outside" = "Outside corridor"),
    rate_pct = rate * 100
  )

ggplot(rates_pct_long,
       aes(x = p_detect,
           y = rate_pct,
           color = where)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1) +
  facet_wrap(~ n_samples, nrow = 1) +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  labs(
    title = "% of visits that result in a detection",
    x = "Detection probability",
    y = "Detection rate (%)"
  ) +
  theme_minimal()


#Bar plot
ggplot(rates_pct_long,
       aes(x = factor(p_detect),
           y = rate_pct,
           fill = where)) +
  geom_col(position = position_dodge(width = 0.7)) +
  facet_wrap(~ n_samples, nrow = 1) +
  scale_fill_brewer(palette = "Dark2", name = NULL) +
  labs(
    title = "% of visits that result in a detection (grouped bars)",
    x = "Detection probability",
    y = "Detection rate (%)"
  ) +
  theme_minimal()

head(scenario_perf)

#plot fail rates:
library(scales)


scenario_perf2 <- scenario_perf %>%
  mutate(
    species = factor(species),
    scenario_ID = factor(scenario_ID),
    p_detect = as.numeric(p_detect),
    n_samples = as.integer(n_samples),
    total_visits = as.integer(total_visits),
    n_failed = as.integer(n_failed),
    fail_rate = as.numeric(fail_rate)
  )

# ---- 1) Fail rate plot ----
p_fail_rate <- ggplot(scenario_perf2,
                      aes(x = total_visits, y = fail_rate, color = factor(p_detect))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(species ~ scenario_ID, scales = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Model fail rate by sampling effort",
    x = "total_visits",
    y = "Fail rate",
    color = "p_detect"
  ) +
  theme_minimal(base_size = 12)

# ---- 2) Number failed plot ----
p_n_failed <- ggplot(scenario_perf2,
                     aes(x = total_visits, y = n_failed, color = factor(p_detect))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(species ~ scenario_ID, scales = "free_x") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title = "Number of failed fits by sampling effort",
    x = "total_visits",
    y = "n_failed",
    color = "p_detect"
  ) +
  theme_minimal(base_size = 12)

p_fail_rate
p_n_failed

