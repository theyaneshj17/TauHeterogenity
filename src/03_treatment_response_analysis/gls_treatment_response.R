################################################################################
# GLS with Natural Cubic Splines + CAR(1)
# Following Carolyn et al.'s longitudinal modeling approach
# Author: theyaneshj17
# Last updated: 2025
################################################################################

rm(list = ls())

# Set working directory - update this path to your local repo
# setwd("path/to/your/TauHeterogenity/")

library(dplyr)
library(tidyr)
library(nlme)
library(splines)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(cowplot)
library(boot)

################################################################################
# 0. DATA PREPARATION
################################################################################

# Load data - files should be placed in data/ folder at repo root
pacc     <- read.csv("../../data/PACC.csv")
subjinfo <- read.csv("../../data/SUBJINFO.csv")
clusters <- read.csv("../../data/Brainplotk=2.csv")

# Standardize subject IDs to uppercase
pacc$BID           <- toupper(pacc$BID)
subjinfo$BID       <- toupper(subjinfo$BID)
clusters$SubjectID <- toupper(clusters$SubjectID)

# Merge datasets
subjinfo_clean <- subjinfo %>% select(BID, TX)
clusters_clean <- clusters %>% select(SubjectID, Subtype)

df <- pacc %>%
  left_join(subjinfo_clean, by = "BID") %>%
  left_join(clusters_clean, by = c("BID" = "SubjectID"))

# Keep complete cases
df_complete <- df %>%
  filter(!is.na(Subtype) & !is.na(TX) & !is.na(PACC))

# Create model variables
df_model <- df_complete %>%
  mutate(
    Subtype   = factor(Subtype, levels = c("S1", "S2")),
    Treatment = factor(ifelse(TX == "Placebo", "Placebo", "Solanezumab"),
                       levels = c("Placebo", "Solanezumab")),
    T_time    = as.numeric(as.character(VISCODE)),
    person_id = as.factor(BID),
    Group     = interaction(Subtype, Treatment)
  ) %>%
  filter(T_time <= 108)

# Extract baseline PACC (visit 1)
baseline_pacc <- df_model %>%
  filter(T_time == 1) %>%
  select(BID, PACC_baseline = PACC)

################################################################################
# 1. PREPARE FULL COVARIATE DATASET
################################################################################

# Extract covariates from subject info
subjinfo_full <- subjinfo %>%
  mutate(BID = toupper(BID)) %>%
  select(BID, AGEYR, SEX, EDCCNTU, RACE, SUVRCER, AMYLCENT, APOEGN, APOEGNPRSNFLG)

# Merge all covariates
df_gls <- df_model %>%
  left_join(subjinfo_full, by = "BID") %>%
  left_join(baseline_pacc, by = "BID") %>%
  mutate(
    T_weeks   = T_time * 4.33,         # Convert visit months to weeks
    SEX       = factor(SEX),
    RACE      = factor(RACE),
    APOE4     = factor(APOEGNPRSNFLG), # APOE-e4 carrier status (0 = non-carrier, 1 = carrier)
    person_id = factor(person_id)
  ) %>%
  filter(!is.na(PACC_baseline) & !is.na(AGEYR) & !is.na(SUVRCER) & !is.na(APOE4))

cat("Sample size after filtering:\n")
cat("  Observations:", nrow(df_gls), "\n")
cat("  Patients:",     length(unique(df_gls$person_id)), "\n")

################################################################################
# 2. FIT GLS MODEL (Carolyn's approach)
################################################################################
# Model specification:
#   - Fixed effects: natural cubic splines (df=3) x Group interaction
#   - Covariates: age, sex, education, race, APOE-e4, baseline amyloid, baseline PACC
#   - Correlation: CAR(1) within subjects
#   - Variance: heteroscedastic by time point

# Full model (with heteroscedastic variance)
gls_model <- gls(
  PACC ~ ns(T_weeks, df = 3) * Group +
    AGEYR + SEX + EDCCNTU + RACE + APOE4 + SUVRCER + PACC_baseline,
  data        = df_gls,
  correlation = corCAR1(form = ~ T_weeks | person_id),
  weights     = varIdent(form = ~ 1 | T_weeks),
  method      = "REML"
)
summary(gls_model)

# Simplified model (CAR(1) only, no heteroscedasticity)
gls_model_simple <- gls(
  PACC ~ ns(T_weeks, df = 3) * Group +
    AGEYR + SEX + EDCCNTU + RACE + APOE4 + SUVRCER + PACC_baseline,
  data        = df_gls,
  correlation = corCAR1(form = ~ T_weeks | person_id),
  method      = "REML"
)
summary(gls_model_simple)

################################################################################
# 3. SPLINE DEGREE SELECTION (AIC comparison)
################################################################################

gls_df2 <- gls(
  PACC ~ ns(T_weeks, df = 2) * Group +
    AGEYR + SEX + EDCCNTU + RACE + APOE4 + SUVRCER + PACC_baseline,
  data = df_gls, correlation = corCAR1(form = ~ T_weeks | person_id), method = "REML"
)
gls_df3 <- gls(
  PACC ~ ns(T_weeks, df = 3) * Group +
    AGEYR + SEX + EDCCNTU + RACE + APOE4 + SUVRCER + PACC_baseline,
  data = df_gls, correlation = corCAR1(form = ~ T_weeks | person_id), method = "REML"
)
gls_df4 <- gls(
  PACC ~ ns(T_weeks, df = 4) * Group +
    AGEYR + SEX + EDCCNTU + RACE + APOE4 + SUVRCER + PACC_baseline,
  data = df_gls, correlation = corCAR1(form = ~ T_weeks | person_id), method = "REML"
)

cat("AIC comparison:\n")
cat("  df = 2:", AIC(gls_df2), "\n")
cat("  df = 3:", AIC(gls_df3), "\n")
cat("  df = 4:", AIC(gls_df4), "\n")

################################################################################
# 4. HELPER: Build prediction grid
################################################################################

make_pred_grid <- function(df, n_points = 100) {
  expand.grid(
    T_weeks       = seq(min(df$T_weeks), max(df$T_weeks), length.out = n_points),
    Group         = levels(df$Group),
    AGEYR         = mean(df$AGEYR),
    SEX           = factor("1", levels = levels(df$SEX)),
    EDCCNTU       = mean(df$EDCCNTU),
    RACE          = factor(names(which.max(table(df$RACE))), levels = levels(df$RACE)),
    APOE4         = factor("0", levels = levels(df$APOE4)),
    SUVRCER       = mean(df$SUVRCER),
    PACC_baseline = mean(df$PACC_baseline)
  )
}

add_labels <- function(df) {
  df %>% mutate(
    Subtype       = factor(ifelse(grepl("S1", Group), "S1", "S2"), levels = c("S1", "S2")),
    Treatment     = factor(ifelse(grepl("Placebo", Group), "Placebo", "Solanezumab"),
                           levels = c("Placebo", "Solanezumab")),
    T_months      = T_weeks / 4.33,
    Subtype_label = ifelse(Subtype == "S1", "S1 (Typical)", "S2 (Cortical)")
  )
}

# Add T_months to observed data
df_gls <- df_gls %>%
  mutate(T_months = T_weeks / 4.33,
         Subtype_label = ifelse(Subtype == "S1", "S1 (Typical)", "S2 (Cortical)"))

################################################################################
# 5. VISUALIZATION: df=3, simple model (no CI)
################################################################################

pred_data_gls <- make_pred_grid(df_gls) %>%
  mutate(fit = predict(gls_model_simple, newdata = .)) %>%
  add_labels()

p_gls <- ggplot() +
  geom_point(data = df_gls,
             aes(x = T_months, y = PACC, color = Treatment),
             alpha = 0.15, size = 0.8) +
  geom_line(data = pred_data_gls,
            aes(x = T_months, y = fit, color = Treatment),
            linewidth = 1.2) +
  facet_wrap(~ Subtype_label, ncol = 2) +
  labs(x = "Time (Months)", y = "PACC Score",
       title    = "GLS with Natural Cubic Splines (df=3)",
       subtitle = "CAR(1) correlation; Covariates: age, sex, education, race, APOE-e4, amyloid, baseline PACC",
       color    = "Treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position    = "bottom",
        plot.title         = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle      = element_text(hjust = 0.5, size = 10)) +
  scale_color_manual(values = c("Placebo" = "#E74C3C", "Solanezumab" = "#3498DB"))

print(p_gls)
ggsave("../../outputs/figures/PACC_gls_natural_spline.png", p_gls, width = 10, height = 5, dpi = 300)

################################################################################
# 6. VISUALIZATION: df=4 with 95% CI band
################################################################################

pred_data_df4 <- make_pred_grid(df_gls) %>%
  mutate(fit = predict(gls_df4, newdata = .))

X4  <- model.matrix(~ ns(T_weeks, df = 4) * Group + AGEYR + SEX + EDCCNTU + RACE + APOE4 + SUVRCER + PACC_baseline,
                    data = pred_data_df4)
se4 <- sqrt(diag(X4 %*% vcov(gls_df4) %*% t(X4)))
pred_data_df4 <- pred_data_df4 %>%
  mutate(se = se4, lower = fit - 1.96 * se, upper = fit + 1.96 * se) %>%
  add_labels()

p_df4 <- ggplot() +
  geom_ribbon(data = pred_data_df4,
              aes(x = T_months, ymin = lower, ymax = upper, fill = Treatment), alpha = 0.2) +
  geom_point(data = df_gls,
             aes(x = T_months, y = PACC, color = Treatment), alpha = 0.15, size = 0.8) +
  geom_line(data = pred_data_df4,
            aes(x = T_months, y = fit, color = Treatment), linewidth = 1.2) +
  facet_wrap(~ Subtype_label, ncol = 2) +
  labs(x = "Time (Months)", y = "PACC Score",
       title    = "GLS with Natural Cubic Splines (df=4)",
       subtitle = "CAR(1) correlation; Shaded area = 95% CI",
       color = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        plot.title      = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, size = 10)) +
  scale_color_manual(values = c("Placebo" = "#E74C3C", "Solanezumab" = "#3498DB")) +
  scale_fill_manual(values  = c("Placebo" = "#E74C3C", "Solanezumab" = "#3498DB"))

print(p_df4)
ggsave("../../outputs/figures/PACC_gls_df4_CI.png", p_df4, width = 10, height = 5, dpi = 300)

################################################################################
# 7. VISUALIZATION: df=3 with 95% CI band
################################################################################

pred_data_df3 <- make_pred_grid(df_gls) %>%
  mutate(fit = predict(gls_df3, newdata = .))

X3  <- model.matrix(~ ns(T_weeks, df = 3) * Group + AGEYR + SEX + EDCCNTU + RACE + APOE4 + SUVRCER + PACC_baseline,
                    data = pred_data_df3)
se3 <- sqrt(diag(X3 %*% vcov(gls_df3) %*% t(X3)))
pred_data_df3 <- pred_data_df3 %>%
  mutate(se = se3, lower = fit - 1.96 * se, upper = fit + 1.96 * se) %>%
  add_labels()

p_df3 <- ggplot() +
  geom_ribbon(data = pred_data_df3,
              aes(x = T_months, ymin = lower, ymax = upper, fill = Treatment), alpha = 0.2) +
  geom_point(data = df_gls,
             aes(x = T_months, y = PACC, color = Treatment), alpha = 0.15, size = 0.8) +
  geom_line(data = pred_data_df3,
            aes(x = T_months, y = fit, color = Treatment), linewidth = 1.2) +
  facet_wrap(~ Subtype_label, ncol = 2) +
  labs(x = "Time (Months)", y = "PACC Score",
       title    = "GLS with Natural Cubic Splines (df=3)",
       subtitle = "CAR(1) correlation; Shaded area = 95% CI",
       color = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        plot.title      = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, size = 10)) +
  scale_color_manual(values = c("Placebo" = "#E74C3C", "Solanezumab" = "#3498DB")) +
  scale_fill_manual(values  = c("Placebo" = "#E74C3C", "Solanezumab" = "#3498DB"))

print(p_df3)
ggsave("../../outputs/figures/PACC_gls_df3_CI.png", p_df3, width = 10, height = 5, dpi = 300)

################################################################################
# 8. PREDICTION INTERVAL (df=3)
#    PI = sqrt(mean SE^2 + residual variance^2); wider than CI, covers individual observations
################################################################################

pred_data_pi <- make_pred_grid(df_gls) %>%
  mutate(fit = predict(gls_df3, newdata = .))

X_pi       <- model.matrix(~ ns(T_weeks, df = 3) * Group + AGEYR + SEX + EDCCNTU + RACE + APOE4 + SUVRCER + PACC_baseline,
                            data = pred_data_pi)
se_mean_pi <- sqrt(diag(X_pi %*% vcov(gls_df3) %*% t(X_pi)))
sigma_resid <- summary(gls_df3)$sigma

pred_data_pi <- pred_data_pi %>%
  mutate(se_pred = sqrt(se_mean_pi^2 + sigma_resid^2),
         lower   = fit - 1.96 * se_pred,
         upper   = fit + 1.96 * se_pred) %>%
  add_labels()

p_pi <- ggplot() +
  geom_ribbon(data = pred_data_pi,
              aes(x = T_months, ymin = lower, ymax = upper, fill = Treatment), alpha = 0.2) +
  geom_point(data = df_gls,
             aes(x = T_months, y = PACC, color = Treatment), alpha = 0.15, size = 0.8) +
  geom_line(data = pred_data_pi,
            aes(x = T_months, y = fit, color = Treatment), linewidth = 1.2) +
  facet_wrap(~ Subtype_label, ncol = 2) +
  labs(x = "Time (Months)", y = "PACC Score",
       title    = "GLS with Natural Cubic Splines (df=3)",
       subtitle = "CAR(1) correlation; Shaded area = 95% Prediction Interval",
       color = "Treatment", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        plot.title      = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, size = 10)) +
  scale_color_manual(values = c("Placebo" = "#E74C3C", "Solanezumab" = "#3498DB")) +
  scale_fill_manual(values  = c("Placebo" = "#E74C3C", "Solanezumab" = "#3498DB"))

print(p_pi)
ggsave("../../outputs/figures/PACC_gls_df3_PredictionInterval.png", p_pi, width = 10, height = 5, dpi = 300)

################################################################################
# 9. BOOTSTRAP CI (n = 5000, cluster bootstrap by patient)
#    Resamples whole patients to preserve within-subject correlation structure
################################################################################

set.seed(123)
n_boot          <- 5000
unique_patients <- unique(df_gls$person_id)
n_patients      <- length(unique_patients)
pred_times      <- seq(min(df_gls$T_weeks), max(df_gls$T_weeks), length.out = 50)
groups          <- levels(df_gls$Group)

# Initialize storage
boot_results <- setNames(
  lapply(groups, function(g) matrix(NA, nrow = n_boot, ncol = length(pred_times))),
  groups
)

cat("Running cluster bootstrap (n =", n_boot, ")...\n")
pb <- txtProgressBar(min = 0, max = n_boot, style = 3)

for (b in 1:n_boot) {
  # Sample patients with replacement (cluster bootstrap)
  sampled_patients <- sample(unique_patients, n_patients, replace = TRUE)

  boot_data <- do.call(rbind, lapply(seq_along(sampled_patients), function(i) {
    pd <- df_gls[df_gls$person_id == sampled_patients[i], ]
    pd$boot_id <- i
    pd
  }))
  boot_data$boot_id <- factor(boot_data$boot_id)

  tryCatch({
    boot_model <- gls(
      PACC ~ ns(T_weeks, df = 3) * Group +
        AGEYR + SEX + EDCCNTU + RACE + APOE4 + SUVRCER + PACC_baseline,
      data        = boot_data,
      correlation = corCAR1(form = ~ T_weeks | boot_id),
      method      = "REML",
      control     = glsControl(maxIter = 100, opt = "optim")
    )

    for (g in groups) {
      pred_boot <- data.frame(
        T_weeks       = pred_times,
        Group         = factor(g, levels = levels(df_gls$Group)),
        AGEYR         = mean(df_gls$AGEYR),
        SEX           = factor("1", levels = levels(df_gls$SEX)),
        EDCCNTU       = mean(df_gls$EDCCNTU),
        RACE          = factor(names(which.max(table(df_gls$RACE))), levels = levels(df_gls$RACE)),
        APOE4         = factor("0", levels = levels(df_gls$APOE4)),
        SUVRCER       = mean(df_gls$SUVRCER),
        PACC_baseline = mean(df_gls$PACC_baseline)
      )
      boot_results[[g]][b, ] <- predict(boot_model, newdata = pred_boot)
    }
  }, error = function(e) {}) # Skip non-converged iterations

  setTxtProgressBar(pb, b)
}
close(pb)

# Summarize bootstrap results
boot_summary <- do.call(rbind, lapply(groups, function(g) {
  valid      <- apply(boot_results[[g]], 1, function(x) !any(is.na(x)))
  mat        <- boot_results[[g]][valid, ]
  cat("\nGroup", g, ":", sum(valid), "successful bootstrap iterations\n")
  data.frame(
    T_weeks = pred_times,
    Group   = g,
    fit     = colMeans(mat),
    lower   = apply(mat, 2, quantile, probs = 0.025),
    upper   = apply(mat, 2, quantile, probs = 0.975)
  )
}))

boot_summary <- boot_summary %>% add_labels()

# Save bootstrap results
saveRDS(boot_summary, "../../outputs/boot_summary.rds")

################################################################################
# 10. PUBLICATION-READY PLOT: Bootstrap CI, split by subtype
################################################################################

# Publication theme
theme_Publication <- function(base_size = 14, base_family = "Helvetica") {
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title        = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      panel.background  = element_rect(colour = NA),
      plot.background   = element_rect(colour = NA),
      panel.border      = element_rect(colour = NA),
      axis.title        = element_text(face = "bold", size = rel(1)),
      axis.title.y      = element_text(angle = 90, vjust = 2),
      axis.title.x      = element_text(vjust = -1),
      axis.line         = element_line(colour = "black"),
      axis.ticks        = element_line(),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.key        = element_rect(colour = NA),
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.title      = element_blank(),
      plot.margin       = unit(c(10, 5, 10, 5), "mm"),
      strip.background  = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text        = element_text(face = "bold"),
      axis.text         = element_text(face = "bold")
    )
}

# Color palette
pal <- c("Placebo" = "#DF0069", "Solanezumab" = "#0000E3")

make_subtype_plot <- function(subtype_code, title_label) {
  ggplot() +
    geom_ribbon(data = boot_summary %>% filter(Subtype == subtype_code),
                aes(x = T_months, ymin = lower, ymax = upper, fill = Treatment),
                alpha = 0.3) +
    geom_line(data = boot_summary %>% filter(Subtype == subtype_code),
              aes(x = T_months, y = fit, color = Treatment), linewidth = 0.8) +
    geom_point(data = df_gls %>% filter(Subtype == subtype_code),
               aes(x = T_months, y = PACC, color = Treatment),
               alpha = 0.2, size = 0.8) +
    coord_cartesian(ylim = c(-25, 10)) +
    labs(x = "Time (Months)", y = "PACC Score (95% CI)", title = title_label) +
    theme_Publication(base_size = 16) +
    theme(panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white"),
          legend.position  = "none",
          axis.line        = element_line(colour = "black")) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values  = pal) +
    scale_x_continuous(breaks = seq(0, 100, by = 25))
}

p_W1 <- make_subtype_plot("S1", "W1")
p_W2 <- make_subtype_plot("S2", "W2")

# Extract shared legend
p_legend_src <- ggplot() +
  geom_line(data = boot_summary, aes(x = T_months, y = fit, color = Treatment), linewidth = 0.8) +
  scale_color_manual(values = pal) +
  theme_Publication(base_size = 16) +
  theme(legend.position = "bottom", legend.title = element_blank())
legend_grob <- cowplot::get_legend(p_legend_src)

# Combine panels
p_combined <- (p_W1 | p_W2) / legend_grob +
  plot_layout(heights = c(10, 1)) +
  plot_annotation(
    title = "Cognitive Decline by AD Subtype and Treatment",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

print(p_combined)
ggsave("../../outputs/figures/PACC_gls_bootstrap_CI_publication.pdf",
       p_combined, width = 10, height = 5)
ggsave("../../outputs/figures/PACC_gls_bootstrap_CI_publication.png",
       p_combined, width = 10, height = 5, dpi = 300)

cat("\nAll figures saved to outputs/figures/\n")
