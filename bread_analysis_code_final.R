# ==============================================================================
# Analysis code: does a glass of water in the microwave reduce moisture
# loss when defrosting sourdough bread?
#
# Design: small-scale randomized experiment, 14 blocks of two frozen sourdough
# slices (same-numbered slice from pack A and pack B), one treated (microwaved
# with a glass of 140 g of room-temperature water) and one control (no glass),
# 35 seconds at default power. Outcome: percentage weight loss.
#
# This script does three things, in this order:
#   1. Randomization (coin flips -> assignment)
#   2. Cleaning and derived variables
#   3. Main analysis (estimate, randomization inference, secondary models, figures)
#
# FINAL VERSION, run on the real data. It is the registered script
# (bread_analysis_code.R) with these changes made after data collection:
#   - the real coin flips replace the placeholders in section 2,
#   - bread_data.csv holds the measured weights (columns slice_id, pre_weight,
#     post_weight; slice_id = pack letter + two-digit block, e.g. A07),
#   - section 8 adds a one-row summary table of the main result (results_table),
#   - Figure 1 (section 10a) is redrawn with 0.5-point bars and a whole-number
#     y-axis, and is saved to a file,
#   - section 11 adds a Coppock-style figure for the write-up, which was NOT in
#     the pre-analysis plan.
# The estimates, tests and models are unchanged.
# ==============================================================================

# ==========================================================================
# 0. PACKAGES
# ==========================================================================
# Install any missing packages (only needed the first time)
pkgs <- c("tidyverse", "stargazer")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

# 'ri' is no longer on CRAN, so it is installed from the read-only CRAN mirror
# on GitHub.
if (!"ri" %in% rownames(installed.packages())) {
  install.packages("https://codeload.github.com/cran/ri/tar.gz/refs/heads/master",
                   repos = NULL, type = "source")
}

library(tidyverse)   # dplyr, tidyr, ggplot2, readr, tibble
library(ri)          # randomization inference
library(stargazer)   # regression tables


# ==========================================================================
# 1. FIXED PARAMETERS
# ==========================================================================
n_blocks     <- 14
T_sec        <- 35          # microwave time (seconds), same in both arms
power        <- "default"   # microwave setting, same in both arms
post_delay_s <- 10          # post-weight taken this many seconds after the microwave stops


# ==========================================================================
# 2. RANDOMIZATION (coin flips -> assignment)
# ==========================================================================
# PHYSICAL COIN FLIPS, in block order (H = heads, T = tails).
# (Recovered from the recorded treatment and order of each block.)
# flip 1: T means the pack A slice is treated
# flip 2: T means the treated slice is microwaved first
flip1_raw <- c("T","T","H","T","H","H","T","H","T","H","H","T","T","T")
flip2_raw <- c("T","H","T","H","T","H","T","H","H","H","H","T","H","T")
using_placeholder_flips <- FALSE  # real flips entered

stopifnot(length(flip1_raw) == n_blocks, length(flip2_raw) == n_blocks,
          all(c(flip1_raw, flip2_raw) %in% c("H", "T")))

blocks <- tibble(
  block = 1:n_blocks,
  flip1 = ifelse(flip1_raw == "T", "tails", "heads"),
  flip2 = ifelse(flip2_raw == "T", "tails", "heads")
)
using_placeholder_flips

assignment <- blocks %>%
  crossing(pack = c("A", "B")) %>%
  mutate(
    treat = as.integer((flip1 == "tails" & pack == "A") |
                       (flip1 == "heads" & pack == "B")),
    treated_first = as.integer(flip2 == "tails"),
    order = case_when(
      treat == 1 & treated_first == 1 ~ 1L,  # treated slice goes first
      treat == 0 & treated_first == 0 ~ 1L,  # control slice goes first
      TRUE                            ~ 2L
    ),
    slice_id = paste0(pack, sprintf("%02d", block))
  ) %>%
  arrange(block, pack)

# Checks that the procedure was followed
stopifnot(
  nrow(assignment) == 2 * n_blocks,
  all(tapply(assignment$treat, assignment$block, sum) == 1),   # one treated per block
  all(tapply(assignment$order, assignment$block, sum) == 3)    # one first, one second
)
table(assignment$treat)
table(assignment$pack, assignment$treat)
table(assignment$treated_first[assignment$treat == 1])


# ==========================================================================
# 3. SPREADSHEETS (optional; run once if wanted)
# ==========================================================================
## write_csv(assignment, "assignment_researcher.csv")
## 
## assignment %>%
##   arrange(block, order) %>%
##   transmute(block, order, slice_id, pack,
##             pre_weight = NA_real_, post_weight = NA_real_) %>%
##   write_csv("sheet_for_RA.csv")


# ==========================================================================
# 4. DATA
# ==========================================================================
if (file.exists("bread_data.csv") && using_placeholder_flips) {
  stop("bread_data.csv exists but the flips in section 2 are still placeholders. Enter the real flips and set using_placeholder_flips <- FALSE.")
}

if (file.exists("bread_data.csv")) {
  using_simulated <- FALSE
  raw <- read_csv("bread_data.csv") %>% left_join(assignment, by = "slice_id")
} else {
  using_simulated <- TRUE
  set.seed(2026)
  block_size <- rnorm(n_blocks, mean = 42, sd = 5)
  block_eff  <- rnorm(n_blocks, mean = 0, sd = 0.7)
  raw <- assignment %>%
    mutate(
      pre_weight  = round(block_size[block] + rnorm(n(), 0, 0.8), 1),
      true_loss   = 3 + block_eff[block] - 0.8 * treat + rnorm(n(), 0, 0.5),
      post_weight = round(pre_weight * (1 - true_loss / 100), 1)
    ) %>%
    select(-true_loss)
}
using_simulated


# ==========================================================================
# 5. CLEANING AND DERIVED VARIABLES
# ==========================================================================
dat <- raw %>%
  mutate(
    loss_pct  = (pre_weight - post_weight) / pre_weight * 100,
    condition = factor(ifelse(treat == 1, "glass", "no glass"),
                       levels = c("no glass", "glass"))
  )

# Block-level data: one row per block, treated minus control
bd <- dat %>%
  group_by(block) %>%
  summarise(
    d_pct         = loss_pct[treat == 1] - loss_pct[treat == 0],
    pre_diff      = pre_weight[treat == 1] - pre_weight[treat == 0],
    treated_first = order[treat == 1] == 1,
    .groups = "drop"
  )
n <- nrow(bd)
n


# ==========================================================================
# 6. RANDOMIZATION CHECKS (randomization inference)
# ==========================================================================
Z        <- dat$treat
Y        <- dat$loss_pct
blockvar <- dat$block

set.seed(1234567)                                  # so the sampled assignments are reproducible
prob  <- genprobexact(Z, blockvar)                 # probability of treatment (0.5 in every block)
perms <- genperms(Z, blockvar, maxiter = 10000)    # possible assignments, within blocks
dim(perms)                                         # rows = slices, columns = assignments

# F-statistic in the realized data and under every sampled assignment
Fstat  <- summary(lm(treat ~ pre_weight, data = dat))$fstatistic[1]
Fstore <- apply(perms, 2, function(z) summary(lm(z ~ pre_weight, data = dat))$fstatistic[1])
mean(Fstore >= Fstat)   # randomization inference p-value

# Descriptive balance: pre-weight by arm, order split, and pack split
dat %>% group_by(condition) %>% summarise(mean_pre_weight = mean(pre_weight))
table(bd$treated_first)
table(dat$pack, dat$treat)


# ==========================================================================
# 7. MAIN ANALYSIS: ESTIMATE
# ==========================================================================
ate <- mean(bd$d_pct)   # mean of the 14 within-block differences
ate


# ==========================================================================
# 8. RANDOMIZATION INFERENCE (sharp null, two-sided)
# ==========================================================================
ate_ri <- estate(Y, Z, prob = prob)               # difference in means

# Sharp null of no effect: potential outcomes are the same under both conditions
Ys_null   <- genouts(Y, Z, ate = 0)
dist_null <- gendist(Ys_null, perms, prob = prob, HT = FALSE)
p_two_sided <- mean(abs(dist_null) >= abs(ate_ri) - 1e-12)

# Confidence interval: distribution of the estimates if the true effect equals the estimate
Ys_est   <- genouts(Y, Z, ate_ri)
dist_est <- gendist(Ys_est, perms, prob = prob, HT = FALSE)
disp     <- dispdist(dist_est, ate_ri, display.plot = FALSE)
ci_ri    <- disp$quantile   # 95% confidence interval
se_ri    <- disp$sd         # standard error of the randomization distribution

c(ATE = ate_ri, SE = se_ri, p_two_sided = p_two_sided)
ci_ri

# Summary of the main result (randomization inference)
results_table <- tibble(
  ATE         = as.numeric(ate_ri),
  SE          = se_ri,
  CI_low      = unname(ci_ri[1]),
  CI_high     = unname(ci_ri[2]),
  p_two_sided = p_two_sided
) %>% mutate(across(everything(), ~ round(.x, 4)))
results_table


# ==========================================================================
# 9. SECONDARY MODELS
# ==========================================================================
m1 <- lm(loss_pct ~ treat, data = dat)                              # naive difference in means
m2 <- lm(loss_pct ~ treat + factor(block), data = dat)              # block fixed effects (primary)
m3 <- lm(loss_pct ~ treat + factor(block) + pre_weight, data = dat) # + pre-weight

stargazer(m1, m2, m3, header = FALSE, type = "text",
          keep = "treat", covariate.labels = "Glass of water",
          dep.var.labels = "Weight loss (%)",
          title = "Effect of a glass of water on weight loss when microwaving frozen sourdough",
          omit.stat = c("f", "ser"),
          add.lines = list(c("Block fixed effects", "No", "Yes", "Yes")))


# ==========================================================================
# 10a. FIGURE: distribution of outcomes by condition
# ==========================================================================
fig1 <- ggplot(dat, aes(loss_pct)) +
  geom_histogram(binwidth = 0.5, boundary = 0, closed = "left") +   # bars of 0.5 percentage points
  facet_wrap(~ condition, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = 0:7, limits = c(0, 7.5)) +
  scale_y_continuous(breaks = function(lims) seq(0, floor(lims[2]), by = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  labs(x = "Weight loss (% of initial weight)", y = "Number of slices")
fig1
ggsave("figure1_distribution_by_condition.png", fig1, width = 5, height = 4, dpi = 300)


# ==========================================================================
# 10b. FIGURE: distribution of estimates under the sharp null
# ==========================================================================
ggplot(tibble(est = as.numeric(dist_null)), aes(est)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = ate_ri, colour = "red", linetype = "dashed") +
  theme_bw() +
  labs(x = "Estimated ATE under the sharp null", y = "Assignments")


# ==========================================================================
# 11. FIGURE FOR THE WRITE-UP: Coppock-style graph (NOT in the pre-analysis plan)
# ==========================================================================
# Individual slices (points), group means with 95% confidence intervals, and
# lines joining the two slices of each block, so the design (blocking) is visible.
summary_loss <- dat %>%
  group_by(condition) %>%
  summarise(mean_loss = mean(loss_pct),
            se = sd(loss_pct) / sqrt(n()),
            lo = mean_loss - qt(0.975, n() - 1) * se,
            hi = mean_loss + qt(0.975, n() - 1) * se,
            .groups = "drop")
summary_loss

fig2 <- ggplot(summary_loss, aes(condition, mean_loss)) +
  geom_line(data = dat, aes(y = loss_pct, group = block), alpha = 0.25) +
  geom_point(data = dat, aes(y = loss_pct), alpha = 0.35,
             position = position_jitter(width = 0.03, height = 0)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ylab("Weight loss (% of initial weight)") +
  ggtitle("Weight loss by condition")
fig2
ggsave("figure2_coppock_style.png", fig2, width = 5, height = 4, dpi = 300)
