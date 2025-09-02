library(dplyr)
library(tidyr)
library(mirt)
library(ggplot2)
library(readr)
library(tibble)
library(scales)
library(ggplot2)
library(patchwork)

sem_imv_irw <- function(file,
                        model_type = "2PL",
                        sem_type   = c("test","person"),
                        bin_width  = 0.5,
                        baseline   = c("global","bin"),
                        trim_b_abs = 6) {
  
  sem_type <- match.arg(sem_type)
  baseline <- match.arg(baseline)
  
  # 1) Read IRW, keep only needed cols
  raw <- read_csv(file, show_col_types = FALSE) %>%
    select(id, item, resp) %>%
    filter(!is.na(resp))
  
  # 2) Wide for mirt
  wide <- raw %>%
    pivot_wider(names_from = item, values_from = resp) %>%
    arrange(id)
  resp_mat <- wide %>% select(-id)
  
  # 3) Fit IRT
  mod <- mirt(resp_mat, 1, itemtype = model_type, SE = TRUE)
  
  # 4) Item params
  it <- as.data.frame(coef(mod, IRTpars = TRUE, simplify = TRUE)$items)
  a_hat <- it[,"a"]; b_hat <- it[,"b"]
  
  # 5) Person scores (theta, person-level SE)
  th <- fscores(mod, full.scores.SE = TRUE)
  colnames(th) <- c("theta", "sem_person")
  th_df <- as_tibble(th) %>% mutate(id = wide$id)
  
  # 6) Join back to long
  long <- raw %>%
    left_join(th_df, by = "id") %>%
    left_join(it %>% rownames_to_column("item") %>% select(item, a = a, b = b),
              by = "item") %>%
    mutate(
      theta = as.numeric(theta),
      a = as.numeric(a),
      b = as.numeric(b),
      p_full = plogis(a * (theta - b)),
      ll_full = -(resp * log(p_full + 1e-12) + (1 - resp) * log(1 - p_full + 1e-12))
    )
  
  # 7) Baseline (per-item)
  if (baseline == "global") {
    item_base <- long %>% group_by(item) %>% summarise(p_base = mean(resp, na.rm=TRUE), .groups="drop")
    long <- long %>% left_join(item_base, by="item")
  } else {
    # bin baseline: mean within theta bin later (computed per bin below)
  }
  
  # 8) Binning
  theta_min <- floor(min(long$theta, na.rm = TRUE))
  theta_max <- ceiling(max(long$theta, na.rm = TRUE))
  breaks <- seq(theta_min, theta_max, by = bin_width)
  long <- long %>% mutate(theta_bin = cut(theta, breaks = breaks, include.lowest = TRUE))
  
  # 9) Test-level SEM curve on θ-grid, then map to bins
  theta_grid <- seq(theta_min, theta_max, length.out = 201)
  Pg_grid <- sapply(seq_along(a_hat), function(j) plogis(a_hat[j] * (theta_grid - b_hat[j])))
  Ij_grid <- sapply(seq_along(a_hat), function(j) a_hat[j]^2 * Pg_grid[,j] * (1 - Pg_grid[,j]))
  Itest_grid <- rowSums(Ij_grid)
  SEM_test_grid <- 1 / sqrt(Itest_grid)
  
  grid_bins <- cut(theta_grid, breaks = breaks, include.lowest = TRUE)
  sem_test_by_bin <- tapply(SEM_test_grid, grid_bins, mean, na.rm = TRUE)
  
  # 10) IMV per bin
  # model log-loss
  eps <- 1e-12
  clamp <- function(p) pmin(pmax(p, eps), 1 - eps)
  
  if (baseline == "global") {
    long <- long %>%
      mutate(ll_base = -(resp * log(clamp(p_base)) + (1 - resp) * log(clamp(1 - p_base))))
  }
  
  imv_df <- long %>%
    group_by(theta_bin) %>%
    {
      if (baseline == "bin") {
        mutate(., p_base_bin = mean(resp, na.rm=TRUE),
               ll_base = -(resp * log(clamp(p_base_bin)) + (1 - resp) * log(clamp(1 - p_base_bin))))
      } else . 
    } %>%
    summarise(
      ll_full = mean(ll_full, na.rm = TRUE),
      ll_base = mean(ll_base, na.rm = TRUE),
      IMV_abs = ll_base - ll_full,
      IMV_rel = IMV_abs / ll_base,
      .groups = "drop"
    )
  
  # 11) Assemble SEM per bin (both flavors)
  sem_person_by_bin <- long %>%
    group_by(theta_bin) %>%
    summarise(SEM_person = mean(sem_person, na.rm = TRUE), .groups = "drop")
  
  sem_df <- tibble(theta_bin = levels(long$theta_bin)) %>%
    left_join(imv_df, by = "theta_bin") %>%
    mutate(SEM_test = as.numeric(sem_test_by_bin[match(theta_bin, names(sem_test_by_bin))])) %>%
    left_join(sem_person_by_bin, by = "theta_bin")
  
  # 12) Choose which SEM to plot
  sem_col <- if (sem_type == "test") "SEM_test" else "SEM_person"
  
  # 13) Plot: IMV (absolute) vs chosen SEM across bins
  plot_lines <- sem_df %>%
    mutate(
      InvSEM = 1 / !!as.name(sem_col),
      IMV_s  = rescale(IMV_abs),
      InvSEM_s = rescale(InvSEM)
    ) %>%
    pivot_longer(cols = c(IMV_s, InvSEM_s), names_to = "metric", values_to = "value") %>%
    ggplot(aes(theta_bin, value, color = metric, group = metric)) +
    geom_line() + geom_point() +
    labs(title = paste0("IMV vs 1/SEM across θ (", sem_type, "-level SEM, ", baseline, " baseline)"),
         x = "θ bin", y = "scaled value", color = NULL) +
    theme_minimal()
  
  plot_scatter <- sem_df %>%
    mutate(SEM = !!as.name(sem_col)) %>%
    ggplot(aes(x = SEM, y = IMV_abs)) +
    geom_point() + geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste0("IMV (abs) vs SEM per bin (", sem_type, "-level)"),
         x = "SEM", y = "IMV (log-loss improvement)") +
    theme_minimal()
  
  p_theta <- th_df %>%
    ggplot(aes(theta)) +
    geom_histogram(bins = 40, color = "black") +
    geom_vline(aes(xintercept = mean(theta, na.rm = TRUE)), linetype = 2) +
    labs(title = "Distribution of θ (EAP)", x = expression(theta), y = "Count") +
    theme_minimal()
  
  # (B) a distribution
  p_a <- it %>%
    ggplot(aes(x = a)) +
    geom_histogram(bins = 40, color = "black") +
    labs(title = "Distribution of item discrimination (a)", x = "a", y = "Count") +
    theme_minimal()
  
  # (C) b distribution (optionally trim extreme values for readability)
  it_plot <- it %>%
    mutate(b_plot = ifelse(abs(b) > trim_b_abs, NA, b))
  p_b <- it_plot %>%
    ggplot(aes(x = b_plot)) +
    geom_histogram(bins = 40, color = "black", na.rm = TRUE) +
    labs(title = paste0("Distribution of item difficulty (b)",
                        ifelse(is.finite(trim_b_abs),
                               paste0(" | |b| ≤ ", trim_b_abs),
                               "")),
         x = "b", y = "Count") +
    theme_minimal()
  
  list(
    sem_type   = sem_type,
    baseline   = baseline,
    summary    = sem_df,
    plot_lines = plot_lines,
    plot_scatter = plot_scatter,
    model      = mod,
    p_theta = p_theta,
    p_a = p_a,
    p_b = p_b
  )
}



files <- c(
  "koreannursing_park_2017.csv","ptam_kretzschmar_2017_reasoning.csv","icar_sapa.csv",
  "geiser_tam.csv","chess_lnirt.csv","roar_lexical.csv",
  "hmcdm_spatialreasoning.csv","gcbs_brotherton_2013_vcl.csv","art.csv"
)

sem_choice <- "test"
base_choice <- "bin"

results <- lapply(files, function(f) {
  sem_imv_irw(
    file      = f,
    sem_type  = sem_choice,
    baseline  = base_choice,
    bin_width = 0.5,
    model_type = "2PL"
  )
})

## ---- build a 3x3 grid for the line plots ----
line_plots <- Map(function(res, f) {
  res$plot_lines + ggtitle(paste0("SEM vs IMV — ", basename(f)))
}, results, files)

grid_lines <- wrap_plots(plotlist = line_plots, ncol = 3) +
  plot_annotation(
    title = paste0("IMV vs 1/SEM across θ (", sem_choice,
                   "-level SEM, ", base_choice, " baseline)"),
    theme = theme(plot.title = element_text(hjust = 0.5))
  )
print(grid_lines)

ggsave("sem_imv_lines_grid_3x3.png", grid_lines, width = 16, height = 12, dpi = 300)

## ---- build a 3x3 grid for the scatter plots ----
scatter_plots <- Map(function(res, f) {
  res$plot_scatter + ggtitle(paste0("IMV vs SEM — ", basename(f)))
}, results, files)

grid_scatter <- wrap_plots(plotlist = scatter_plots, ncol = 3) +
  plot_annotation(
    title = paste0("IMV (abs) vs SEM per bin (", sem_choice,
                   "-level SEM, ", base_choice, " baseline)"),
    theme = theme(plot.title = element_text(hjust = 0.5))
  )
print(grid_scatter)
