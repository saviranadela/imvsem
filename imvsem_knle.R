library(readr)
library(dplyr)
library(tidyr)
library(mirt)
library(ggplot2)
library(tibble)


df <- read_csv("koreannursing_park_2017.csv")

df <- df %>%
  pivot_wider(names_from = item, values_from = resp)

df <- df %>%
  select(-id)

mod <- mirt(df, 1, itemtype = "2PL")

theta_scores <- fscores(mod, full.scores.SE = TRUE)
colnames(theta_scores) <- c("theta", "sem")

theta_df <- cbind(id = rownames(df), theta_scores)


## calculated predicted probabilities

item_params <- coef(mod, IRTpars = TRUE, simplify = TRUE)$items

predict_probs <- function(theta, a, b) {
  1 / (1 + exp(-a * (theta - b)))
}

theta_df <- as.data.frame(theta_df)

item_params <- as.data.frame(item_params)


long_df <- df %>%
  mutate(id = rownames(.)) %>%
  pivot_longer(-id, names_to = "item", values_to = "resp") %>%
  left_join(theta_df, by = "id") %>%
  left_join(
    item_params %>% rownames_to_column("item") %>%
      select(item, a = a, b = b),
    by = "item"
  ) %>%
  mutate(
    theta = as.numeric(theta),
    a = as.numeric(a),
    b = as.numeric(b),
    p_full = predict_probs(theta, a, b),
    logloss_full = - (resp * log(p_full) + (1 - resp) * log(1 - p_full))
  )


item_baseline <- long_df %>%
  group_by(item) %>%
  summarise(p_baseline = mean(resp, na.rm = TRUE))

long_df <- long_df %>%
  left_join(item_baseline, by = "item") %>%
  mutate(
    logloss_baseline = - (resp * log(p_baseline) + (1 - resp) * log(1 - p_baseline))
  )

## theta bins

range(long_df$theta, na.rm = TRUE)

min_theta <- floor(min(long_df$theta, na.rm = TRUE))      # -3
max_theta <- ceiling(max(long_df$theta, na.rm = TRUE))    # +4


long_df <- long_df %>%
  mutate(theta_bin = cut(theta, breaks = seq(min_theta, max_theta, by = 0.5)))

long_df <- long_df %>%
  mutate(sem = as.numeric(sem))

summary_df <- long_df %>%
  group_by(theta_bin) %>%
  summarise(
    avg_sem = mean(sem, na.rm = TRUE),
    avg_logloss_full = mean(logloss_full, na.rm = TRUE),
    avg_logloss_baseline = mean(logloss_baseline, na.rm = TRUE),
    IMV = (avg_logloss_baseline - avg_logloss_full) / avg_logloss_baseline
  )


ggplot(summary_df, aes(x = theta_bin)) +
  geom_line(aes(y = avg_sem, group = 1, color = "SEM")) +
  geom_line(aes(y = IMV, group = 1, color = "IMV")) +
  labs(
    x = "Theta Bin",
    y = "Value",
    title = "SEM and IMV Across Theta",
    color = "Metric"
  ) +
  theme_minimal()


## check b distribution

item_params_df <- as.data.frame(item_params) %>%
  rownames_to_column(var = "item") %>%
  select(item, a, b) %>% 
  filter(abs(b) < 6)

ggplot(item_params_df, aes(x = b)) +
  geom_histogram(binwidth = 0.5, fill = "#619CFF", color = "black", boundary = 0) +
  labs(
    title = "Distribution of Item Difficulties (b)",
    x = "Item Difficulty (b)",
    y = "Number of Items"
  ) +
  theme_minimal()

ggplot(item_params_df, aes(x = b, y = a)) +
  geom_point(alpha = 0.6, color = "#0072B2") +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 1) +
  labs(
    title = "Relationship Between Item Difficulty (b) and Discrimination (a)",
    x = "Item Difficulty (b)",
    y = "Item Discrimination (a)"
  ) +
  theme_minimal()

long_df %>%
  distinct(id, theta_bin) %>%  # one row per student per bin
  count(theta_bin) %>%         # count students per bin
  ggplot(aes(x = theta_bin, y = n)) +
  geom_col(fill = "#56B4E9", color = "black") +
  labs(
    title = "Number of Students per Theta Bin",
    x = "Theta Bin",
    y = "Number of Students"
  ) +
  theme_minimal()

