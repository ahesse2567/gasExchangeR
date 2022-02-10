library(tidyverse)

test <- read_csv("data/processed/fixed_speeds_grades/mar22_105_pre.csv")

ggplot(data = test[test[["stage"]] >= 6,], aes(x = ex_time, y = vo2)) +
  geom_line() +
  # geom_point(alpha = 0.5) +
  theme_bw() +
  xlab("Time (min)") +
  ylab("VO2") +
  ggtitle("Unaveraged Gas Exchange Data Can be Highly Variable")
