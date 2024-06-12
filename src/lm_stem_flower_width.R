# libraries --------------------------------------------------------------------
library(dplyr)     # for manipulating data
library(ggplot2)   # for visualizing data 
library(here)      # for creating relative file-paths 
library(emmeans)   # for testing pairwise comparisons

# import -----------------------------------------------------------------------
perc_viab_raw <- read.csv(here("data", "raw_data.csv"))

# data clean -------------------------------------------------------------------

perc_viab_tidy <- perc_viab_raw %>%
  rename(Mow = Group_d) %>%
  mutate(
    Mow = factor(Mow, levels = c("Day 1", "Day 7", "Day 14", "Day 21")),
    Treatment = factor(Treatment, levels = c("Inside", "Outside")), 
    Pheno = factor(Pheno, levels = c("P2", "P3", "P4"))
  )

# exploratory data analysis ----------------------------------------------------

(plot_flower_width <- perc_viab_tidy %>%
  dplyr::filter(!is.na(Pheno)) %>%
  ggplot(aes(x = Pheno, y = Width_flower)) + 
  geom_boxplot() 
)

(plot_stem_width <- perc_viab_tidy %>%
  dplyr::filter(!is.na(Pheno)) %>%
  ggplot(aes(x = Pheno, y = Width_stem)) + 
  geom_boxplot() 
)

# statistical analysis ---------------------------------------------------------

## stem width ------------------------------------------------------------------

# model fit
lm_stem_pheno <- lm(Width_stem ~ Pheno, data = perc_viab_tidy)
plot(lm_stem_pheno)

car::Anova(lm_stem_pheno, type = "II")

# pairwise comparisons
emm_stem_width <- emmeans(lm_stem_pheno, specs = "Pheno")
pairs_stem_width <- pairs(emm_stem_width)

## flower width ----------------------------------------------------------------

# model fit
lm_flower_pheno <- lm(Width_stem ~ Pheno, data = perc_viab_tidy)
plot(lm_stem_pheno)

car::Anova(lm_flower_pheno, type = "II")

# pairwise comparisons
emm_flower_width <- emmeans(lm_flower_pheno, specs = "Pheno")
pairs_flower_width <- pairs(emm_flower_width)

# data visualization -----------------------------------------------------------

## stem width ------------------------------------------------------------------

emm_stem_width_df <- emm_stem_width %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  mutate(
    emmean = round(emmean, digits = 2),
    lower_cl = round(lower_cl, digits = 2),
    upper_cl = round(upper_cl, digits = 2)
  ) 

(plot_stem_width <- perc_viab_tidy %>%
  dplyr::filter(!is.na(Pheno)) %>%
  ggplot(aes(x = Pheno, y = Width_stem)) + 
  geom_point(alpha = 0.1) + 
  geom_point(aes(x = pheno, y = emmean), data = emm_stem_width_df) +
  ylim(0, 2) + 
  labs(x = "Phenophase", y = "Stem width (mm)") + 
  theme_bw()
)

## flower width ----------------------------------------------------------------

emm_flower_width_df <- emm_flower_width %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  mutate(
    emmean = round(emmean, digits = 2),
    lower_cl = round(lower_cl, digits = 2),
    upper_cl = round(upper_cl, digits = 2)
  ) 

(plot_flower_width <- emm_flower_width_df %>%
    ggplot(aes(x = pheno, y = emmean)) + 
    geom_point() + 
    geom_pointrange(aes(ymin = lower_cl, ymax = upper_cl)) + 
    ylim(0, 0.5) + 
    labs(x = "Phenophase", y = "Flower width (mm)") + 
    theme_bw()
)

## multi-panel figure ----------------------------------------------------------

