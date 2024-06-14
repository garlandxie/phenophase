# libraries -----
library(car)       # for testing ANOVA type II models 
library(dplyr)     # for manipulating data  
library(ggplot2)   # for visualizing data
library(here)      # for creating relative file-paths
library(emmeans)   # for calculating least-square means 

# import ----

perc_viab_raw <- read.csv(here("data", "raw_data.csv"))

# data clean ----

perc_viab_tidy <- perc_viab_raw %>%
  rename(Mow = Group_d) %>%
  mutate(
    Mow = factor(Mow, levels = c("Day 1", "Day 7", "Day 14", "Day 21")),
    Treatment = factor(Treatment, levels = c("Inside", "Outside")), 
    Pheno = factor(Pheno, levels = c("P2", "P3", "P4"))
  )

# exploratory data analysis: percent seed viability ----

## environmental condition ----

# sample size for percent seed viability between inside + outside environment
perc_viab_raw %>%
  group_by(Treatment) %>%
  summarize(sample_size = sum(!is.na(Perc_viability))) 

(plot_env <- perc_viab_raw %>%

  # remove missing values from treatment (inside and outside enviromment)
  # since "NAs" can show up in the x-axis label 
  filter(!is.na(Treatment)) %>%
   
  ggplot(aes(x = Treatment, y = Perc_viability, fill = Treatment)) + 
  
  # tried to use boxplots but visually difficult to interpret
  # try violin plots (for now) since there appears to be many values 
  # close to zero for the outside environmental conditions
  geom_violin() + 
  
  # show mean differences between percent seed viability between
  # inside and outside environmental conditions 
  # cross bar is a bit extravagant, but it works for now
   stat_summary(fun = "mean",
                geom = "crossbar",
                color = "black") + 
 
  labs(
    x = "Environmental Conditions",
    y = "Percent seed viability (%)"
    ) + 
  
  # range of percent seed viability should be 0 to 1 
  ylim(0,1) + 
   
  theme_bw()
)

## phenophase  -----

(plot_phen_viab <- perc_viab_raw %>%
  filter(!is.na(Pheno)) %>%
  ggplot(aes(x = Pheno, y = Perc_viability)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  xlab("Phenophase") + 
  ylab("Seed viability") + 
  scale_x_discrete(labels = c("Flowering", "Flowering complete", "Seeds matured")) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x = element_text(size = 14), 
    axis.title.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14), 
    axis.title.y = element_text(size = 14)
    ) 
)

## day -----

(plot_day_viab <- perc_viab_raw %>%
  filter(!is.na(Group_n)) %>%
  mutate(Group_n = factor(Group_n)) %>%
  ggplot(aes(x = Group_n, y = Perc_viability)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white")) +
  xlab("Testing day (days since cut)") + 
  ylab("Seed viability") +
  scale_x_discrete(labels = c("Day 1", "Day 7", "Day 14", "Day 21")) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x = element_text(size = 14), 
    axis.title.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14), 
    axis.title.y = element_text(size = 14)
    ) 
)

# statistical analyses: percent seed viability ----

## model fit ----

# fit a binomial regression that models proportional data 
# (i.e., the fraction of viable seeds per 20 seeds from a given flower head)
glm_viab <- glm(
  Perc_viability ~ Pheno*Treatment*Mow, 
  family = binomial(link="logit"),
  data = perc_viab_tidy
)

# calculate dispersion parameter to determine evidence of 
# overdispersion or underdispersion
sum(residuals(glm_viab, "pearson")^2) / glm_viab$df.residual

# refit binomial regression to account for underdispersion
# using the quasibinomial family 
glm_viab_und <- glm(
  Perc_viability ~ Pheno*Treatment*Mow, 
  family = quasibinomial(link="logit"),
  data = perc_viab_tidy
)

# check diagnostics for assumptions of binomial regression
# independence of errors
# linearity in the logit for continuous variables
# lack of strongly influential outliers
plot(glm_viab)
plot(glm_viab_und)

## anova ----

# use type III ANOVA adjusted sum of squares for the unbalanced design
# would like to see global significance of main effects and interactions
# prior to running pairwise comparisons 
anova_viab <- car::Anova(glm_viab_und, type = "II", test.statistic = "F")

# check reference grid for unbalanced sample sizes 
ref_grid_viab <- emmeans::ref_grid(glm_viab_und)
ref_grid_viab@grid

## contrast coding ----

# visualize the nature of interactions before doing any statistical comparisons
emmeans::emmip(glm_viab_und, Mow ~ Pheno | Treatment)
emmeans::emmip(glm_viab_und, Treatment ~ Pheno | Mow)

# calculate estimate marginal means of each environmental condition separately
# by conditioning on phenophase and mowing frequency
viab_emm_pheno_t <- emmeans(glm_viab, pairwise ~ Pheno | Treatment)
viab_emm_mow_t <- emmeans(glm_viab, pairwise ~ Mow | Treatment)

# calculate estimate marginal means for the interaction between

# mowing and phenophase
emm_mow_pheno <- emmeans(glm_viab_und, ~ Mow*Pheno)

pheno_over_time <- emmeans::contrast(
  emm_mow_pheno, 
  interaction = c("poly", "del.eff")
  )

pheno_per_week <- emmeans::contrast(
  emm_mow_pheno, 
  interaction = c("consec", "poly")
)

# pheno and treatment
emm_pheno_trt <- emmeans(glm_viab_und, ~ Pheno*Treatment)
pairs_pheno_trt <- pairs(emm_pheno_trt) 
pairs_pheno_trt_tidy <- pairs_pheno_trt %>%
  as.data.frame() %>%
  dplyr::filter(
    contrast %in% c(
      "P2 Inside - P2 Outside",
      "P3 Inside - P3 Outside",
      "P4 Inside - P4 Outside"
    )
  ) %>%
  mutate(
    back_estimate = boot::inv.logit(estimate), 
    back_estimate = round(back_estimate, digits = 2), 
    
    back_se = boot::inv.logit(SE), 
    back_se = round(back_se, digits = 2)
  ) 

# visualize the effects of mowing and phenophase
# to confirm results provided by the interaction contrasts  
(perc_viab_tidy %>%
  filter(!is.na(Mow)) %>%
  ggplot(aes(x = Mow, y = Perc_viability)) + 
  geom_boxplot() + 
  geom_point(alpha = 0.1) + 
  facet_wrap(~Pheno) + 
  ylim(0,1) +
  labs(
    x = "Phenophase",
    y = "Seed Viability (%) per flower head") + 
  theme_bw()
)

(perc_viab_tidy %>%
    filter(!is.na(Mow)) %>%
    ggplot(aes(x = Pheno, y = Perc_viability)) + 
    geom_point(alpha = 0.1) + 
    stat_summary(fun = "mean", colour = "red", geom = "line") + 
    facet_wrap(~Mow) + 
    ylim(0,1) +
    labs(
      x = "Phenophase",
      y = "Seed Viability (%) per flower head") + 
    theme_bw()
)

# data visualization ----

## interaction contrast - phenophase per week ----

emm_mow_pheno_df <- emm_mow_pheno %>%
  as.data.frame() %>%
  mutate(
    back_emmean = boot::inv.logit(emmean), 
    back_emmean = round(back_emmean, digits = 2),
    
    back_lcl = boot::inv.logit(asymp.LCL), 
    back_ucl = boot::inv.logit(asymp.UCL)
  ) %>% 
  dplyr::filter(Mow %in% c("Day 7", "Day 14", "Day 21")) %>%
  mutate(
    Mow = case_when(
      Mow == "Day 7"  ~ "Week 1 (Days 1 -7)", 
      Mow == "Day 14" ~ "Week 2 (Days 7 - 14)", 
      Mow == "Day 21" ~ "Week 3 (Days 14 - 21)", 
      TRUE ~ Mow),
    
    Pheno = case_when(
      Pheno == "P2" ~ "Flower open (P2)", 
      Pheno == "P3" ~ "Flower maturation (P3)", 
      Pheno == "P4" ~ "Dehiscence (P4)", 
      TRUE ~ Pheno), 
    
    Pheno = factor(
      Pheno, levels = c(
        "Flower open (P2)", 
        "Flower maturation (P3)",
        "Dehiscence (P4)")
      )
  )

(plot_pheno_per_week <- perc_viab_tidy %>%
  mutate(Mow = as.character(Mow)) %>%
  filter(!is.na(Mow) & Mow %in% c("Day 7", "Day 14", "Day 21")) %>%
  mutate(
    Mow = case_when(
      Mow == "Day 7"  ~ "Week 1 (Days 1 -7)", 
      Mow == "Day 14" ~ "Week 2 (Days 7 - 14)", 
      Mow == "Day 21" ~ "Week 3 (Days 14 - 21)", 
     TRUE ~ Mow),
    Mow = factor(Mow),
    Pheno = case_when(
      Pheno == "P2" ~ "Flower open (P2)", 
      Pheno == "P3" ~ "Flower maturation (P3)", 
      Pheno == "P4" ~ "Dehiscence (P4)", 
      TRUE ~ Pheno),
    Pheno = factor(
      Pheno, levels = c(
        "Flower open (P2)", 
        "Flower maturation (P3)",
        "Dehiscence (P4)")
      )
    )
   %>%
  ggplot(aes(x = Pheno, y = Perc_viability)) + 
  geom_jitter(width = 0.1, alpha = 0.2) + 
  geom_point(
    aes(x = Pheno, y = back_emmean), 
    size = 2.5, 
    shape = "triangle", 
    color = "red", 
    data = emm_mow_pheno_df
      ) +  
  facet_wrap(~Mow, nrow = 3, ncol = 1) + 
  labs(
    x = "Phenophase",
    y = "Viability (%) of 20 seeds per flower head") + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(
    margin = margin(t = 10, r = 0, b = 0, l = 0),
    size = 12),
  
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(
    margin = margin(t = 0, r = 10, b = 0, l = 0),
    size = 12)
    )
)

## interaction contrast - phenophase and storage ----
emm_pheno_trt_df <- emm_pheno_trt %>%
  as.data.frame() %>%
  mutate(
    back_emmean = boot::inv.logit(emmean), 
    back_emmean = round(back_emmean, digits = 2),
    
    back_lcl = boot::inv.logit(asymp.LCL), 
    back_ucl = boot::inv.logit(asymp.UCL)
  ) %>%
  mutate(
    Pheno = case_when(
      Pheno == "P2" ~ "Flowering", 
      Pheno == "P3" ~ "Flowering completed", 
      Pheno == "P4" ~ "Seeds matured", 
      TRUE ~ Pheno)
  )

(plot_pheno_trt_int <- perc_viab_tidy %>%
  mutate(
    Pheno = case_when(
      Pheno == "P2" ~ "Flowering", 
      Pheno == "P3" ~ "Flowering completed", 
      Pheno == "P4" ~ "Seeds matured", 
      TRUE ~ Pheno)
    ) %>% 
  dplyr::filter(!is.na(Pheno)) %>%
  ggplot(aes(x = Treatment, y = Perc_viability)) +
  geom_jitter(width = 0.1, alpha = 0.1) + 
  geom_point(
    aes(x = Treatment, y = back_emmean), 
      size = 2.5, 
      shape = "triangle", 
      color = "red", 
      data = emm_pheno_trt_df
    ) + 
  facet_wrap(~Pheno) + 
  labs(
    x = "Storage Environment", 
    y = "Viability (%) of 20 seeds per flower head"
    ) + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(
      margin = margin(t = 10, r = 0, b = 0, l = 0),
      size = 12),
    
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(
      margin = margin(t = 0, r = 10, b = 0, l = 0),
      size = 12)
  )
)

# save to disk -----

anova_viab %>%
  janitor::clean_names() %>%
  tibble::rownames_to_column(var = "variables") %>%
  mutate(pr_chisq = round(pr_f, digits = 2)) %>%
  write.csv(
  file = here("output", "anova_viab.csv"), 
  row.names = FALSE
)

pheno_per_week %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  mutate(
    estimate = round(estimate, digits = 2),
    se = round(se, digits = 2),
    z_ratio = round(z_ratio, digits = 2)
    ) %>%
  write.csv(
    file = here("output", "pheno_per_week.csv"), 
    row.names = FALSE
  )

pheno_over_time %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  mutate(
    estimate = round(estimate, digits = 2),
    se = round(se, digits = 2),
    z_ratio = round(z_ratio, digits = 2), 
    p_value = round(p_value, digits = 2)
  ) %>%
  write.csv(
    file = here("output", "pheno_over_time.csv"), 
    row.names = FALSE
  )

ggsave(
  filename = here("output", "plot_pheno_per_week.png"), 
  plot = plot_pheno_per_week, 
  device = "png", 
  units = "in", 
  width = 5, 
  height = 5
)

ggsave(
  filename = here("output", "plot_pheno_trt_int.png"), 
  plot = plot_pheno_trt_int, 
  device = "png", 
  units = "in", 
  width = 7, 
  height = 4
)