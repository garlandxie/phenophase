# libraries -----
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(here)

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

# Data visualization: percent seed viability ----

## Environmental condition ----

# sample size for percent seed viability between inside + outside environment
perc_viab_raw %>%
  group_by(Treatment) %>%
  summarize(sample_size = sum(!is.na(Perc_viability))) %>%

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

## Phenophase  -----

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

## Day -----

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

# reviewer concern ----

# fit a binomial regression that models proportional data 
# (i.e., the fraction of viable seeds per 20 seeds from a given flower head)
glm_viab <- glm(
  Perc_viability ~ Pheno*Treatment*Mow, 
  family = binomial(link="logit"),
  data = perc_viab_tidy
)

# check diagnostics for assumptions of binomial regression
# independence of errors
# linearity in the logit for continuous variables
# lack of strongly influential outliers
plot(glm_viab)

# use type II ANOVA adjusted sum of squares for the unbalanced design
# would like to see global significance of main effects and interactions
# prior to running pairwise comparisons 
anova_viab <- car::Anova(glm_viab, type = "II")

# check reference grid for unbalanced sample sizes 
ref_grid_viab <- ref_grid(glm_viab)
ref_grid_viab@grid

# visualize the nature of interactions before doing any statistical comparisons
emmeans::emmip(glm_viab, Mow ~ Pheno | Treatment)
emmeans::emmip(glm_viab, Mow ~ Pheno | Treatment, CIs = TRUE)

# calculate estimate marginal means of each environmental condition separately
# by conditioning on phenophase and mowing frequency
viab_emm_pheno_t <- emmeans(glm_viab, pairwise ~ Pheno | Treatment)
viab_emm_mow_t <- emmeans(glm_viab, pairwise ~ Mow | Treatment)

# calculate estimate marginal means for the interaction between
# mowing and phenophase
viab_emm_int <- emmeans(glm_viab, ~ Mow*Pheno)
pairs_viab <- contrast(viab_emm_int, interaction = c("consec", "poly"))

# visualize the effects of mowing and phenophase
# to confirm results provided by the interaction contrasts  
perc_viab_tidy %>%
  filter(!is.na(Mow)) %>%
  ggplot(aes(x = Pheno, y = Perc_viability)) + 
  geom_boxplot() + 
  geom_point(alpha = 0.1) + 
  facet_wrap(~Mow) + 
  ylim(0,1) +
  labs(
    x = "Phenophase",
    y = "Seed Viability (%) per flower head") + 
  theme_bw()
