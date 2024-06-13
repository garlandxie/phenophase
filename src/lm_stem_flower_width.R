# libraries --------------------------------------------------------------------
library(dplyr)     # for manipulating data
library(ggplot2)   # for visualizing data 
library(here)      # for creating relative file-paths 
library(emmeans)   # for testing pairwise comparisons
library(ggsignif)  # for creating significance bars in plots

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

# check diagnostic plots for the following assumptions
# (1) heteroscedasticity
# (2) normality
# (3) influential outliers
plot(lm_stem_pheno)

# check global significance of phenophase 
# since there are three levels (P2, P3, and P4)
summary(aov(lm_stem_pheno))

# check sample design before calculating estimated marginal means
ref_grid_sp <- emmeans::ref_grid(lm_stem_pheno)
ref_grid_sp@grid

emm_stem_width <- emmeans(lm_stem_pheno, specs = "Pheno")
pairs_stem_width <- emm_stem_width %>%
  pairs() %>%
  as.data.frame()

## flower width ----------------------------------------------------------------

# model fit
lm_flower_pheno <- lm(Width_flower ~ Pheno, data = perc_viab_tidy)

# check diagnostic plots for the following assumptions
# (1) heteroscedasticity
# (2) normality
# (3) influential outliers
plot(lm_flower_pheno)

# check global significance of phenophase
# since there are three levels (P2, P3, and P4)
summary(aov(lm_flower_pheno))

# check sample design before calculating estimated marginal means
ref_grid_fp <- emmeans::ref_grid(lm_flower_pheno)
ref_grid_fp@grid

# pairwise comparisons
emm_flower_width <- emmeans(lm_flower_pheno, specs = "Pheno")
pairs_flower_width <- emm_flower_width %>%
  pairs() %>%
  as.data.frame()

# data visualization -----------------------------------------------------------

## stem width ------------------------------------------------------------------

emm_stem_width_df <- emm_stem_width %>%
  as.data.frame() %>%
  mutate(
    emmean = round(emmean, digits = 2),
    Pheno = case_when(
      Pheno == "P2" ~ "Flowering", 
      Pheno == "P3" ~ "Flowering Completed", 
      Pheno == "P4" ~ "Seeds matured", 
      TRUE ~ Pheno)
  ) 

pairs_sw_p2_p3 <- pairs_stem_width %>%
  filter(contrast == "P2 - P3") %>%
  pull(p.value) %>%
  round(digits = 2)

pairs_sw_p2_p4 <- pairs_stem_width %>%
  filter(contrast == "P2 - P4") %>%
  pull(p.value) %>%
  round(digits = 2)

pairs_sw_p3_p4 <- pairs_stem_width %>%
  filter(contrast == "P3 - P4") %>%
  pull(p.value) %>%
  round(digits = 2)

(plot_stem_width <- perc_viab_tidy %>%
  dplyr::filter(!is.na(Pheno)) %>%
  mutate(
    Pheno = case_when(
      Pheno == "P2" ~ "Flowering", 
      Pheno == "P3" ~ "Flowering Completed", 
      Pheno == "P4" ~ "Seeds matured", 
      TRUE ~ Pheno)
    ) %>%
  ggplot(aes(x = Pheno, y = Width_stem)) + 
  geom_jitter(width = 0.2, alpha = 0.1) + 
  geom_point(
    aes(x = Pheno, y = emmean), 
    color = "red", 
    size = 2.5, 
    shape = "triangle", 
    data = emm_stem_width_df
    ) +
    
  geom_signif(
    y_position = 0.3, 
    xmin = 1, 
    xmax = 2, 
    annotation = paste("p = ", pairs_sw_p2_p3),
    alpha = 0.5
    ) +   
  
  geom_signif(
    y_position = 0.325, 
    xmin = 2, 
    xmax = 3, 
    annotation = paste("p = ", pairs_sw_p3_p4),
    alpha = 0.5
    ) +  
    
   geom_signif(
    y_position = 0.35, 
    xmin = 1, 
    xmax = 3, 
    annotation = paste("p = ", pairs_sw_p2_p4),
    alpha = 0.5
    ) +    
    
  ylim(0.05, 0.4) + 
  labs(x = "Phenophase", y = "Stem width per flower head (mm)") + 
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

## flower width ----------------------------------------------------------------

emm_flower_width_df <- emm_flower_width %>%
  as.data.frame() %>%
  mutate(
    emmean = round(emmean, digits = 2),
    Pheno = case_when(
      Pheno == "P2" ~ "Flowering", 
      Pheno == "P3" ~ "Flowering Completed", 
      Pheno == "P4" ~ "Seeds matured", 
      TRUE ~ Pheno
    )
  )
  
pairs_fw_p2_p3 <- pairs_flower_width %>%
  filter(contrast == "P2 - P3") %>%
  # use p < 0.001 if the p-value is really small 
  mutate(p.value = case_when(
    p.value < 0.001 ~ "<0.001")
  ) %>%
  pull(p.value) 

pairs_fw_p2_p4 <- pairs_flower_width %>%
  filter(contrast == "P2 - P4") %>%
  # use p < 0.001 if the p-value is really small 
  mutate(p.value = case_when(
    p.value < 0.001 ~ "<0.001")
  ) %>%
  pull(p.value) 

pairs_fw_p3_p4 <- pairs_flower_width %>%
  filter(contrast == "P3 - P4") %>%
  # use p < 0.001 if the p-value is really small 
  mutate(p.value = case_when(
    p.value < 0.001 ~ "<0.001")
  ) %>%
  pull(p.value) 

# check for duplicate values before creating jitter plots 
table(perc_viab_tidy$Width_flower,perc_viab_tidy$Pheno)

(plot_flower_width <- perc_viab_tidy %>%
    dplyr::filter(!is.na(Pheno)) %>%
    mutate(Pheno = case_when(
      Pheno == "P2" ~ "Flowering", 
      Pheno == "P3" ~ "Flowering Completed", 
      Pheno == "P4" ~ "Seeds matured", 
      TRUE ~ Pheno) 
      ) %>%
    ggplot(aes(x = Pheno, y = Width_flower)) + 
    geom_jitter(width = 0.2, alpha = 0.1) + 
    geom_point(
      aes(x = Pheno, y = emmean), 
      color = "red", 
      size = 2.5, 
      shape = "triangle", 
      data = emm_flower_width_df
    ) +
  
    geom_signif(
      y_position = 2.0, 
      xmin = 1, 
      xmax = 2, 
      annotation = paste("p = ", as.character(pairs_fw_p2_p3)),
      alpha = 0.5
    ) + 
    
    geom_signif(
      y_position = 2.2, 
      xmin = 2, 
      xmax = 3, 
      annotation = paste("p = ", as.character(pairs_fw_p3_p4)),
      alpha = 0.5
    ) +
    
    geom_signif(
      y_position = 2.4, 
      xmin = 1, 
      xmax = 3, 
      annotation = paste("p = ", as.character(pairs_fw_p2_p4)),
      alpha = 0.5
    ) +
  
    ylim(0.2, 2.5) + 
    labs(x = "Phenophase", y = "Flower width per flower head (mm)") + 
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

# save to disk -----------------------------------------------------------------

ggsave(
  plot = plot_flower_width, 
  filename = here("output", "plot_fw_pheno.png"),
  device = "png", 
  units = "in",
  width = 6, 
  height = 5
)

ggsave(
  plot = plot_stem_width, 
  filename = here("output", "plot_sw_pheno.png"),
  device = "png", 
  units = "in",
  width = 6, 
  height = 5
)
