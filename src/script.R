# libraries -----
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(here)

# import ----

perc_viab_raw <- read.csv(here("data", "raw_data.csv"))

# data clean ----

perc_viab_tidy <- rename(perc_viab_raw, Mow = Group_d)

# data visualization ----

## Flower head width ----

(plot_flower_width <- perc_viab_raw %>%
   
  # remove missing values from Pheno since "NA" can show up as a x-axis label
  filter(!is.na(Pheno)) %>%
   
  # show median differences of flower head width across phenophase 2,3,and 4
  # phenophase 1 was excluded since flower phenology did not occur here
  ggplot(aes(x = Pheno, y = Width_flower)) + 
  geom_boxplot() +
  labs(x ="Phenophase", y ="Flower head width (cm)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14)
    )
)
  
## Stem head width -----

(plot_stem_width <- perc_viab_raw %>%
   
  # remove missing values from Pheno since "NA" can show up as a x-axis label
  filter(!is.na(Pheno)) %>%
   
  # show median differences of stem width across phenophase 2,3,and 4
  # phenophase 1 was excluded since flower phenology did not occur here
  ggplot(aes(x = Pheno, y = Width_stem)) +
  geom_boxplot() +
  labs(x ="Phenophase", y ="Stem width (cm)") + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )
)


## Percent seed viability ----

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
  
## Predator -----

# visualize the effects of predator on the percent seed viability
(plot_predator <- perc_viab_raw %>%
  filter(!is.na(Predator)) %>% 
  ggplot(aes(x = Predator, y = Perc_viability)) +
  geom_boxplot() +
  labs(x="Predator", y="Seed Viability (%)") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14)
    ) 
)

## Decomposed ----

# visualize the effects of decomposition on percent seed viability
(decomposed <- perc_viab_raw %>%
  filter(!is.na(Predator)) %>%
  ggplot(aes(x=Predator)) +
  geom_bar() +
  facet_wrap(~Decomposed) +
  labs(x="Predator", y="Frequency") + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )
)

## Phenophase on seed viability -----

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

## Effect of day on seed viability -----

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


# statistical analysis -----

## effect of phenophase on flower width ----

# flower width significantly different between T2, T3, T4
aov_flower <- aov(Width_flower ~ Pheno, data = perc_viab_raw) 
summary(aov_flower)
TukeyHSD(aov_flower)

## effect of phenophase on stem width ----

# stem width does not vary too much
aov_stem <- aov(Width_stem ~ Pheno, data = perc_viab_raw)
summary(aov_stem)
TukeyHSD(aov_stem)

## effect of environmental conditions on viability ----

# inside or outside effects; inside had higher viability
lm_environment <- glm(Perc_viability ~ Treatment, data = perc_viab_ra)
tt_environment


# reviewer concern ----
# proportional data -> binomal glm
# three-way interaction
# marginal means to account for unbalanced sample size
# sum of squares II or III? 

Full_model <- aov(Measurement~ Pheno.x + Group_d + Treatment + width_flower + Predator, data=df_filtered)
summary(Full_model)
TukeyHSD(Full_model)

tt_predator <- t.test(Decomposed ~ Predator, data = sidra3)

#effect of environment across all phenophases and start dates
tt_environment <- t.test(Measurement ~ Treatment, data = data2)
tt_environment



#effect of phenophase
phenophase_p <- aov(Measurement~ Phenophase, data=data)
summary(phenophase_p)
TukeyHSD(phenophase_p)
data$Phenophase <- as.factor(data$Phenophase)



#effect of group (Day)
data$Group_n <- as.factor(data$Group_n)
Group_day <- aov(Measurement~ Group_n, data=data)
summary(Group_day)
TukeyHSD(Group_day)



#effect of phenophase, testing date, and treatment
data_final$Phenophase <- as.factor(data_final$Phenophase)
Full_model <- aov(Measurement~ Pheno.x + Group_d + Treatment + width_flower, data=data_final)
summary(Full_model)
TukeyHSD(Full_model)

residuals <- residuals(Full_model)
hist(residuals, main = "Histogram of Residuals", xlab = "Residuals")
shapiro.test(residuals)
qqnorm(residuals)
qqline(residuals)

scatterplot <- lm(Measurement~width_flower, data=df_filtered)
summary(scatterplot)
scatty <- ggplot(df_filtered, aes(y = Measurement, x = width_flower)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  geom_jitter() +
  labs(y = "Seed viability", x = "Flower width (cm)")
scatty + theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14), 
                axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)) 


data_p4 <- data_final[data_final$Pheno.x != "P2",]
data_p2 <- data_final[data_final$Pheno.x != "P3",]
data_p2real <- data_p2[data_p2$Pheno.x != "P4",]
data_p4real <- data_p4[data_p4$Pheno.x != "P3",]
data_p3real <- data_p4[data_p4$Pheno.x != "P4",]

Full_model_p2 <- aov(Measurement~ Group_d + Treatment + width_flower, data=data_p2real)
summary(Full_model_p2)
custom_colors <- c("Inside" = "white", "Outside" = "lightgrey")
boxplot <- ggplot(data_p2real, aes(x = Group, y = Measurement, fill= Treatment)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white")) +
  xlab("Testing day (days since cut)") + 
  ylab("Seed viability") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_x_discrete(labels = c("Day 1", "Day 7", "Day 14", "Day 21")) +
  scale_fill_manual(values = custom_colors)
boxplot
boxplot + theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14), 
                axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)) 

Full_model_p3 <- aov(Measurement~ Group_d + Treatment + width_flower, data=data_p3real)
summary(Full_model_p3)
boxplot <- ggplot(data_p3real, aes(x = Group_d, y = Measurement, fill= Treatment)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white")) +
  xlab("Testing day (days since cut)") + 
  ylab("Seed viability") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_fill_manual(values = custom_colors)
boxplot + theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14), 
                axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)) 

Full_model_p4 <- aov(Measurement~ Group + Treatment + width_flower, data=data_p4real)
summary(Full_model_p4)
boxplot <- ggplot(data_p4real, aes(x = Group_d, y = Measurement, fill= Treatment)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white")) +
  xlab("Testing day (days since cut)") + 
  ylab("Seed viability") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_fill_manual(values = custom_colors)
boxplot + theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14), 
                axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)) 


# Filter for Phenophase 3
data_phenophase_3 <- filter(data, Phenophase == Third)
# Convert ‘Group’ to a factor with levels in the desired order
data_phenophase_3$Group <- factor(data_phenophase_3$Group,
                                  levels = c(“July_25_2023", “August_1_2023”, “August_7_2023", “August_14_2023”))
# Rename levels of ‘Group’ to “Day 1", “Day 8”, “Day 14", “Day 21”
data_phenophase_3$Group <- factor(data_phenophase_3$Group,
                                  levels = c(“July_25_2023", “August_1_2023”, “August_7_2023", “August_14_2023”),
                                  labels = c(“Day 1", “Day 8”, “Day 14", “Day 21”))
# Assuming ‘Treatment’ is also a factor.
data_phenophase_3$Treatment <- as.factor(data_phenophase_3$Treatment)
# Create the box plot
ggplot(data_phenophase_3, aes(x=Group, y=Measurement, fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title=“Phenophase 3", x=“Days Post Pruning”, y=“Seed Viability(%)“)
# Perform ANOVA
Aov.test <- aov(Measurement ~ Group * Treatment, data = data_phenophase_3)
# To view the summary of the ANOVA
summary(Aov.test)