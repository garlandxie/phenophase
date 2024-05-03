citation()


library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
# Load the data
data2 <- read.csv("~/Desktop/CD.csv")
head(data)
sidra2 <- read.csv("~/Desktop/Sidra2.csv")
head(sidra2)

sidra2 <- read.csv("~/Desktop/Fallaha ms/Sidra_final.csv")
sidra3 <- read.csv("~/Desktop/Fallaha ms/Sidra_predator2.csv")

merged_dataset <- merge(sidra2, sidra3, by = "Individual", all = TRUE)
write <- write.csv(merged_dataset, file = "~/Desktop/Sidra_combined")
write <- write.csv(sidra3, file = "~/Desktop/Sidra_combined.csv")
sidra3 <- read.csv("~/Desktop/Fallaha ms/Sidra_combined.csv")

write <- write.csv(df_filtered, file = "~/Desktop/Sidra_ms_supp.csv")


install.packages("dplyr")
library(dplyr)

# Assuming your original data frame is 'df'
# Replace 'x_column' with the name of your x-axis column
# Replace 'y_column' with the name of your y-axis column
df_filtered <- sidra3 %>% filter(!is.na(Pheno.x))

#plot flower head width
width <- ggplot(df_filtered, aes(x=Pheno.x, y=width_flower)) +
  geom_boxplot() +
  labs(x="Phenophase", y="Flower head width (cm)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(panel.background = element_rect(fill = "white"))
width + theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14), 
              axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)) 

#flower width significantly different between T2, T3, T4
flower <- aov(width_flower~ Pheno.x, data=df_filtered) 
summary(flower)
TukeyHSD(flower)

#plot stem head width
stems <- ggplot(df_filtered, aes(x=Pheno.x, y=width_stem)) +
  geom_boxplot() +
  labs(x="Phenophase", y="Stem width (cm)") +
theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
 theme(panel.background = element_rect(fill = "white"))
stems + theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14), 
              axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)) 

#stem width does not vary too much
stem <- aov(width_stem~ Pheno.x, data=df_filtered)
summary(stem)
TukeyHSD(stem)

#inside or outside effects; inside had higher viability
tt_environment <- t.test(Measurement ~ Treatment, data = df_filtered)
tt_environment
environment <- ggplot(df_filtered, aes(x=Treatment, y=Measurement, fill=Treatment)) +
  geom_boxplot(color="black", fill="white") +
  theme_bw() +
  labs(x="Environmental effects post cutting", y="Seed Viability (%)")
theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(panel.background = element_rect(fill = "white"))
environment + theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14), 
              axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)) 

df_filtered <- sidra3 %>% filter(!is.na(Predator))
predator <- ggplot(df_filtered, aes(x=Predator, y=Measurement)) +
  geom_boxplot() +
  labs(x="Predator", y="Seed Viability (%)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(panel.background = element_rect(fill = "white"))
predator + theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14), 
              axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)) 

decomposed <- ggplot(df_filtered, aes(x=Predator)) +
  geom_bar() +
  facet_wrap(~Decomposed) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(panel.background = element_rect(fill = "white")) +
  labs(x="Predator", y="Frequency") +
decomposed + theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14), 
                 axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)) 

#
Full_model <- aov(Measurement~ Pheno.x + Group_d + Treatment + width_flower + Predator, data=df_filtered)
summary(Full_model)
TukeyHSD(Full_model)

tt_predator <- t.test(Decomposed ~ Predator, data = sidra3)





#effect of environment across all phenophases and start dates
tt_environment <- t.test(Measurement ~ Treatment, data = data2)
tt_environment
ggplot(data, aes(x=Treatment, y=Measurement, fill=Treatment)) +
  geom_boxplot() +
  labs(x="Environmental effects post cutting", y="Seed Viability (%)")

boxplot <- ggplot(data, aes(x = Treatment, y = Measurement)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white")) +
  xlab("Environmental effects post cutting") + 
  ylab("Seed viability") +
  title("Decomposed")
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
boxplot + theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14), 
                axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)) 


#effect of phenophase
phenophase_p <- aov(Measurement~ Phenophase, data=data)
summary(phenophase_p)
TukeyHSD(phenophase_p)
data$Phenophase <- as.factor(data$Phenophase)
boxplot <- ggplot(data, aes(x = Phenophase, y = Measurement)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white")) +
  xlab("Phenophase") + 
  ylab("Seed viability") + 
  scale_x_discrete(labels = c("Flowering", "Flowering complete", "Seeds matured")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
boxplot + theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14), 
                axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)) 


#effect of group (Day)
data$Group_n <- as.factor(data$Group_n)
Group_day <- aov(Measurement~ Group_n, data=data)
summary(Group_day)
TukeyHSD(Group_day)
boxplot <- ggplot(data, aes(x = Group_n, y = Measurement)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white")) +
  xlab("Testing day (days since cut)") + 
  ylab("Seed viability") +
  scale_x_discrete(labels = c("Day 1", "Day 7", "Day 14", "Day 21")) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
boxplot + theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14), 
                axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14)) 


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