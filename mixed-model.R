#### Mixed-Effects Model on data from Study One of Ackerlund Brandt et al. (2015) ####

#### Authors: William Brady DeHart & Brent A. Kaplan ####

#### Read in library packages ####

library(tidyverse) #Nice package for data manipulation/preparation. Includes ggplot for graphs
library(broom) #Data manipulation
library(lme4) #Mixed-model package. Several packages exist but this allows for specifically selecting auto-regressive correlations among repeated data.
library(emmeans)
library(ggpubr)
library(MuMIn)
library(lmerTest)
library(psych)

#### Functions ####
## https://stackoverflow.com/questions/16249570/uppercase-the-first-letter-in-data-frame
capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}

#### Read in data ####
Final <- read.csv("data/extracted-data.csv", header = T, sep = ",")
Final <- Final %>%
  select(child, condition, session, frequency)

# Re-order session so that control is the comparison group
Final$condition <- factor(Final$condition, 
                          levels = c("control", "experimenter", "child"))
Final$child <- capFirst(Final$child)
Final$condition <- capFirst(Final$condition)
Final$Condition <- Final$condition

#### Plot of data by subject ####
png("Single Subject.png", width = 20, height = 12, units = "in", res = 300)
ggplot(Final, aes(x = session, y = frequency, group = Condition)) + 
  geom_line(size = 1) +
  geom_point(aes(shape = Condition, fill = Condition), size = 5, stroke = 1.5) +
  scale_shape_manual(values = c(22, 21, 21)) + 
  scale_fill_manual(values = c("white", "white", "black")) +
  scale_x_continuous(limits = c(1,14), breaks = c(seq(1,14,2)), 
                     labels = c( "1", "3", "5", "7", "9", "11", "13")) +
  scale_y_continuous(limits = c(-1,17), 
                     breaks = c(0, 5, 10, 15), 
                     labels = c("0", "5", "10", "15")) +
  theme_bw(base_size = 28) +
  ylab("Selections") +
  xlab("Session") +
  ggtitle("Single-Subject Results") +
  geom_text(data=Final,x=14.5,y=16.5,aes(label=child),
            vjust = "inward", hjust = "inward", size = 7)+
  theme(axis.text = element_text(color = "black"), 
        panel.grid = element_blank(), 
        axis.line = element_line(color = "black", size = 1), 
        plot.title = element_text(hjust = 0.5), 
        panel.border = element_rect(size = 1), 
        axis.ticks = element_line(size = 1), 
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  facet_wrap(~child, ncol = 3)
dev.off()

#### Mixed-Model ####

### Model Fit
model.1 <- lmerTest::lmer(frequency ~ condition*session + 
                            (1 + session | child/condition), 
                          data = Final, REML = TRUE) 
# REML = TRUE necessary for Kenward_Roger adjustment.
summary(model.1, ddf = "Kenward-Roger")

### Pairwise Comparisons
emmeans(model.1, pairwise ~ condition)

### R2
r.squaredGLMM(model.1)

#### Prediction Graphs ####

### Extract predicted values
Final <- Final %>% 
  mutate(pred_dist = fitted(model.1)) 

png("Predictions.png", width = 20, height = 12, units = "in", res = 300)
ggplot(Final, aes(x=session, y=pred_dist, group=Condition)) + 
  geom_smooth(aes(group = Condition, linetype = Condition), method='lm', 
              fullrange = TRUE, color = "black", size = 1.5, se=FALSE) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_x_continuous(limits = c(1,14), breaks = c(seq(1,14,2)), 
                     labels = c( "1", "3", "5", "7", "9", "11", "13")) +
  scale_y_continuous(limits = c(-1,17), breaks = c(0, 5, 10, 15), 
                     labels = c("0", "5", "10", "15")) +
  geom_text(data=Final,x=.75,y=16.75,aes(label=child),
            vjust = "inward", hjust = "inward", size = 7)+
  xlab("Session") +
  ylab("Predicted Behavior") +
  ggtitle("Individual Participant Model Predictions") +
  theme_bw(base_size = 24) +
  theme(axis.text = element_text(color = "black"), 
        panel.grid = element_blank(), 
        axis.line = element_line(color = "black", size = 1), 
        plot.title = element_text(hjust = 0.5), 
        panel.border = element_rect(size = 1), 
        axis.ticks = element_line(size = 1), 
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        legend.key.width = unit(3, "line")) +
  facet_wrap(~child, ncol = 3)
dev.off()

#### Diagnostic Plots ####
Plot1 <- ggplot(model.1, aes(x = frequency, y = .fitted)) +
  geom_point(size = 4) +
  stat_smooth(method = "loess", size = 1.5, color = "black", lty = "dashed") +
  scale_x_continuous(limits = c(0,15), breaks = c(0, 3, 6, 9, 12, 15), 
                     labels = c("0", "3", "6", "9", "12", "15")) +
  scale_y_continuous(limits = c(-2,15), breaks = c(0, 5, 10, 15), 
                     labels = c("0", "5", "10", "15")) +
  xlab("Frequency") +
  ylab("Fitted") +
  ggtitle("Test of Linearity") +
  theme_bw(base_size = 24) +
  theme(axis.text = element_text(color = "black"), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(size = 1), 
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(hjust = 0.5))

Plot2 <- ggplot(model.1, aes(.fitted, .resid)) +
  geom_point(size = 4) +
  stat_smooth(method="loess", size = 1.5, color = "black", lty = 2) +
  geom_hline(yintercept=0, col="black", linetype="solid", size = 1) +
  scale_x_continuous(limits = c(-.5,15), breaks = c(0, 3, 6, 9, 12, 15), 
                     labels = c("0", "3", "6", "9", "12", "15")) +
  xlab("Fitted values")+ylab("Residuals") +
  ggtitle("Residual vs Fitted Plot") +
  theme_bw(base_size = 24) +
  theme(axis.text = element_text(color = "black"), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(size = 1), 
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(hjust = 0.5))

Plot3 <- ggplot(model.1, aes(sample = .scresid)) +
  stat_qq(size = 4, shape = 1, stroke = 1.5) +
  stat_qq_line(size = 1.5, linetype = "dashed") +
  xlab("Theoretical Quantiles") +
  ylab("Standardized Residuals") +
  ggtitle("Normal Q-Q") + 
  theme_bw(base_size = 24) +
  theme(axis.text = element_text(color = "black"), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(size = 1), 
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(hjust = 0.5))

Plot4 <- ggplot(model.1,(aes(x = .resid))) + 
  geom_histogram(color = "black", fill = "lightgray") +
  scale_x_continuous(limits = c(-8,8), breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8), 
                     labels = c("-8", "-6", "-4", "-2", "0", "2", "4", "6", "8")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme_bw(base_size = 24) +
  ggtitle("Histogram of Residuals") +
  xlab("Residuals") +
  ylab("Count") +
  theme(axis.text = element_text(color = "black"), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(hjust = 0.5))

png("Diagnostics.png", width = 20, height = 12, units = "in", res = 300)
ggarrange(Plot1, Plot2, Plot3, Plot4)
dev.off()

