#Necessary Packages
library(tidyverse)
library(cocor)
library(ggplot2)
library(car)
library(emmeans)
library(patchwork)
library(lme4)
library(lmerTest)
library(ggsignif)
library(multcomp)
#Set WD and correct contrast mode for ANOVA's
setwd("C:/Users/Marce/OneDrive/Desktop/491/Rwork")
options(contrasts = c("contr.sum","contr.poly"))

#Tile heating: load file
dftiles <- read.csv("Tiletemps.csv") %>%
  na.omit() %>%
  mutate(
    Ht = factor(Ht, levels = c(0,1), labels = c("Low","High")),
    Clr = factor(Clr, levels = c(0,1), labels = c("White","Black"))
  )

#Make LM for Tile temperature
lm_TempF <- lm(TempF ~ Ht * Clr, data = dftiles)
qqnorm(resid(lm_TempF))
qqline(resid(lm_TempF))
plot(fitted(lm_TempF), resid(lm_TempF))
#Levene and shapiro tests are unreliable at this small sample size
Anova(lm_TempF, type = 3) #Ht: p = 3.582e-05, Clr: p = 0.0006441, Ht:Clr: p = 0.0397927
#Pairwise Holm testing for posthoc (see notes for why this instead of stratifying and doing one-way ANOVA's)
emmtempf <- emmeans(lm_TempF, ~ Clr | Ht)
pairs((emmtempf), adjust = "holm") # Low: p = 0.0704, High: p = 0.0005
summary(emmtempf)


extract_anova_lm <- function(model, response_name) {
  tab <- as.data.frame(Anova(model, type = 3))
  tab <- tab[!rownames(tab) %in% c("(Intercept)", "Residuals"), ]
  tab$Term <- rownames(tab)
  tab %>%
    transmute(
      Response  = response_name,
      Term,
      Statistic = sprintf("F(%.0f,%.0f) = %.2f", Df, df.residual(model), `F value`),
      p.value   = `Pr(>F)`
    )
}

tile_stats <- extract_anova_lm(lm_TempF, "Tile surface temperature")

saveRDS(tile_stats, "stats_tiletemp.rds")
print(tile_stats)

#SUMMARY STATS ---> RAW MEANS AND STDEVS
summary_stats <- dftiles %>%
  group_by(Ht, Clr) %>%
  summarize(
    Mean = mean(TempF),
    SD = sd(TempF)
  )
print(summary_stats)

######################################################################

dftiles$CodeLabel <- paste(dftiles$Ht,
                           dftiles$Clr,
                           sep = "-")
dftiles$Ht <- factor(dftiles$Ht, levels = c("High","Low"), labels = c("High-Shore\n(11:48 AM)", "Low-Shore\n(12:30 PM)"))

#Boxplot of temperatures
sig_df <- data.frame(
  Ht = factor(c("High-Shore\n(11:48 AM)", "Low-Shore\n(12:30 PM)"), levels = levels(dftiles$Ht)),
  xmin = c(1, 1), xmax = c(2, 2),
  y_position = c(37, 26.5),
  annotations = c("***", "ns")   # p=.0005 -> ***, p=.070 -> ns
)

Tileboxplot <- ggplot(dftiles, aes(x = Clr, y = TempF, fill = Clr)) +
  geom_boxplot(outlier.shape = NA, colour = "black") + 
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~Ht, scales = "free_x") +
  scale_fill_manual(values = c("White" = "white", "Black" = "grey70")) +
  labs(x = "Tile Colour",
       y = "Surface Temperature (\u00B0C)",
       tag = "A") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(18, 41)) +   # <- ADD THIS: headroom for the bracket at y=37
  geom_signif(
    data = sig_df,
    aes(xmin = xmin, xmax = xmax, y_position = y_position, annotations = annotations),
    manual = TRUE,
    tip_length = 0.02,
    textsize = 3,
    inherit.aes = FALSE
  )

Tileboxplot
saveRDS(Tileboxplot, "manufig_tiletemp.rds")









#################################################################################





dftilesraw <- read.csv("Tiletemps2.csv") %>%
  na.omit() %>%
  mutate(
    Ht = factor(Ht, levels = c(0,1), labels = c("Low","High")),
    Clr = factor(Clr, levels = c(0,1), labels = c("White","Black"))
  )

dftilesraw$Ht <- factor(dftilesraw$Ht, levels = c("High","Low"), labels = c("High-Shore\n(11:48 AM)", "Low-Shore\n(12:30 PM)"))

Tileprecplot <- ggplot(dftilesraw, aes(x = factor(Tile), y = TempF)) + 
  #The 3 Raw Measurements per Tile:
  geom_point(alpha = 0.3, size = 3, color = "blue") +
  # A tiny line showing the spread (Range or SD) of those 3 points:
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.2, color = "gray40") +
  facet_grid(. ~ Ht + Clr, scales = "free_x", space = "free") +
  theme_bw() +
  labs(title = "Supplementary Figure: Measurement Precision",
       subtitle = "Raw replicate measurements (n = 3) for each tile (N = 12)",
       y = "Temperature (°C)", 
       x = "Individual Tile ID") +
  theme(strip.text = element_text(size = 10))
Tileprecplot



  