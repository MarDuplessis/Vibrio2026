library(tidyverse)
library(cocor)
library(ggplot2)
library(car)
library(emmeans)
library(patchwork)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggsignif)

#Set WD and correct contrast mode for ANOVA's
setwd("C:/Users/Marce/OneDrive/Desktop/491/Rwork")
options(contrasts = c("contr.sum","contr.poly"))

#import, clean, convert 0/1 treatments to factors
dftemps <- read.csv("Oystertemps.csv") %>%
  na.omit() %>%
  mutate(
    Tide_ht = factor(Tide_ht, levels = c(0,1), labels = c("Low","High")),
    Tile_clr = factor(Tile_clr, levels = c(0,1), labels = c("White","Black")),
    ShellClr = factor(ShellClr, levels = c(0,1), labels = c("Light","Dark"))
  )

#NTERNAL TEMPERATURE Primary LMM
lmm_tempin <- lmer(Temp_in ~ Tide_ht * Tile_clr * ShellClr + (1|Tile_no), data = dftemps)
qqnorm(resid(lmm_tempin))
qqline(resid(lmm_tempin))
plot(fitted(lmm_tempin), resid(lmm_tempin))
qqnorm(ranef(lmm_tempin)$Tile_no[[1]], main = "Random Effects QQ")
qqline(ranef(lmm_tempin)$Tile_no[[1]])
isSingular(lmm_tempin)
anova(lmm_tempin, type = 3, ddf = "Satterthwaite")

#EXTERNAL TEMPERATURE Primary LMM
lmm_tempex <- lmer(Temp_ex ~ Tide_ht * Tile_clr * ShellClr + (1|Tile_no), data = dftemps)
qqnorm(resid(lmm_tempex))
qqline(resid(lmm_tempex))
plot(fitted(lmm_tempex), resid(lmm_tempex))
qqnorm(ranef(lmm_tempex)$Tile_no[[1]], main = "Random Effects QQ")
qqline(ranef(lmm_tempex)$Tile_no[[1]])
isSingular(lmm_tempex)
anova(lmm_tempex, type = 3, ddf = "Satterthwaite") 

#EMMeans: INTERNAL
emmtempin_ht <- emmeans(lmm_tempin, ~ Tide_ht, lmer.df = "satterthwaite")
summary(emmtempin_ht)
pairs(emmtempin_ht)
emmtempin_shell <- emmeans(lmm_tempin, ~ ShellClr, lmer.df = "satterthwaite")
summary(emmtempin_shell)
pairs(emmtempin_shell)
#EMMeans: EXTERNAL
emmtempex_ht <- emmeans(lmm_tempex, ~ Tide_ht, lmer.df = "satterthwaite")
summary(emmtempex_ht)
pairs(emmtempex_ht)
emmtempex_shell <- emmeans(lmm_tempex, ~ ShellClr, lmer.df = "satterthwaite")
summary(emmtempex_shell)
pairs(emmtempex_shell)



#REGRESSION: SIMPLE INTERNAL:EXTERNAL 
# Data preparation - add a combined group variable for colouring
dftemps_plot <- dftemps %>%
  mutate(Group = interaction(Tide_ht, ShellClr, sep = " / "))

# Fit simple lm — used ONLY for the plotted regression line/confidence band and R²
lm_temps_plot <- lm(Temp_in ~ Temp_ex, data = dftemps)

# Fit LMM — used for the REPORTED, statistically defensible p-value
lmm_temps_plot <- lmer(Temp_in ~ Temp_ex + (1|Tile_no), data = dftemps)
lmm_temps_stats <- anova(lmm_temps_plot, type = 3, ddf = "Satterthwaite")

r2 <- summary(lm_temps_plot)$r.squared
pval_lmm <- lmm_temps_stats$`Pr(>F)`[1]
pval_label <- ifelse(pval_lmm < .001, "p < 0.001", paste0("p = ", round(pval_lmm, 3)))

# Build the plot
TempRegression_plot <- ggplot(dftemps_plot, aes(x = Temp_ex, y = Temp_in)) +
  geom_smooth(method = "lm", color = "grey40", fill = "grey85", 
              linewidth = 0.8, alpha = 0.4) +
  geom_point(size = 3, alpha = 0.85, color = "grey30") +
  annotate("text", 
           x = min(dftemps$Temp_ex, na.rm = TRUE) + 0.3,
           y = max(dftemps$Temp_in, na.rm = TRUE) - 0.3,
           label = paste0(" R\u00B2 = ", round(r2, 2), "\n", pval_label),
           hjust = 0, vjust = 1, size = 3, color = "grey20") +   # 4 -> 3
  labs(x = "External Shell Temperature (\u00B0C)",
       y = "Internal Tissue Temperature (\u00B0C)",
       caption = NULL, 
       tag = "E") +
  theme_bw(base_size = 11) +   # 20 -> 11
  theme(legend.position = "none",
        plot.caption = element_text(hjust = 0.5, face = "italic", 
                                    color = "grey30", size = 7))   # matches other panels' caption size
#View the plot
TempRegression_plot
saveRDS(TempRegression_plot, "manu_tempregplot.rds")


#to report in results TEXT:
lmm_temps_stats
lmm_temps_plot

#lmm_tempregress <- lmer(Temp_in ~ Temp_ex + Tide_ht + Tile_clr + ShellClr + (1|Tile_no), 
                        #data = dftemps)
#summary(lmm_tempregress)
#anova(lmm_tempregress, type = 3)

#VIOLIN PLOT: HEIGHT EFFECT
emmtempin_ht_df <- as.data.frame(emmtempin_ht)
Effectplot_ht <- ggplot(dftemps, aes(x = Tide_ht, y = Temp_in, fill = Tide_ht)) +
  # 1. The Violin & Jitter Layers
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_jitter(width = 0.05, alpha = 0.3, color = "grey35") +
  # 2. MODEL EMM ± 95% CI
  geom_errorbar(data = emmtempin_ht_df, inherit.aes = FALSE,
                aes(x = Tide_ht, y = emmean, ymin = lower.CL, ymax = upper.CL), 
                width = 0.05, linewidth = 0.8, color = "grey20") +
  geom_point(data = emmtempin_ht_df, inherit.aes = FALSE,
             aes(x = Tide_ht, y = emmean), 
             size = 3, color = "grey20", shape = 16) +
  # 3. Aesthetics and Formatting — print-scale base, no oversized overrides
  scale_fill_manual(values = rep("grey80", length(unique(dftemps$Tide_ht)))) +
  labs(x = "Shore Height", 
       y = "Internal Temperature (\u00B0C)",
       caption = "Pink Circles: Model-Estimated Marginal Geometric Mean \u00B1 95% CI") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        plot.caption = element_text(hjust = 0.5, face = "italic", color = "grey30")) + 
  coord_cartesian(ylim = c(16,31.5))

#VIOLIN Plot: Colour effect
emmtempin_shell_df <- as.data.frame(emmtempin_shell)
Effectplot_clr <- ggplot(dftemps, aes(x = ShellClr, y = Temp_in, fill = ShellClr)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_jitter(width = 0.05, alpha = 0.3, color = "grey35") + 
  geom_errorbar(data = emmtempin_shell_df, inherit.aes = FALSE,
                aes(x = ShellClr, y = emmean, ymin = lower.CL, ymax = upper.CL), 
                width = 0.05, linewidth = 0.8, color = "grey20") +
  geom_point(data = emmtempin_shell_df, inherit.aes = FALSE,
             aes(x = ShellClr, y = emmean), 
             size = 3, color = "grey20", shape = 16) +
  scale_fill_manual(values = rep("grey80", length(unique(dftemps$ShellClr)))) +
  labs(x = "Shell Colour", 
       y = "Internal Temperature (\u00B0C)",
       caption = "Pink Circles: Model EMM \u00B1 95% CI") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        plot.caption = element_text(hjust = 0.5, face = "italic", color = "grey30")) + 
  coord_cartesian(ylim = c(16,31.5))

#GLUE TOGETHER CLR and HT VIOLIN PLOTS
Plot_A <- Effectplot_ht + 
  geom_signif(comparisons = list(c("Low", "High")),
              annotations = "***",   # was "p < 0.001"
              y_position = 30, textsize = 3, vjust = -0.1) +
  labs(caption = NULL, tag = "A")

Plot_B <- Effectplot_clr + 
  geom_signif(comparisons = list(c("Light", "Dark")),
              annotations = "***",   # was "p < 0.001"
              y_position = 30, textsize = 3, vjust = -0.1) +
  labs(caption = NULL, tag = "B")

Plot_A
Plot_B


#EXTERNAL TEMP NOW:
emmtempex_ht_df <- as.data.frame(emmtempex_ht)
Effectplot_ht_ex <- ggplot(dftemps, aes(x = Tide_ht, y = Temp_ex, fill = Tide_ht)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_jitter(width = 0.05, alpha = 0.3, color = "grey35") +
  geom_errorbar(data = emmtempex_ht_df, inherit.aes = FALSE,
                aes(x = Tide_ht, y = emmean, ymin = lower.CL, ymax = upper.CL),
                width = 0.05, linewidth = 0.8, color = "grey20") +
  geom_point(data = emmtempex_ht_df, inherit.aes = FALSE,
             aes(x = Tide_ht, y = emmean),
             size = 3, color = "grey20", shape = 16) +
  scale_fill_manual(values = rep("grey80", length(unique(dftemps$Tide_ht)))) +
  labs(x = "Shore Height", 
       y = "Shell Temperature (\u00B0C)",
       caption = "Pink: Model-Estimated Marginal Mean \u00B1 95% CI") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        plot.caption = element_text(hjust = 0.5, face = "italic", color = "grey30")) +
  coord_cartesian(ylim = c(16,31.5))

emmtempex_shell_df <- as.data.frame(emmtempex_shell)
Effectplot_clr_ex <- ggplot(dftemps, aes(x = ShellClr, y = Temp_ex, fill = ShellClr)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_jitter(width = 0.05, alpha = 0.3, color = "grey35") +
  geom_errorbar(data = emmtempex_shell_df, inherit.aes = FALSE,
                aes(x = ShellClr, y = emmean, ymin = lower.CL, ymax = upper.CL),
                width = 0.05, linewidth = 0.8, color = "grey20") +
  geom_point(data = emmtempex_shell_df, inherit.aes = FALSE,
             aes(x = ShellClr, y = emmean),
             size = 3, color = "grey20", shape = 16) +
  scale_fill_manual(values = rep("grey80", length(unique(dftemps$ShellClr)))) +
  labs(x = "Shell Colour", y = "Shell Temperature (\u00B0C)",
       caption = "Pink: Model EMM \u00B1 95% CI") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        plot.caption = element_text(hjust = 0.5, face = "italic", color = "grey30")) +
  coord_cartesian(ylim = c(16,31.5))

Plot_A_ex <- Effectplot_ht_ex + 
  geom_signif(comparisons = list(c("Low", "High")),
              annotations = "**",   # p = 0.002
              y_position = 30, textsize = 3, vjust = -0.1) +
  labs(caption = NULL, tag = "C")

Plot_B_ex <- Effectplot_clr_ex + 
  geom_signif(comparisons = list(c("Light", "Dark")),
              annotations = "***",   # was "p < 0.001"
              y_position = 30, textsize = 3, vjust = -0.1) +
  labs(caption = NULL, tag = "D")

Plot_A_ex
Plot_B_ex


#emmeans for reporting in results section
print(emmtempin_ht_df)
print(emmtempin_shell_df)
print(emmtempex_ht_df)
print(emmtempex_shell_df)
#######################################
# ============================================================
# Load the tile temperature panel (built in a separate script)
# ============================================================
Tileboxplot <- readRDS("manufig_tiletemp.rds")

#quick check of each component panel individually:
Tileboxplot
Plot_A
Plot_B
Plot_A_ex
Plot_B_ex
TempRegression_plot
# ------------------------------------------------------------
# Re-assign tags now that tile temp is the new first row
# (A shifts everything else down one letter)
# ------------------------------------------------------------
Tileboxplot      <- Tileboxplot      + labs(tag = "A")
Plot_A           <- Plot_A           + labs(tag = "D")
Plot_B           <- Plot_B           + labs(tag = "E")
Plot_A_ex        <- Plot_A_ex        + labs(tag = "B")
Plot_B_ex        <- Plot_B_ex        + labs(tag = "C")
TempRegression_plot <- TempRegression_plot + labs(tag = "F")



# ============================================================
# Full Figure 1 assembly: Tile temp -> Internal -> External -> Regression
# ============================================================

Figure1 <- Tileboxplot /
  (Plot_A_ex + Plot_B_ex) / (Plot_A + Plot_B) /
  TempRegression_plot +
  plot_layout(heights = c(1, 1, 1, 1.2)) &
  theme(
    text = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.caption = element_text(size = 7, face = "italic", color = "grey30"),
    plot.margin = margin(t = 14, r = 6, b = 14, l = 6)   # <- THE MISSING LINE
  )

Figure1

ggsave("Figure1_ThermalCascade.tiff", Figure1,
       width = 25, height = 33, units = "cm", dpi = 600, compression = "lzw")   # 26 -> 32
saveRDS(Figure1, "manu_Figure1_full.rds")



#############################################################
# ============================================================
# EXPORT: pull Type III stats for the master results table
# Append this at the very end of the temperature script,
# after all models (including lmm_temps_plot) are fit.
# ============================================================

extract_anova <- function(model, response_name) {
  tab <- as.data.frame(anova(model, type = 3, ddf = "Satterthwaite"))
  tab$Term <- rownames(tab)
  tab %>%
    transmute(
      Response  = response_name,
      Term,
      Statistic = sprintf("F(%.0f,%.1f) = %.2f", NumDF, DenDF, `F value`),
      p.value   = `Pr(>F)`
    )
}

temp_stats <- bind_rows(
  extract_anova(lmm_tempin,      "Internal tissue temperature"),
  extract_anova(lmm_tempex,      "External shell temperature"),
  extract_anova(lmm_temps_plot,  "Internal ~ External temperature (buffering)")
)

saveRDS(temp_stats, "stats_temperature.rds")
print(temp_stats)  
