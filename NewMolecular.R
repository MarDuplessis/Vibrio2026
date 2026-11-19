library(tidyverse)
library(cocor)
library(ggplot2)
library(car)
library(emmeans)
library(patchwork)
library(lme4)
library(lmerTest)
library(scales)
library(ggsignif)
library(dplyr)

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
extract_lmer_coef <- function(model, response_name) {
  tab <- as.data.frame(summary(model)$coefficients)
  tab$Term <- rownames(tab)
  tab %>%
    filter(Term != "(Intercept)") %>%
    transmute(
      Response  = response_name,
      Term,
      Estimate  = round(Estimate, 3),
      Statistic = sprintf("t(%.1f) = %.2f", df, `t value`),
      p.value   = `Pr(>|t|)`
    )
}

# Set WD & correct contrast mode for Type III ANOVA
setwd("C:/Users/Marce/OneDrive/Desktop/491/Rwork")
options(contrasts = c("contr.sum","contr.poly"))

# Read the CSV. 
df_masternew <- read.csv("Vib16S2.csv") %>%
  mutate(
    # Ensure factors are set correctly
    Tide_ht = factor(Tide_ht, levels = c(0,1), labels = c("Low","High")),
    Tile_clr = factor(Tile_clr, levels = c(0,1), labels = c("White","Black")),
    ShellClr = factor(ShellClr, levels = c(0,1), labels = c("Light","Dark")),
    log_Vib = log10(Copy_Vib),
    log_16S = log10(Copy_16S),
    log_TLH_Vib = log10(TLH.Vib),
    log_TDH_Vib = log10(TDH.Vib),
    log_TRH_Vib = log10(TRH.Vib),
    log_Vib_16S = log10(Vib.16S)
    
  )
#Make a raw, unaltered dataframe 
df_masternewraw <- df_masternew
# 1. Impute Cq = 40 and overwrite the existing ratio columns
df_masternew <- df_masternew %>%
  mutate(
    # Impute the raw Cq non-detects to 40
    Cq_TLH = replace_na(Cq_TLH, 40),
    Cq_TDH = replace_na(Cq_TDH, 40),
    Cq_TRH = replace_na(Cq_TRH, 40),
    
    # Overwrite the existing CSV columns to eliminate the NAs
    TLH.Vib = (2^-Cq_TLH) / (1.81^-Cq_Vib),
    TDH.Vib = (2^-Cq_TDH) / (1.81^-Cq_Vib),
    TRH.Vib = (2^-Cq_TRH) / (1.81^-Cq_Vib),
    
    # Log10 transform the newly completed ratio columns for the model
    log_TLH_Vib = log10(TLH.Vib),
    log_TDH_Vib = log10(TDH.Vib),
    log_TRH_Vib = log10(TRH.Vib)
  )


#Test 1: effect on total bacteria (Response = logged copy no of total bacterial 16S)
lmm_log16S <- lmer(log_16S ~ Tide_ht * Tile_clr * ShellClr + (1|Tile_no), data = df_masternew)
#Diagnostics
qqnorm(resid(lmm_log16S))
qqline(resid(lmm_log16S))
plot(fitted(lmm_log16S), resid(lmm_log16S))
qqnorm(ranef(lmm_log16S)$Tile_no[[1]], main = "Random Effects QQ")
qqline(ranef(lmm_log16S)$Tile_no[[1]])
isSingular(lmm_log16S) #TRUE: random effect variance is zero. so all this means is you couldve gotten away without a random effect
anova(lmm_log16S, type = 3, ddf = "Satterthwaite") # no significance of anything (All F [1,18] < 2.05, All p > 0.1)
summary(lmm_log16S)
#Test 2: effect on total Vibrio (Response = logged copy no of Vibrio 16S)
lmm_logVib <- lmer(log_Vib ~ Tide_ht * Tile_clr * ShellClr + (1|Tile_no), data = df_masternew)
qqnorm(resid(lmm_logVib)) 
qqline(resid(lmm_logVib))
plot(fitted(lmm_logVib), resid(lmm_logVib))
qqnorm(ranef(lmm_logVib)$Tile_no[[1]], main = "Random Effects QQ")
qqline(ranef(lmm_logVib)$Tile_no[[1]])
isSingular(lmm_logVib)
anova(lmm_logVib, type = 3, ddf = "Satterthwaite") #Tide_ht: F [1,18] = 5.3007, p = 0.03347*, ShellClr: F [1,18] = 7.4064, p = 0.01400*

      # --- 2. GENERATE THE EMMEANS OBJECTS ---
      
      # Extract EMMs for Tidal Height main effect
      emmvib_ht <- emmeans(lmm_logVib, ~ Tide_ht, lmer.df = "satterthwaite")
      
      # Extract EMMs for Shell Colour main effect
      emmvib_shell <- emmeans(lmm_logVib, ~ ShellClr, lmer.df = "satterthwaite")
      
      
      # --- 3. VIEW THE RESULTS (Optional but recommended) ---
      # This will print the log-scale means, SE, df, and confidence intervals to your console
      summary(emmvib_ht)
      summary(emmvib_shell)
      
      # For shore height
      pairs(emmvib_ht, reverse = TRUE)
      
      # For shell colour  
      pairs(emmvib_shell, reverse = TRUE)
      confint(pairs(emmvib_ht, reverse = TRUE))
      confint(pairs(emmvib_shell, reverse = TRUE))
      
      
#######################################################
      # Inferential test - LMM without treatments (mediation model)
      lmm_logVib_temp <- lmer(log_Vib ~ Temp_in + (1|Tile_no), data = df_masternew)
      summary(lmm_logVib_temp)
      lmm_logVib_temp_stats <- anova(lmm_logVib_temp, type = 3, ddf = "Satterthwaite")
      pval_vib_lmm <- lmm_logVib_temp_stats$`Pr(>F)`[1]
      pval_vib_label <- ifelse(pval_vib_lmm < .001, "p < 0.001", paste0("p = ", round(pval_vib_lmm, 4)))
      
      # Data preparation - add group variable for colouring
      df_masternew_plot <- df_masternew %>%
        mutate(Group = interaction(Tide_ht, ShellClr, sep = " / "))
      
      # Fit simple lm — ONLY for the plotted regression line/confidence band and R²
      lm_vib_plot <- lm(log_Vib ~ Temp_in, data = df_masternew)
      r2_vib <- summary(lm_vib_plot)$r.squared
      
      # Multiple regression with treatment terms (mediation check)
      lmm_vib_multiple <- lmer(log_Vib ~ Temp_in + Tide_ht + Tile_clr + ShellClr + (1|Tile_no), 
                               data = df_masternew)
      summary(lmm_vib_multiple)
      vif(lmm_vib_multiple)
      # VIF: Temp_in = 2.56, Tide_ht = 2.48, Tile_clr = 1.07, ShellClr = 1.29 —
      # no evidence of problematic multicollinearity (all VIF < 5)
      
#######################################################

      VibTemp_plot <- ggplot(df_masternew_plot, aes(x = Temp_in, y = log_Vib)) +
        geom_smooth(method = "lm", color = "grey40", fill = "grey85",
                    linewidth = 0.8, alpha = 0.4) +
        geom_point(size = 3, alpha = 0.85, colour = "grey30") +
        annotate("text",
                 x = min(df_masternew$Temp_in, na.rm = TRUE) + 0.3,
                 y = max(df_masternew$log_Vib, na.rm = TRUE) - 0.3,
                 label = paste0("R\u00B2 = ", round(r2_vib, 2), "\n", pval_vib_label),
                 hjust = 0, vjust = 1, size = 3, color = "grey20") +
        scale_y_continuous(
          breaks = c(0, 1, 2, 3, 4),
          labels = c("1", "10", "100", "1,000", "10,000")
        ) + 
        labs(x = "Internal Tissue Temperature (\u00B0C)",
             y = "Vibrio 16S rRNA (copies/\u00B5L)",
             caption = "Grey band: 95% CI (simple linear fit). p-value from mixed model accounting for tile clustering.") +
        theme_bw(base_size = 11) +
        theme(legend.position = "none",
              plot.caption = element_text(hjust = 0.5, face = "italic",
                                          color = "grey30", size = 7))
      
      VibTemp_plot
      saveRDS(VibTemp_plot, "manu_vibtempreg.rds")   
      
      
###############################################################################

lmm_logTLHVib <- lmer(log_TLH_Vib ~ Tide_ht * Tile_clr * ShellClr + (1|Tile_no), data = df_masternew)
qqnorm(resid(lmm_logTLHVib)) 
qqline(resid(lmm_logTLHVib))
plot(fitted(lmm_logTLHVib), resid(lmm_logTLHVib))
qqnorm(ranef(lmm_logTLHVib)$Tile_no[[1]], main = "Random Effects QQ")
qqline(ranef(lmm_logTLHVib)$Tile_no[[1]])
isSingular(lmm_logTLHVib)
anova(lmm_logTLHVib, type = 3, ddf = "Satterthwaite")  # No Significance (All F < 1.8, all p > 0.2)

lmm_logTDHVib <- lmer(log_TDH_Vib ~ Tide_ht * Tile_clr * ShellClr + (1|Tile_no), data = df_masternew)
qqnorm(resid(lmm_logTDHVib)) 
qqline(resid(lmm_logTDHVib))
plot(fitted(lmm_logTDHVib), resid(lmm_logTDHVib))
qqnorm(ranef(lmm_logTDHVib)$Tile_no[[1]], main = "Random Effects QQ")
qqline(ranef(lmm_logTDHVib)$Tile_no[[1]])
isSingular(lmm_logTDHVib)
anova(lmm_logTDHVib, type = 3, ddf = "Satterthwaite") #No Significance (All F < 1.75, all p > 0.2)

lmm_logTRHVib <- lmer(log_TRH_Vib ~ Tide_ht * Tile_clr * ShellClr + (1|Tile_no), data = df_masternew)
qqnorm(resid(lmm_logTRHVib)) 
qqline(resid(lmm_logTRHVib))
plot(fitted(lmm_logTRHVib), resid(lmm_logTRHVib))
qqnorm(ranef(lmm_logTRHVib)$Tile_no[[1]], main = "Random Effects QQ")
qqline(ranef(lmm_logTRHVib)$Tile_no[[1]])
isSingular(lmm_logTRHVib)
anova(lmm_logTRHVib, type = 3, ddf = "Satterthwaite") #No significance (All F < 2.38, all p > 0.1)

#Stats export, MINUS multiple regression


molecular_stats <- bind_rows(
  extract_anova(lmm_log16S,     "Universal bacterial 16S rRNA"),
  extract_anova(lmm_logVib,     "Vibrio-specific 16S rRNA"),
  extract_anova(lmm_logVib_temp, "Vibrio-specific 16S rRNA ~ Internal temperature"),
  extract_anova(lmm_logTLHVib,  "tlh:Vibrio 16S ratio"),
  extract_anova(lmm_logTDHVib,  "tdh:Vibrio 16S ratio"),
  extract_anova(lmm_logTRHVib,  "trh:Vibrio 16S ratio")
)

saveRDS(molecular_stats, "stats_molecular.rds")
print(molecular_stats)

#Mediation stats export (multiple regression)
mediation_stats <- bind_rows(
  extract_lmer_coef(lmm_vib_multiple, "Vibrio ~ Temperature + Treatments")
)

saveRDS(mediation_stats, "stats_mediation.rds")
print(mediation_stats)


# --- 1. PREP THE DATA & BACK-TRANSFORM ---
# Back-transform the Model EMMs for Tidal Height (10^x)
# Replace 'emmvib_ht' with the actual name of your emmeans object
emm_vib_ht_df <- as.data.frame(emmvib_ht) %>%
  mutate(
    emmean_back = 10^emmean,
    lower.CL_back = 10^lower.CL,
    upper.CL_back = 10^upper.CL
  )
# Back-transform the Model EMMs for Shell Colour (10^x)
# Replace 'emmvib_shell' with the actual name of your emmeans object
emm_vib_shell_df <- as.data.frame(emmvib_shell) %>%
  mutate(
    emmean_back = 10^emmean,
    lower.CL_back = 10^lower.CL,
    upper.CL_back = 10^upper.CL
  )


# --- 2. BUILD THE PLOTS ---

# Panel A: Tidal Height
Effectplot_vib_ht <- ggplot(df_masternew, aes(x = Tide_ht, y = Copy_Vib, fill = Tide_ht)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_jitter(width = 0.05, alpha = 0.3, color = "grey35") + 
  geom_errorbar(data = emm_vib_ht_df, inherit.aes = FALSE,
                aes(x = Tide_ht, y = emmean_back, ymin = lower.CL_back, ymax = upper.CL_back), 
                width = 0.05, linewidth = 0.8, color = "grey20") +
  geom_point(data = emm_vib_ht_df, inherit.aes = FALSE,
             aes(x = Tide_ht, y = emmean_back), 
             size = 3, color = "grey20", shape = 16) +
  scale_y_log10(labels = label_comma(), limits = c(0.1, 25000)) + 
  scale_fill_manual(values = rep("grey80", length(unique(df_masternew$Tide_ht)))) +
  labs(x = "Shore Height", y = "Vibrio 16S rRNA (copies/\u00B5L)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# Panel B: Shell Colour
Effectplot_vib_clr <- ggplot(df_masternew, aes(x = ShellClr, y = Copy_Vib, fill = ShellClr)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_jitter(width = 0.05, alpha = 0.3, color = "grey35") + 
  geom_errorbar(data = emm_vib_shell_df, inherit.aes = FALSE,
                aes(x = ShellClr, y = emmean_back, ymin = lower.CL_back, ymax = upper.CL_back), 
                width = 0.05, linewidth = 0.8, color = "grey20") +
  geom_point(data = emm_vib_shell_df, inherit.aes = FALSE,
             aes(x = ShellClr, y = emmean_back), 
             size = 3, color = "grey20", shape = 16) +
  scale_y_log10(labels = label_comma(), limits = c(0.1, 25000)) + 
  scale_fill_manual(values = rep("grey80", length(unique(df_masternew$ShellClr)))) +
  labs(x = "Shell Colour", y = "Vibrio 16S rRNA (copies/\u00B5L)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# --- 3. ADD SIGNIFICANCE BRACKETS & GLUE TOGETHER ---

Plot_Vib_A <- Effectplot_vib_ht + 
  geom_signif(
    comparisons = list(c("Low", "High")), 
    annotations = "*",              #
    y_position = 3.9,
    textsize = 3, vjust = -0.1
  ) +
  labs(tag = "A")

Plot_Vib_B <- Effectplot_vib_clr + 
  geom_signif(
    comparisons = list(c("Light", "Dark")), 
    annotations = "*",              #
    y_position = 3.9,
    textsize = 3, vjust = -0.1
  ) +
  labs(tag = "B")

# Combine with Patchwork
# Load/confirm panels are ready (already built above in this script)
# Plot_Vib_A, Plot_Vib_B  -> Vibrio 16S by treatment (A, B)
# VibTemp_plot            -> Vibrio 16S ~ internal temperature regression (C)

# Strip any leftover in-panel captions (these move to the master caption instead)
Plot_Vib_A <- Plot_Vib_A + labs(caption = NULL)
Plot_Vib_B <- Plot_Vib_B + labs(caption = NULL)
VibTemp_plot <- VibTemp_plot + labs(caption = NULL, tag = "C")

Figure2 <- (Plot_Vib_A + Plot_Vib_B) /
  VibTemp_plot +
  plot_layout(heights = c(1, 1.2)) &
  theme(
    text = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.margin = margin(t = 14, r = 6, b = 14, l = 6)
  )

Figure2

ggsave("Figure2_VibrioAbundance.tiff", Figure2,
       width = 18, height = 20, units = "cm", dpi = 600, compression = "lzw")

saveRDS(Figure2, "manu_Figure2_full.rds")




