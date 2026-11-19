#Necessary Packages
library(tidyverse)
library(cocor)
library(ggplot2)
library(car)
library(emmeans)
library(patchwork)
library(lme4)
library(lmerTest)
library(ARTool)
library(scales)
library(survival)
library(logistf)

setwd("C:/Users/Marce/OneDrive/Desktop/491/Rwork")
options(contrasts = c("contr.sum","contr.poly"))


dfMPN <- read.csv("MPN.csv") %>%
  na.omit %>%
  mutate(
    # 1101 = confirmed censored (no upper conf. limit in your raw data)
    # 1100 = has conf. limits (90-4100), so treat as exact estimate
    censored = ifelse(MPN.g == 1101, 1, 0),
    logMPN.g = log10(MPN.g),
    Tide_ht  = factor(Tide_ht,  levels = c(0,1), labels = c("Low-Shore","High-Shore")),
    Tile_clr = factor(Tile_clr, levels = c(0,1), labels = c("White","Black")),
    Shellclr = factor(Shellclr, levels = c(0,1), labels = c("Light","Dark"))
  )


# Asks: do treatments predict whether a sample hits the detection limit?
# Used because: heavy censoring (83%) causes separation in standard logistic
firth_detect <- logistf(censored ~ Tide_ht * Tile_clr * Shellclr,
                        data = dfMPN,
                        firth = TRUE)
summary(firth_detect)
table(dfMPN$Tide_ht, dfMPN$Tile_clr)

lrt_chisq <- 2 * (firth_detect$loglik["full"] - firth_detect$loglik["null"])
lrt_df    <- firth_detect$df
lrt_p     <- 1 - pchisq(lrt_chisq, lrt_df)

mpn_stats <- tibble::tibble(
  Response  = "Culturable Vibrio (MPN assay saturation)",
  Term      = "Overall model (Firth LRT)",
  Statistic = sprintf("\u03c7\u00b2(%d) = %.2f", lrt_df, lrt_chisq),
  p.value   = lrt_p
)

saveRDS(mpn_stats, "stats_mpn_supplementary.rds")
print(mpn_stats)

# --- MODEL 2: Tobit Regression (SECONDARY) ---
# Asks: do treatments predict MPN magnitude, accounting for censoring properly?
# Used because: formally models censored observations rather than treating as exact
tobit_MPN <- survreg(
  Surv(logMPN.g, censored, type = "right") ~ Tide_ht * Tile_clr * Shellclr,
  dist = "gaussian",
  data = dfMPN
)
summary(tobit_MPN)

# Main effects only - no interactions
tobit_reduced <- survreg(
  Surv(logMPN.g, censored, type = "right") ~ Tide_ht + Tile_clr + Shellclr,
  dist = "gaussian",
  data = dfMPN
)
summary(tobit_reduced)

# --- MODEL 3: Log-linear LM (SUPPLEMENTARY) ---
# Asks: same question as Tobit but ignoring censoring
# Used because: comparable to prior literature, allows emmeans back-calculation
# Caveat: treats censored values as exact — interpret cautiously
lm_logMPN <- lm(logMPN.g ~ Tide_ht * Tile_clr * Shellclr, data = dfMPN)

# Diagnostics
par(mfrow=c(1,2))
qqnorm(resid(lm_logMPN)); qqline(resid(lm_logMPN))
plot(lm_logMPN, which=1)
par(mfrow=c(1,1))

Anova(lm_logMPN, type = 3) 

# EMmeans for plotting (back-transformed to raw MPN scale)
emmMPN_df <- as.data.frame(emmeans(lm_logMPN, ~ Tide_ht * Shellclr)) %>%
  mutate(
    MPN    = 10^emmean,
    lower  = 10^lower.CL,
    upper  = 10^upper.CL
  )

glm_detect <- glm(censored ~ Tide_ht * Tile_clr * Shellclr,
                  family = binomial(link = "logit"),
                  data = dfMPN)
Anova(glm_detect, type = 3)


glm_gamma <- glm(MPN.g ~ Tide_ht * Tile_clr * Shellclr,
                 family = Gamma(link = "log"),
                 data = dfMPN)
qqnorm(resid(glm_gamma))

Anova(glm_gamma, type = 3)


lm_logMPN <- lm(logMPN.g ~ Tide_ht * Tile_clr * Shellclr, data = dfMPN)

qqnorm(resid(lm_logMPN))
qqline(resid(lm_logMPN))
Anova(lm_logMPN, type = 3) #tide height x shell colour interaction (p = 0.003348)
emmlogMPN <- emmeans(lm_logMPN, ~ Tide_ht * Tile_clr * Shellclr)
pairs(emmlogMPN, adjust = "holm")
emmlogMPN_df <- as.data.frame(emmlogMPN)
emmlogMPN_df$MPN <- 10^emmlogMPN_df$emmean
emmlogMPN_df$lower <- 10^emmlogMPN_df$lower.CL
emmlogMPN_df$upper <- 10^emmlogMPN_df$upper.CL

emmMPN <- emmeans(lm_logMPN, ~ Tide_ht * Tile_clr * Shellclr)
summary(emmMPN, type = "response")
pairs(emmMPN, adjust = "holm")
emmMPN_df <- as.data.frame(emmMPN)
emmMPN_df$MPN <- 10^emmMPN_df$emmean
emmMPN_df$lower <- 10^emmMPN_df$lower.CL
emmMPN_df$upper <- 10^emmMPN_df$upper.CL






# 1. Back-calculate the raw assay ceiling (10^x)
raw_upper_limit <- 10^max(dfMPN$logMPN.g, na.rm = TRUE)

# 2. Set the height for the letters 
# We multiply by 1.5 to push the letters visually higher on the log scale
letter_height <- raw_upper_limit * 1.5  

# 4. Generate the Plot
MPNPLOT <- ggplot(dfMPN %>% filter(!is.na(Shellclr)), 
       aes(x = Shellclr, y = 10^logMPN.g, fill = Tile_clr)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), 
              alpha = 0.5) +
  geom_hline(yintercept = raw_upper_limit, linetype = "dashed",
             color = "#CC0000", linewidth = 1) +
  facet_wrap(~ Tide_ht) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = c("White" = "grey95", "Black" = "grey30")) +
  labs(x = "Shell Colour", y = "Presumptive Culturable Vibrio(MPN/g)",
       fill = "Tile Colour",
       caption = "Red dashed line: Assay Upper Limit of Quantification") +
  theme_bw(base_size = 14) +
  theme(plot.caption = element_text(color = "#CC0000", face = "italic", hjust = 0.5))
MPNPLOT



