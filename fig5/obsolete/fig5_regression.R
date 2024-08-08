library(readr)
library(dplyr)
library(ggplot2)

df <- read_tsv("~/bioinf/msms/res_ratio_hc.tsv")

# convert both intensity ratio and conc ration in log10 scale
df$log_intensity_ratio <- log10(df$intensity_ratio)
df$log_conc_ratio <- log10(df$conc_ratio)

df_median <- df %>% 
  group_by(Protease, run1, match_ig_type, log_conc_ratio) %>%
  summarise(med = median(log_intensity_ratio), 
            N = n()) 
  

coefs <- df_median %>%
  group_by(Protease, run1, match_ig_type) %>%
  do(regr = lm(med ~ log_conc_ratio, weights = N, .)) %>%
  summarise(Protease = Protease, 
            run1 = run1, 
            match_ig_type = match_ig_type,
            intercept = coef(regr)[1], 
            slope = coef(regr)[2], 
            r2 = summary(regr)$r.squared)


coefs$eq <- paste0("italic(y) == ", format(coefs$intercept, digits = 2), 
                   " + ", format(coefs$slope, digits = 2), " %.% italic(x)")

coefs$r2 <- paste0("italic(R)^2 == ", format(coefs$r2, digits = 2))

ggplot(df, aes(x = log_conc_ratio, y = log_intensity_ratio, color = as.factor(log_conc_ratio))) +
  facet_grid(match_ig_type ~ Protease+run1) +
  geom_boxplot(outlier.size = 0.3) +
  geom_abline(data = coefs, aes(intercept = intercept, slope = slope), color = "black", linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, color = "black") + # ground truth
  geom_text(data = coefs, aes(0, -4, label = eq), parse = TRUE, size = 1.5, color = "black", hjust = 0) +
  geom_text(data = coefs, aes(0, -4.5, label = r2), parse = TRUE, size = 1.5, color = "black", hjust = 0) +
  theme_bw() +
  coord_fixed(ratio=1) 
  


