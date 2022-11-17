Evolbarriers
================
RP Bhatia
16 September 2022

# Data analysis performed in paper "Evolutionary barriers to horizontal gene transfer in macrophage associated Salmonella" doi: <https://doi.org/10.1101/2022.04.01.486712>

# library load

``` r
library(reshape)
library(reshape2)
```

``` r
library(plyr)
```

``` r
library(dplyr)
```

``` r
library(broom)
library(tidyr)
```

``` r
library(gapminder)
library(tidyverse)
```

``` r
library(ggplot2)
library(ggforce)
```

``` r
library(gridExtra)
```

``` r
library(ggpubr)
```

``` r
library(viridis)
```

``` r
library(ggsignif)
```

``` r
library(hrbrthemes)
library(knitr)
```

``` r
library(wesanderson)
library(RColorBrewer)
library(paletteer)
```

``` r
library(scales)
```

``` r
library(tsiMisc)
library(MASS)
```

``` r
library(EnvStats)
```

``` r
library(goft)
```

``` r
library(paletteer)
library(awtools)
library(emmeans)
```

# Mandel's test for Linearity (testing effect of protein dosage)

# load data obtained from GrowthProfiler960

``` r
df = read.csv("inspi2_mlm_data.csv", header = T, stringsAsFactors = T)
head(df)
```

    ##     Genes   aTc        W
    ## 1 b0119_1 aTc_0 1.095270
    ## 2 b0119_1 aTc_2 1.215448
    ## 3 b0119_1 aTc_4 1.085465
    ## 4 b0119_1 aTc_6 1.141068
    ## 5 b0119_1 aTc_8 1.079213
    ## 6 b0119_2 aTc_0 1.101868

# arrange data

``` r
df$Genes <- gsub("_1", "", df$Genes)
df$Genes <- gsub("_2", "", df$Genes)
df$Genes <- gsub("_3", "", df$Genes)
df$Genes <- gsub("_4", "", df$Genes)


df_mlm = df %>% dplyr::group_by(Genes, aTc) %>% 
  dplyr::summarise_all(purrr::discard, is.na)

df_mlm$aTc <- gsub("aTc_", "", df_mlm$aTc)
head(df_mlm)
```

    ## # A tibble: 6 x 3
    ## # Groups:   Genes, aTc [2]
    ##   Genes aTc       W
    ##   <chr> <chr> <dbl>
    ## 1 b0119 0     1.10 
    ## 2 b0119 0     1.10 
    ## 3 b0119 0     0.932
    ## 4 b0119 0     0.923
    ## 5 b0119 2     1.22 
    ## 6 b0119 2     1.07

# perform mandel's test (code is run using data for one infection relevant environment (InSPI2))

``` r
df_mlm$aTc <- as.numeric(df_mlm$aTc)
bgene <- unique(df_mlm$Genes)

res_SE_lm = rep(0,44)
res_SE_poly = rep(0,44)
df_lm = rep(0,44)
df_poly = rep(0,44)
n = rep(0,44)
DSsq = rep(0,44)
TV = rep(0,44)
F_stat = rep(0,44)
p_val = rep(0,44)
cal_func = rep(0,44)

for (i in 1:length(bgene)){
  current_bgene <- bgene[i]
  df <- filter(df_mlm, Genes==current_bgene) # subset the dataframe
  fit_l <- lm(W ~ aTc, data = df)
  fit_p <- lm(W ~ aTc + I(aTc^2), data = df)
  res_SE_lm[i] = summary(fit_l)$sigma
  res_SE_poly[i] = summary(fit_p)$sigma
  df_lm[i] = summary(fit_l)$df[2]
  df_poly[i] = summary(fit_p)$df[2]
  n[i] = nrow(df)
  DSsq[i] = (df_lm[i])*res_SE_lm[i]^2 - (df_poly[i])*res_SE_poly[i]^2
  TV[i] = DSsq[i]/res_SE_poly[i]^2
  F_stat[i] = qf(0.95, 1, df_poly[i])
  p_val[i] = pf(TV[i], 1, (n[i]-3), lower.tail = F)
  cal_func[i] = ifelse(TV[i] < F_stat[i], "linear", ifelse(TV[i] > F_stat[i], "polynomial", NA))
  dat = matrix(nrow = 45, ncol= 6)
  dat[,1] = c("Gene", bgene)
  dat[,2] = c("DSsq", DSsq)
  dat[,3] = c("TV", TV)
  dat[,4] = c("F_statistic", F_stat)
  dat[,5] = c("p_val", p_val)
  dat[,6] = c("best_fit", cal_func)
  filename <- paste("mandels_result_inspi2.csv")
  write.table(dat,filename,row.names= FALSE, col.names= FALSE, 
              sep=",", quote = FALSE)
}
head(knitr::kable(dat))
```

    ## [1] "|      |                     |                     |                 |                     |           |"
    ## [2] "|:-----|:--------------------|:--------------------|:----------------|:--------------------|:----------|"
    ## [3] "|Gene  |DSsq                 |TV                   |F_statistic      |p_val                |best_fit   |"
    ## [4] "|b0119 |0.000144477153556311 |0.0146190304201364   |4.45132177246813 |0.90518056369526     |linear     |"
    ## [5] "|b0127 |0.00774982220670625  |0.998214468751948    |4.45132177246813 |0.331752782216723    |linear     |"
    ## [6] "|b0179 |0.00114628565136329  |0.123415452100135    |4.45132177246813 |0.729675561329717    |linear     |"

# plot data

``` r
lm_plt = list()
poly_plt = list()

for (i in 1:length(bgene)){
  current_bgene <- bgene[i]
  df <- filter(df_mlm, Genes==current_bgene)
  l<- ggplot(df, aes(y = W, x = aTc)) + 
    geom_point(shape=16, color="springgreen") +
    geom_smooth(method="lm", formula= y ~ x , se = FALSE,
                size = 1, colour = "Black") +
    xlab("aTc (ng/ml)") + 
    ylab("Relative_fitness(w)") +
    scale_y_continuous(limits = c(0, 3.0), breaks = seq(0, 3.0, by = 0.2)) +
    ggtitle(paste(bgene[i]))
  lm_plt[[i]] = l
  p<- ggplot(df, aes(y = W, x = aTc)) + 
    geom_point(shape=16, color="springgreen") +
    geom_smooth(method="lm", formula= y ~ x + I(x^2) , se = FALSE,
                size = 1, colour = "Black") +
    xlab("aTc (ng/ml)") + 
    ylab("Relative_fitness(w)") +
    scale_y_continuous(limits = c(0, 3.0), breaks = seq(0, 3.0, by = 0.2)) +
    ggtitle(paste(bgene[i]))
  poly_plt[[i]] = p
}

ggsave(
  filename = "lm_plots_inspi2.pdf", 
  plot = marrangeGrob(lm_plt, nrow=3, ncol=3), 
  width = 15, height = 9)

ggsave(
  filename = "poly_plots_inspi2.pdf", 
  plot = marrangeGrob(poly_plt, nrow=3, ncol=3), 
  width = 15, height = 9)
```

# plotting gene dosage data

``` r
hm = read.csv("genedosage_data.csv", header = T)

hmp = ggplot(hm, aes(x = factor(env, level = c('InSPI2', 'InSPI2 Hypoxic', 'InSPI2 Ciprofloxacin', 'InSPI2 Low Mg')), y=gene)) + 
  geom_tile(aes(fill=dd), colour="black", size=1) + 
  coord_equal() + 
  theme(text = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(colour = "white", fill = "white",
                                        size = 1, linetype = "solid"),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.1)) + 
  theme(legend.position = "bottom") +
  scale_x_discrete(position = 'top') +
  scale_fill_manual(values = c("darkslategray3", "dodgerblue4"),
                    name = "", labels = c("No effect of dosage", "Dosage effect")) +
  xlab("") +
  ylab("Dosage dependent genes")

ggsave("genedosage.png", hmp, height = 7, width = 7)
hmp
```

![](HGT_barriers_files/figure-markdown_github/unnamed-chunk-6-1.png)

# load data obtained from NGS

``` r
df1 = read.csv("fitness_deepseq.csv", header = T, stringsAsFactors = T)

head(df)
```

    ## # A tibble: 6 x 3
    ## # Groups:   Genes, aTc [2]
    ##   Genes   aTc     W
    ##   <chr> <dbl> <dbl>
    ## 1 b4373     0 1.14 
    ## 2 b4373     0 1.14 
    ## 3 b4373     0 0.976
    ## 4 b4373     0 0.970
    ## 5 b4373     2 1.21 
    ## 6 b4373     2 1.16

``` r
df_dpsq1 = df1 %>% dplyr::group_by(gene,Env) %>% 
  dplyr::summarise_all(purrr::discard, is.na)

head(df_dpsq1)
```

    ## # A tibble: 6 x 10
    ## # Groups:   gene [2]
    ##   gene  Env            rel_fitness   PPI    RI    MI    bp F_cat  GC_DEV FOP_DEV
    ##   <fct> <fct>                <dbl> <int> <int> <int> <int> <fct>   <dbl>   <dbl>
    ## 1 b0119 InSPI2               0.967    50     2     0   358 N     0.00132  0.0744
    ## 2 b0119 InSPI2_cip           0.801    50     2     0   358 N     0.00132  0.0744
    ## 3 b0119 InSPI2_hypoxic       0.989    50     2     0   358 N     0.00132  0.0744
    ## 4 b0119 InSPI2_lowmg         0.919    50     2     0   358 N     0.00132  0.0744
    ## 5 b0127 InSPI2               0.919     0     1     0   923 O     0.00161  0.0194
    ## 6 b0127 InSPI2_cip           0.788     0     1     0   923 O     0.00161  0.0194

# Run linear regression for all environments and extract p-values (except for gene categories) (code is running regression for only one variable. Variables can be changed in the model in line 169)

``` r
model.lm <- unlist(lapply(split(df_dpsq1,
                                list(df_dpsq1$Env)), function(chunk)
                               {return(summary(lm(rel_fitness ~ PPI, 
                                                  data=chunk))$coefficients[2,4])}))

mo = data.frame(model.lm)
padj = p.adjust(mo$model.lm, method = "fdr", n = length(mo$model.lm))

knitr::kable(padj)
```

|          x|
|----------:|
|  0.7540117|
|  0.7540117|
|  0.7540117|
|  0.7540117|

# Wilcoxon Rank Sum Test for effect of gene categories

``` r
df_fcat = read.csv("fitness_deepseq.csv", header = T, stringsAsFactors = T)

df_fcat = df_fcat %>% dplyr::group_by(gene,Env) %>% 
  dplyr::summarise_all(purrr::discard, is.na)


df_fcat$F_cat = as.factor(df_fcat$F_cat)
df_fcat <- df_fcat[df_fcat$F_cat!="N",]
df_fcat <- df_fcat[df_fcat$F_cat!="IO",]
df_fcat <- df_fcat[df_fcat$F_cat!="",]
levels(df_fcat$F_cat)
```

    ## [1] "I"  "IO" "N"  "O"

``` r
df_fcat$F_cat <- factor(df_fcat$F_cat) # get rid of it...
levels(df_fcat$F_cat)
```

    ## [1] "I" "O"

``` r
head(df_fcat)
```

    ## # A tibble: 6 x 10
    ## # Groups:   gene [2]
    ##   gene  Env            rel_fitness   PPI    RI    MI    bp F_cat  GC_DEV FOP_DEV
    ##   <fct> <fct>                <dbl> <int> <int> <int> <int> <fct>   <dbl>   <dbl>
    ## 1 b0127 InSPI2               0.919     0     1     0   923 O     0.00161  0.0194
    ## 2 b0127 InSPI2_cip           0.788     0     1     0   923 O     0.00161  0.0194
    ## 3 b0127 InSPI2_hypoxic       0.803     0     1     0   923 O     0.00161  0.0194
    ## 4 b0127 InSPI2_lowmg         1.06      0     1     0   923 O     0.00161  0.0194
    ## 5 b0179 InSPI2               0.772    46     0     8  1022 O     0.00826  0.0292
    ## 6 b0179 InSPI2_cip           0.718    46     0     8  1022 O     0.00826  0.0292

``` r
fcat_env = filter(df_fcat, Env == "InSPI2") # data can be filtered using different environments
Inform = filter(fcat_env, F_cat == "I")
Oper = filter(fcat_env, F_cat == "O")

head(fcat_env)
```

    ## # A tibble: 6 x 10
    ## # Groups:   gene [6]
    ##   gene  Env    rel_fitness   PPI    RI    MI    bp F_cat  GC_DEV FOP_DEV
    ##   <fct> <fct>        <dbl> <int> <int> <int> <int> <fct>   <dbl>   <dbl>
    ## 1 b0127 InSPI2       0.919     0     1     0   923 O     0.00161 0.0194 
    ## 2 b0179 InSPI2       0.772    46     0     8  1022 O     0.00826 0.0292 
    ## 3 b0215 InSPI2       0.991    26     0     0   728 I     0.0102  0.00820
    ## 4 b0440 InSPI2       0.877    42     2     0   272 I     0.00545 0.0110 
    ## 5 b0607 InSPI2       0.981    15     0     0   425 O     0.0221  0.0699 
    ## 6 b0623 InSPI2       0.795    30     1     0   206 I     0.00208 0.0714

``` r
wrs.test = wilcox.test(Inform$rel_fitness,Oper$rel_fitness,
                       alternative = "less", paired = F)

#extract median values
summary(Inform$rel_fitness)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.7919  0.8863  0.8045  0.9849  1.0744

``` r
summary(Oper$rel_fitness)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.5506  0.7625  0.9130  0.8634  0.9728  1.0490

# plotting relation between fitness and gene interactions

``` r
#format labels
df_dpsq1$Env = factor(df_dpsq1$Env, levels = c("InSPI2", "InSPI2_hypoxic", 
                                           "InSPI2_cip", "InSPI2_lowmg"), labels = c("InSPI2", "InSPI2 Hypoxic","InSPI2 Ciprofloxacin", "InSPI2 Low Mg"))

#plot
ppi_labels = data.frame(Env = c("InSPI2", "InSPI2 Hypoxic","InSPI2 Ciprofloxacin", "InSPI2 Low Mg"), label = c("R\u00B2 = 0.002 p = 0.75", "R\u00B2 = 0.004 p = 0.75", "R\u00B2 = 0.008 p = 0.75", "R\u00B2 = 0.002 p = 0.75"))

ppi.lm <- lm(rel_fitness ~ PPI, df_dpsq1)

fp = ggplot(df_dpsq1, aes(x = PPI, y = rel_fitness)) +
  geom_point(size = 1, colour = "springgreen4", shape = 18) +
  facet_grid(~ Env) +
  xlab("Protein-protein interactions") +
  ylab("Mean Relative Fitness (w)") +
  geom_hline(yintercept = 1, lty = 3, colour = "grey54") +
  geom_abline(slope = coef(ppi.lm)[["PPI"]], 
              intercept = coef(ppi.lm)[["(Intercept)"]]) +
  theme(text = element_text(size = 9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_text(x = 35, y = 1.12, aes(label = label), data = ppi_labels, size = 3)

ri_labels = data.frame(Env = c("InSPI2", "InSPI2 Hypoxic","InSPI2 Ciprofloxacin", "InSPI2 Low Mg"), label = c("R\u00B2 = 0.010 p = 0.70", "R\u00B2 = 0.002 p = 0.73", "R\u00B2 = 0.011 p = 0.70", "R\u00B2 = 0.009 p = 0.70"))

ri.lm <- lm(rel_fitness ~ RI, df_dpsq1)

fr = ggplot(df_dpsq1, aes(x = RI, y = rel_fitness)) +
  geom_point(size = 1, colour = "springgreen4", shape = 18) +
  geom_abline(slope = coef(ri.lm)[["RI"]], 
              intercept = coef(ri.lm)[["(Intercept)"]]) +
  facet_grid(~ Env) +
  xlab("Regulatory interactions") +
  ylab("Mean Relative Fitness (w)") +
  geom_hline(yintercept = 1, lty = 3, colour = "grey54") +
  theme(text = element_text(size = 9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_text(x = 5, y = 1.12, aes(label = label), data = ri_labels, size = 3)

mi_labels = data.frame(Env = c("InSPI2", "InSPI2 Hypoxic","InSPI2 Ciprofloxacin", "InSPI2 Low Mg"), label = c("R\u00B2 = 0.111 p = 0.05", "R\u00B2 = 0.09 p = 0.05", "R\u00B2 = 0.038 p = 0.20", "R\u00B2 = 0.097 p = 0.05"))

mi.lm <- lm(rel_fitness ~ MI, df_dpsq1)

fm = ggplot(df_dpsq1, aes(x = MI, y = rel_fitness)) +
  geom_point(size = 1, colour = "springgreen4", shape = 18) +
  geom_abline(slope = coef(mi.lm)[["MI"]], 
              intercept = coef(mi.lm)[["(Intercept)"]]) +
  facet_grid(~ Env) +
  xlab("Metabolic interactions") +
  ylab("Mean Relative Fitness (w)") +
  geom_hline(yintercept = 1, lty = 3, colour = "grey54") +
  theme(text = element_text(size = 9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_text(x = 12, y = 1.12, aes(label = label), data = mi_labels, size = 3)

fg = grid.arrange(fp,fr,fm, heights =c(9,9,9))
```

![](HGT_barriers_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
ggsave("mean_fitness_vs_interactions.png", fg, height = 13, width = 18, units = "cm")
grid.arrange(fg)
```

![](HGT_barriers_files/figure-markdown_github/unnamed-chunk-10-2.png)

# plotting relation between fitness and gene length

``` r
bp_labels = data.frame(Env = c("InSPI2", "InSPI2 Hypoxic","InSPI2 Ciprofloxacin", "InSPI2 Low Mg"), label = c("R\u00B2 = 0.25 p = 6.17e-04", "R\u00B2 = 0.3 p = 1.29e-04", "R\u00B2 = 0.13 p = 0.016", "R\u00B2 = 0.13 p = 0.019"))

bp.lm <- lm(rel_fitness ~ bp, df_dpsq1)

f = ggplot(df_dpsq1, aes(x = bp, y = rel_fitness)) +
  geom_point(size = 1, colour = "springgreen4", shape = 18) +
  geom_abline(slope = coef(bp.lm)[["bp"]], 
              intercept = coef(bp.lm)[["(Intercept)"]]) + 
  facet_wrap(~ Env) +
  xlab("Gene Length (bp)") +
  ylab("Mean Relative Fitness (w)") +
  geom_hline(yintercept = 1, lty = 3, colour = "grey54") +
  theme(text = element_text(size = 9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0, 3000), breaks = seq(0, 3000, by = 500)) +
  geom_text(x = 1000, y = 1.20, aes(label = label), data = bp_labels, size = 3)

ggsave("mean_fitness_vs_bp.png", f, height = 15, width = 15, units = "cm")
f
```

![](HGT_barriers_files/figure-markdown_github/unnamed-chunk-11-1.png)

# plotting relation between fitness and GC content and, codon usage

``` r
gc_labels = data.frame(Env = c("InSPI2", "InSPI2 Hypoxic","InSPI2 Ciprofloxacin", "InSPI2 Low Mg"), label = c("R\u00B2 = 0.004 p = 0.79", "R\u00B2 = 0.005 p = 0.79", "R\u00B2 = 0.037 p = 0.79", "R\u00B2 = 0.00 p = 0.79"))

gc.lm <- lm(rel_fitness ~ GC_DEV, df_dpsq1)

fgc = ggplot(df_dpsq1, aes(x = GC_DEV, y = rel_fitness)) +
  geom_point(size = 1, colour = "springgreen4", shape = 18) +
  geom_abline(slope = coef(gc.lm)[["GC_DEV"]], 
              intercept = coef(gc.lm)[["(Intercept)"]]) + 
  facet_grid(~ Env) +
  xlab("Absolute Deviation in GC Content") +
  ylab("Mean Relative Fitness (w)") +
  geom_hline(yintercept = 1, lty = 3, colour = "grey54") +
  theme(text = element_text(size = 9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_text(x = 0.02, y = 1.20, aes(label = label), data = gc_labels, size = 3)


cu_labels = data.frame(Env = c("InSPI2", "InSPI2 Hypoxic","InSPI2 Ciprofloxacin", "InSPI2 Low Mg"), label = c("R\u00B2 = 0.049 p = 0.29", "R\u00B2 = 0.056 p = 0.44", "R\u00B2 = 0.014 p = 0.29", "R\u00B2 = 0.022 p = 0.44"))

cu.lm <- lm(rel_fitness ~ FOP_DEV, df_dpsq1)

fcu = ggplot(df_dpsq1, aes(x = FOP_DEV, y = rel_fitness)) +
  geom_point(size = 1, colour = "springgreen4", shape = 18) +
  geom_abline(slope = coef(cu.lm)[["FOP_DEV"]], 
              intercept = coef(cu.lm)[["(Intercept)"]]) +
  facet_grid(~ Env) +
  xlab("Absolute Deviation in Codon Usage (FOP)") +
  ylab("Mean Relative Fitness (w)") +
  geom_hline(yintercept = 1, lty = 3, colour = "grey54") +
  theme(text = element_text(size = 9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_text(x = 0.05, y = 1.20, aes(label = label), data = cu_labels, size = 3)

a = grid.arrange(fgc,fcu, heights =c(6,6))
```

![](HGT_barriers_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
ggsave("mean_fitness_vs_gc_fop.png", a, height = 13, width = 18, units = "cm")
grid.arrange(a)
```

![](HGT_barriers_files/figure-markdown_github/unnamed-chunk-12-2.png)

# boxplot represenation for effect of gene categories on fitness

``` r
plt_data = read.csv("deepseq_ge_data.csv", header = T, stringsAsFactors = T, check.names = F)

plt_df = plt_data %>%
  pivot_longer(cols = starts_with("inspi2"),
               names_to = "Env",
               values_to = "w",
               values_drop_na = TRUE)

plt_df$Functional_category = as.factor(plt_df$Functional_category)
plt_df <- plt_df[plt_df$Functional_category!="N",]
plt_df <- plt_df[plt_df$Functional_category!="IO",]
levels(plt_df$Functional_category)
```

    ## [1] "I"  "IO" "N"  "O"

``` r
plt_df$Functional_category <- factor(plt_df$Functional_category) # get rid of it...
levels(plt_df$Functional_category)
```

    ## [1] "I" "O"

``` r
plt_df$Env = factor(plt_df$Env, levels = c("inspi2", "inspi2hypoxic", 
                                           "inspi2cip", "inspi2lowmg"), labels = c("InSPI2", "InSPI2 Hypoxic","InSPI2 Ciprofloxacin", "InSPI2 Low Mg"))


b <- ggplot(plt_df) +
  aes(x = Functional_category, y = w, color = Functional_category) +
  geom_boxplot(width = 0.5) +
  scale_color_paletteer_d("yarrr::info") +
  geom_jitter(size=1, alpha=0.9, position = position_jitter(width = .09)) +
  facet_wrap(~ Env) +
  xlab("Functional category") +
  ylab("Mean Relative fitness (w)") +
  scale_x_discrete(labels = c("Informational", "Operational")) +
  theme(legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 1, lty=2) +
  scale_y_continuous(limits = c(0, 1.4), breaks = seq(0, 1.4, by = 0.2))

ggsave("mean_fitness_vs_funct_cat.png", b, height = 6, width = 6)
```

    ## Warning: Removed 11 rows containing missing values (geom_point).

``` r
b
```

    ## Warning: Removed 11 rows containing missing values (geom_point).

![](HGT_barriers_files/figure-markdown_github/unnamed-chunk-13-1.png)

# plotting the distribution of fitness effects (DFEs) with embedded histograms

``` r
hist_df = read.csv("hist_data.csv", header = T)

hist_df = hist_df %>%
  pivot_longer(!Rank, names_to = "Env", values_to = "Mean_fitness")

rank_df1 = filter(hist_df, Env == "InSPI2")

hi = ggplot(rank_df1, aes(Mean_fitness)) +
  geom_histogram(color = "springgreen4", fill = "springgreen4", position="identity", alpha=0.5, binwidth = 0.01) +
  geom_density(alpha=0.2, fill="springgreen4", color = "springgreen4") +
  xlab("Mean Relative fitness (w)") +
  ylab("Frequency") +
  theme(text = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0.0, 1.4), breaks = seq(0.0, 1.4, by = 0.4))

rank_df2 = filter(hist_df, Env == "InSPI2_hypoxic")

ha = ggplot(rank_df2, aes(Mean_fitness)) +
  geom_histogram(color = "springgreen4", fill = "springgreen4", position="identity", alpha=0.5, binwidth = 0.01) +
  geom_density(alpha=0.2, fill="springgreen4", color = "springgreen4") +
  xlab("Mean Relative fitness (w)") +
  ylab("Frequency") +
  theme(text = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0.0, 1.4), breaks = seq(0.0, 1.4, by = 0.4))

rank_df3 = filter(hist_df, Env == "InSPI2_cip")

hc = ggplot(rank_df3, aes(Mean_fitness)) +
  geom_histogram(color = "springgreen4", fill = "springgreen4", position="identity", alpha=0.5, binwidth = 0.01) +
  geom_density(alpha=0.2, fill="springgreen4", color = "springgreen4") +
  xlab("Mean Relative fitness (w)") +
  ylab("Frequency") +
  theme(text = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0.0, 1.4), breaks = seq(0.0, 1.4, by = 0.4))

rank_df4 = filter(hist_df, Env == "InSPI2_lowmg")

hl = ggplot(rank_df4, aes(Mean_fitness)) +
  geom_histogram(color = "springgreen4", fill = "springgreen4", position="identity", alpha=0.5, binwidth = 0.01) +
  geom_density(alpha=0.2, fill="springgreen4", color = "springgreen4") +
  xlab("Mean Relative fitness (w)") +
  ylab("Frequency") +
  theme(text = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0.0, 1.4), breaks = seq(0.0, 1.4, by = 0.4))

df_sort = read.csv("DFE_plot.csv", header = T)

df_sort1 = filter(df_sort, Env == "inspi2")

di = ggplot(df_sort1, aes(Number, mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  geom_hline(yintercept = 1, lty = 2) +
  xlab("Rank order of fitness of transferred genes") +
  ylab("Mean Relative fitness (w)") +
  theme(text = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle("InSPI2") +
  scale_y_continuous(limits = c(0.0, 1.4), breaks = seq(0.0, 1.4, by = 0.2))

df_sort2 = filter(df_sort, Env == "inspi2hypoxic")

da = ggplot(df_sort2, aes(Number, mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  geom_hline(yintercept = 1, lty = 2) +
  xlab("Rank order of fitness of transferred genes") +
  ylab("Mean Relative fitness (w)") +
  theme(text = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle("InSPI2 Hypoxic") +
  scale_y_continuous(limits = c(0.0, 1.4), breaks = seq(0.0, 1.4, by = 0.2))

df_sort3 = filter(df_sort, Env == "inspi2cip")

dc = ggplot(df_sort3, aes(Number, mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  geom_hline(yintercept = 1, lty = 2) +
  xlab("Rank order of fitness of transferred genes") +
  ylab("Mean Relative fitness (w)") +
  theme(text = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle("InSPI2 Ciprofloxacin") +
  scale_y_continuous(limits = c(0.0, 1.4), breaks = seq(0.0, 1.4, by = 0.2))

df_sort4 = filter(df_sort, Env == "inspi2lowmg")

dl = ggplot(df_sort4, aes(Number, mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  geom_hline(yintercept = 1, lty = 2) +
  xlab("Rank order of fitness of transferred genes") +
  ylab("Mean Relative fitness (w)") +
  theme(text = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle("InSPI2 Low Mg") +
  scale_y_continuous(limits = c(0.0, 1.4), breaks = seq(0.0, 1.4, by = 0.2))

ei = di + annotation_custom(ggplotGrob(hi), xmin = 20, xmax = 42, 
                       ymin = 0.1, ymax = 0.75)
```

    ## Warning: Removed 2 rows containing missing values (geom_bar).

``` r
ea = da + annotation_custom(ggplotGrob(ha), xmin = 20, xmax = 42, 
                            ymin = 0.1, ymax = 0.75)
```

    ## Warning: Removed 2 rows containing missing values (geom_bar).

``` r
ec = dc + annotation_custom(ggplotGrob(hc), xmin = 20, xmax = 42, 
                            ymin = 0.1, ymax = 0.75)
```

    ## Warning: Removed 2 rows containing missing values (geom_bar).

``` r
el = dl + annotation_custom(ggplotGrob(hl), xmin = 20, xmax = 42, 
                            ymin = 0.1, ymax = 0.75)
```

    ## Warning: Removed 2 rows containing missing values (geom_bar).

``` r
g = grid.arrange(ei,ea,ec,el, heights = c(7,7),ncol = 2, nrow = 2)
```

![](HGT_barriers_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
ggsave("dfe_embedded_histogram.png",g, height = 6, width = 9)
```

# effect of environment on central tendency of DFEs (Wilcoxon Signed Rank test)

``` r
df1 = read.csv("fitness_deepseq.csv", header = T, stringsAsFactors = T)

df1$Env <- gsub("_", " ", df1$Env)
df1$Env <- gsub("InSPI2", "inspi2", df1$Env)


comp1 = c("inspi2", "inspi2 hypoxic")
comp2 = c("inspi2", "inspi2 cip")
comp3 = c("inspi2", "inspi2 lowmg")
comp4 = c("inspi2 hypoxic", "inspi2 cip")
comp5 = c("inspi2 hypoxic", "inspi2 lowmg")
comp6 = c("inspi2 cip", "inspi2 lowmg")
comparisons = list(comp1, comp2, comp3, comp4, comp5, comp6)

#extract cols needed for the test

df1 = df1[, 1:3]


df_dpsq1 = df1 %>% dplyr::group_by(gene,Env) %>% 
  dplyr::summarise_all(purrr::discard, is.na)

# perform Wilcoxon Rank Sum test in a for loop (use mean fitness values)

wrs_pval = rep(0,6)
wrs_Zstat = rep(0,6)

for(i in 1:length(comparisons)){
  current_comp = comparisons[[i]]
  wrs_df = filter(df_dpsq1, Env %in% current_comp)
  wrs = pivot_wider(wrs_df, names_from = Env, values_from = rel_fitness) %>%
    unnest(cols = everything())
  wrs.test = wilcox.test(unlist(wrs[2]), unlist(wrs[3]),
                         alternative = "two.sided", paired = TRUE)
  wrs_pval[i] = wrs.test$p.value
  wrs_Zstat[i] = qnorm(wrs_pval[i]/2) #calculate Z statistic
}
```

``` r
wrs_padj = p.adjust(wrs_pval, method = "fdr", n = 6)
w = data.frame(wrs_padj, row.names = comparisons)
w
```

    ##                                         wrs_padj
    ## c("inspi2", "inspi2 hypoxic")       4.966615e-02
    ## c("inspi2", "inspi2 cip")           1.416750e-01
    ## c("inspi2", "inspi2 lowmg")         8.891651e-05
    ## c("inspi2 hypoxic", "inspi2 cip")   3.675914e-02
    ## c("inspi2 hypoxic", "inspi2 lowmg") 2.883952e-01
    ## c("inspi2 cip", "inspi2 lowmg")     7.670096e-05

# effect of environment on shape of DFEs (Kolmogorov-Smirnov (K-S) test)

``` r
comp1 = c("inspi2", "inspi2 hypoxic")
comp2 = c("inspi2", "inspi2 cip")
comp3 = c("inspi2", "inspi2 lowmg")
comp4 = c("inspi2 hypoxic", "inspi2 cip")
comp5 = c("inspi2 hypoxic", "inspi2 lowmg")
comp6 = c("inspi2 cip", "inspi2 lowmg")
comparisons = list(comp1, comp2, comp3, comp4, comp5, comp6)

# perform KS test in a for loop (use mean fitness values)

ks_pval = rep(0,6)
ks_Dstat = rep(0,6)

for(i in 1:length(comparisons)){
  current_comp = comparisons[[i]]
  ks_df = filter(df_dpsq1, Env %in% current_comp)
  ks = pivot_wider(ks_df, names_from = Env, values_from = rel_fitness) %>%
    unnest(cols = everything())
  k = ks.test(unlist(ks[2]), unlist(ks[3]),
                         alternative = "two.sided")
  ks_pval[i] = k$p.value
  ks_Dstat[i] = k$statistic
}
```

``` r
ks_padj = p.adjust(ks_pval, method = "fdr", n = 6)
k = data.frame(ks_padj, row.names = comparisons)
k
```

    ##                                        ks_padj
    ## c("inspi2", "inspi2 hypoxic")       0.29477456
    ## c("inspi2", "inspi2 cip")           0.23981804
    ## c("inspi2", "inspi2 lowmg")         0.36359785
    ## c("inspi2 hypoxic", "inspi2 cip")   0.21074191
    ## c("inspi2 hypoxic", "inspi2 lowmg") 0.62490700
    ## c("inspi2 cip", "inspi2 lowmg")     0.00640948

# G X E interactions

``` r
ge = read.csv("deepseq_ge_data.csv", header = T, stringsAsFactors = T, check.names = F)

df = ge%>%
  pivot_longer(cols = starts_with("inspi2"),names_to = "Env", values_to = "W")

head(df)
```

    ## # A tibble: 6 x 4
    ##   gene  Functional_category Env             W
    ##   <fct> <fct>               <chr>       <dbl>
    ## 1 b0119 N                   inspi2      0.972
    ## 2 b0119 N                   inspi2      0.978
    ## 3 b0119 N                   inspi2      0.965
    ## 4 b0119 N                   inspi2      0.960
    ## 5 b0119 N                   inspi2lowmg 0.886
    ## 6 b0119 N                   inspi2lowmg 0.906

``` r
mod2 <- lm(W ~ Env * gene, data=df)

anova(mod2)
```

    ## Analysis of Variance Table
    ## 
    ## Response: W
    ##            Df Sum Sq Mean Sq F value    Pr(>F)    
    ## Env         3  0.804 0.26811  78.755 < 2.2e-16 ***
    ## gene       42 33.275 0.79227 232.720 < 2.2e-16 ***
    ## Env:gene  126  5.715 0.04536  13.324 < 2.2e-16 ***
    ## Residuals 516  1.757 0.00340                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
posthoc2 <- lsmeans(mod2, pairwise ~ Env|gene, adjust="Tukey")

contr2 <- data.frame(posthoc2$contrasts) %>% mutate(sign = ifelse(p.value<0.05, "*", ""))

kable(contr2)
```

<table>
<colgroup>
<col width="30%" />
<col width="7%" />
<col width="12%" />
<col width="11%" />
<col width="5%" />
<col width="13%" />
<col width="11%" />
<col width="6%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">contrast</th>
<th align="left">gene</th>
<th align="right">estimate</th>
<th align="right">SE</th>
<th align="right">df</th>
<th align="right">t.ratio</th>
<th align="right">p.value</th>
<th align="left">sign</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0119</td>
<td align="right">0.1676326</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">4.0630749</td>
<td align="right">0.0003256</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0119</td>
<td align="right">-0.0243456</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.5900879</td>
<td align="right">0.9350665</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0119</td>
<td align="right">0.0563496</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.3658014</td>
<td align="right">0.5214543</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0119</td>
<td align="right">-0.1919782</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-4.6531628</td>
<td align="right">0.0000246</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0119</td>
<td align="right">-0.1112829</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.6972735</td>
<td align="right">0.0362110</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0119</td>
<td align="right">0.0806952</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.9558893</td>
<td align="right">0.2062848</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0127</td>
<td align="right">0.1323797</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">3.2086160</td>
<td align="right">0.0077106</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0127</td>
<td align="right">0.1120085</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">2.7148586</td>
<td align="right">0.0344917</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0127</td>
<td align="right">-0.1496725</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.6277592</td>
<td align="right">0.0017819</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0127</td>
<td align="right">-0.0203712</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.4937574</td>
<td align="right">0.9604685</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0127</td>
<td align="right">-0.2820522</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-6.8363753</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0127</td>
<td align="right">-0.2616810</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-6.3426179</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0179</td>
<td align="right">0.0503458</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.2202803</td>
<td align="right">0.6142221</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0179</td>
<td align="right">-0.2387198</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-5.7860865</td>
<td align="right">0.0000001</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0179</td>
<td align="right">-0.0359897</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.8723179</td>
<td align="right">0.8191867</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0179</td>
<td align="right">-0.2890656</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-7.0063667</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0179</td>
<td align="right">-0.0863355</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.0925982</td>
<td align="right">0.1568500</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0179</td>
<td align="right">0.2027301</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">4.9137685</td>
<td align="right">0.0000071</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0215</td>
<td align="right">0.1794261</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">4.3489250</td>
<td align="right">0.0000969</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0215</td>
<td align="right">0.0518332</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.2563327</td>
<td align="right">0.5912263</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0215</td>
<td align="right">-0.0413646</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.0025935</td>
<td align="right">0.7479532</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0215</td>
<td align="right">-0.1275928</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.0925923</td>
<td align="right">0.0112147</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0215</td>
<td align="right">-0.2207906</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-5.3515186</td>
<td align="right">0.0000008</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0215</td>
<td align="right">-0.0931978</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.2589262</td>
<td align="right">0.1091670</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0423</td>
<td align="right">-0.0221804</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.5376081</td>
<td align="right">0.9498172</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0423</td>
<td align="right">-0.0066827</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.1619754</td>
<td align="right">0.9984865</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0423</td>
<td align="right">-0.0757970</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.8371670</td>
<td align="right">0.2570831</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0423</td>
<td align="right">0.0154977</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.3756327</td>
<td align="right">0.9819180</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0423</td>
<td align="right">-0.0536166</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.2995589</td>
<td align="right">0.5636084</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0423</td>
<td align="right">-0.0691143</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.6751916</td>
<td align="right">0.3378835</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0440</td>
<td align="right">-0.3872076</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-9.3851291</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0440</td>
<td align="right">-0.0228736</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.5544091</td>
<td align="right">0.9453340</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0440</td>
<td align="right">-0.1103517</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.6747025</td>
<td align="right">0.0385251</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0440</td>
<td align="right">0.3643340</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">8.8307199</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0440</td>
<td align="right">0.2768559</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">6.7104265</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0440</td>
<td align="right">-0.0874781</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.1202934</td>
<td align="right">0.1479902</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0607</td>
<td align="right">0.0975210</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">2.3637122</td>
<td align="right">0.0854930</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0607</td>
<td align="right">-0.0111594</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.2704807</td>
<td align="right">0.9930754</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0607</td>
<td align="right">0.0053974</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.1308231</td>
<td align="right">0.9991998</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0607</td>
<td align="right">-0.1086804</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.6341928</td>
<td align="right">0.0429966</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0607</td>
<td align="right">-0.0921236</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.2328891</td>
<td align="right">0.1157805</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0607</td>
<td align="right">0.0165568</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.4013037</td>
<td align="right">0.9781148</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0623</td>
<td align="right">-0.1129096</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.7367002</td>
<td align="right">0.0324546</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0623</td>
<td align="right">-0.0638047</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.5464962</td>
<td align="right">0.4104826</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0623</td>
<td align="right">-0.0591548</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.4337933</td>
<td align="right">0.4787872</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0623</td>
<td align="right">0.0491049</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.1902040</td>
<td align="right">0.6333221</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0623</td>
<td align="right">0.0537548</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.3029070</td>
<td align="right">0.5614705</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0623</td>
<td align="right">0.0046499</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.1127030</td>
<td align="right">0.9994875</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0642</td>
<td align="right">-0.0391414</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.9487077</td>
<td align="right">0.7784986</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0642</td>
<td align="right">-0.0030344</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.0735479</td>
<td align="right">0.9998572</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0642</td>
<td align="right">-0.0478569</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.1599534</td>
<td align="right">0.6524087</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0642</td>
<td align="right">0.0361070</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.8751598</td>
<td align="right">0.8177340</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0642</td>
<td align="right">-0.0087155</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.2112457</td>
<td align="right">0.9966657</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0642</td>
<td align="right">-0.0448224</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.0864055</td>
<td align="right">0.6980187</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0695</td>
<td align="right">-0.0279164</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.6766367</td>
<td align="right">0.9059357</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0695</td>
<td align="right">-0.0466646</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.1310567</td>
<td align="right">0.6704838</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0695</td>
<td align="right">-0.1957725</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-4.7451289</td>
<td align="right">0.0000160</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0695</td>
<td align="right">-0.0187483</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.4544200</td>
<td align="right">0.9687557</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0695</td>
<td align="right">-0.1678561</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-4.0684922</td>
<td align="right">0.0003184</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0695</td>
<td align="right">-0.1491078</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.6140722</td>
<td align="right">0.0018743</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0780</td>
<td align="right">0.1551352</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">3.7601626</td>
<td align="right">0.0010828</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0780</td>
<td align="right">-0.0395080</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.9575943</td>
<td align="right">0.7735583</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0780</td>
<td align="right">0.0141288</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.3424545</td>
<td align="right">0.9861762</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0780</td>
<td align="right">-0.1946432</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-4.7177569</td>
<td align="right">0.0000182</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0780</td>
<td align="right">-0.1410063</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.4177081</td>
<td align="right">0.0037939</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0780</td>
<td align="right">0.0536368</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.3000488</td>
<td align="right">0.5632956</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0785</td>
<td align="right">0.1378547</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">3.3413193</td>
<td align="right">0.0049408</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0785</td>
<td align="right">0.0130861</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.3171797</td>
<td align="right">0.9889480</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0785</td>
<td align="right">-0.0108899</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.2639484</td>
<td align="right">0.9935566</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0785</td>
<td align="right">-0.1247686</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.0241396</td>
<td align="right">0.0138991</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0785</td>
<td align="right">-0.1487446</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.6052677</td>
<td align="right">0.0019361</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0785</td>
<td align="right">-0.0239759</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.5811281</td>
<td align="right">0.9377406</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0812</td>
<td align="right">0.1272226</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">3.0836183</td>
<td align="right">0.0115379</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0812</td>
<td align="right">-0.0024011</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.0581986</td>
<td align="right">0.9999292</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0812</td>
<td align="right">-0.0131767</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.3193771</td>
<td align="right">0.9887226</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0812</td>
<td align="right">-0.1296237</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.1418169</td>
<td align="right">0.0095826</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0812</td>
<td align="right">-0.1403993</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.4029954</td>
<td align="right">0.0039937</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0812</td>
<td align="right">-0.0107756</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.2611785</td>
<td align="right">0.9937540</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0880</td>
<td align="right">0.7412731</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">17.9669615</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0880</td>
<td align="right">0.2082099</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">5.0465870</td>
<td align="right">0.0000037</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0880</td>
<td align="right">-0.1114164</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.7005073</td>
<td align="right">0.0358894</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0880</td>
<td align="right">-0.5330632</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-12.9203745</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0880</td>
<td align="right">-0.8526895</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-20.6674688</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0880</td>
<td align="right">-0.3196263</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-7.7470944</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0882</td>
<td align="right">-0.1770832</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-4.2921392</td>
<td align="right">0.0001241</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0882</td>
<td align="right">0.0775333</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.8792504</td>
<td align="right">0.2382383</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0882</td>
<td align="right">-0.3229414</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-7.8274478</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0882</td>
<td align="right">0.2546165</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">6.1713896</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0882</td>
<td align="right">-0.1458582</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.5353086</td>
<td align="right">0.0024982</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0882</td>
<td align="right">-0.4004747</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-9.7066982</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0891</td>
<td align="right">-0.0275809</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.6685059</td>
<td align="right">0.9089251</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0891</td>
<td align="right">-0.0955451</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.3158205</td>
<td align="right">0.0957451</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0891</td>
<td align="right">-0.0005353</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.0129740</td>
<td align="right">0.9999992</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0891</td>
<td align="right">-0.0679642</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.6473147</td>
<td align="right">0.3530327</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0891</td>
<td align="right">0.0270456</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.6555319</td>
<td align="right">0.9135876</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0891</td>
<td align="right">0.0950098</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">2.3028465</td>
<td align="right">0.0986844</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b0948</td>
<td align="right">-0.0224233</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.5434947</td>
<td align="right">0.9482719</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b0948</td>
<td align="right">-0.0030406</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.0736983</td>
<td align="right">0.9998563</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b0948</td>
<td align="right">-0.0486194</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.1784367</td>
<td align="right">0.6407639</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b0948</td>
<td align="right">0.0193827</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.4697964</td>
<td align="right">0.9656568</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b0948</td>
<td align="right">-0.0261962</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.6349419</td>
<td align="right">0.9207134</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b0948</td>
<td align="right">-0.0455788</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.1047383</td>
<td align="right">0.6867786</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b1000</td>
<td align="right">0.0823084</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.9949893</td>
<td align="right">0.1911628</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b1000</td>
<td align="right">-0.0297895</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.7220375</td>
<td align="right">0.8882992</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b1000</td>
<td align="right">0.0272986</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.6616618</td>
<td align="right">0.9114012</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b1000</td>
<td align="right">-0.1120979</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.7170268</td>
<td align="right">0.0342847</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b1000</td>
<td align="right">-0.0550098</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.3333275</td>
<td align="right">0.5420730</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b1000</td>
<td align="right">0.0570881</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.3836992</td>
<td align="right">0.5101465</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b1094</td>
<td align="right">-0.1167230</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.8291297</td>
<td align="right">0.0249434</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b1094</td>
<td align="right">-0.0754991</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.8299457</td>
<td align="right">0.2604084</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b1094</td>
<td align="right">-0.0800014</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.9390732</td>
<td align="right">0.2130335</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b1094</td>
<td align="right">0.0412239</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.9991839</td>
<td align="right">0.7499260</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b1094</td>
<td align="right">0.0367216</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.8900564</td>
<td align="right">0.8100385</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b1094</td>
<td align="right">-0.0045023</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.1091275</td>
<td align="right">0.9995346</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b1290</td>
<td align="right">0.1558929</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">3.7785296</td>
<td align="right">0.0010092</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b1290</td>
<td align="right">-0.0464038</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.1247335</td>
<td align="right">0.6744145</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b1290</td>
<td align="right">0.0162368</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.3935479</td>
<td align="right">0.9793119</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b1290</td>
<td align="right">-0.2022967</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-4.9032631</td>
<td align="right">0.0000075</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b1290</td>
<td align="right">-0.1396561</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.3849817</td>
<td align="right">0.0042514</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b1290</td>
<td align="right">0.0626406</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.5182815</td>
<td align="right">0.4272285</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b1686</td>
<td align="right">0.0436097</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.0570101</td>
<td align="right">0.7158237</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b1686</td>
<td align="right">-0.0551093</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.3357382</td>
<td align="right">0.5405387</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b1686</td>
<td align="right">0.0032719</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.0793047</td>
<td align="right">0.9998210</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b1686</td>
<td align="right">-0.0987190</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.3927483</td>
<td align="right">0.0797200</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b1686</td>
<td align="right">-0.0403377</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.9777054</td>
<td align="right">0.7622327</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b1686</td>
<td align="right">0.0583812</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.4150429</td>
<td align="right">0.4904688</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b1718</td>
<td align="right">0.0569636</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.3806831</td>
<td align="right">0.5120488</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b1718</td>
<td align="right">-0.0777153</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.8836610</td>
<td align="right">0.2363163</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b1718</td>
<td align="right">0.0156536</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.3794111</td>
<td align="right">0.9813866</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b1718</td>
<td align="right">-0.1346789</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.2643441</td>
<td align="right">0.0064099</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b1718</td>
<td align="right">-0.0413100</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.0012720</td>
<td align="right">0.7487184</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b1718</td>
<td align="right">0.0933688</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">2.2630721</td>
<td align="right">0.1081418</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b1763</td>
<td align="right">0.0000000</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.0000000</td>
<td align="right">1.0000000</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b1763</td>
<td align="right">-0.2788015</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-6.7575855</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b1763</td>
<td align="right">0.0000000</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.0000000</td>
<td align="right">1.0000000</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b1763</td>
<td align="right">-0.2788015</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-6.7575855</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b1763</td>
<td align="right">0.0000000</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.0000000</td>
<td align="right">1.0000000</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b1763</td>
<td align="right">0.2788015</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">6.7575855</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b1913</td>
<td align="right">0.1008367</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">2.4440766</td>
<td align="right">0.0702905</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b1913</td>
<td align="right">-0.0853801</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.0694402</td>
<td align="right">0.1645521</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b1913</td>
<td align="right">-0.0061491</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.1490428</td>
<td align="right">0.9988190</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b1913</td>
<td align="right">-0.1862167</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-4.5135168</td>
<td align="right">0.0000467</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b1913</td>
<td align="right">-0.1069858</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.5931193</td>
<td align="right">0.0479734</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b1913</td>
<td align="right">0.0792309</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.9203975</td>
<td align="right">0.2207016</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b2341</td>
<td align="right">-0.0390382</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.9462065</td>
<td align="right">0.7798819</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b2341</td>
<td align="right">-0.0071191</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.1725521</td>
<td align="right">0.9981726</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b2341</td>
<td align="right">-0.0664240</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.6099843</td>
<td align="right">0.3738383</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b2341</td>
<td align="right">0.0319191</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.7736544</td>
<td align="right">0.8663538</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b2341</td>
<td align="right">-0.0273859</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.6637779</td>
<td align="right">0.9106396</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b2341</td>
<td align="right">-0.0593050</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.4374322</td>
<td align="right">0.4765289</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b2530</td>
<td align="right">-0.0052784</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.1279375</td>
<td align="right">0.9992514</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b2530</td>
<td align="right">-0.0625212</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.5153866</td>
<td align="right">0.4289614</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b2530</td>
<td align="right">-0.0821405</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.9909198</td>
<td align="right">0.1926996</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b2530</td>
<td align="right">-0.0572428</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.3874491</td>
<td align="right">0.5077834</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b2530</td>
<td align="right">-0.0768621</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.8629823</td>
<td align="right">0.2454146</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b2530</td>
<td align="right">-0.0196193</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.4755332</td>
<td align="right">0.9644547</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b2576</td>
<td align="right">-0.1523505</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.6926690</td>
<td align="right">0.0013987</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b2576</td>
<td align="right">0.0771554</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.8700914</td>
<td align="right">0.2422617</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b2576</td>
<td align="right">-0.3052587</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-7.3988543</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b2576</td>
<td align="right">0.2295060</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">5.5627603</td>
<td align="right">0.0000003</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b2576</td>
<td align="right">-0.1529082</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.7061853</td>
<td align="right">0.0013293</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b2576</td>
<td align="right">-0.3824141</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-9.2689456</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b2990</td>
<td align="right">0.0266454</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.6458316</td>
<td align="right">0.9169866</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b2990</td>
<td align="right">-0.0490881</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.1897970</td>
<td align="right">0.6335798</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b2990</td>
<td align="right">0.0024225</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.0587174</td>
<td align="right">0.9999273</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b2990</td>
<td align="right">-0.0757336</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.8356286</td>
<td align="right">0.2577893</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b2990</td>
<td align="right">-0.0242229</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.5871143</td>
<td align="right">0.9359611</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b2990</td>
<td align="right">0.0515107</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.2485143</td>
<td align="right">0.5962191</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b3006</td>
<td align="right">0.0392114</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.9504041</td>
<td align="right">0.7775587</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b3006</td>
<td align="right">-0.0342043</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.8290429</td>
<td align="right">0.8406776</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b3006</td>
<td align="right">-0.0100605</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.2438453</td>
<td align="right">0.9949001</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b3006</td>
<td align="right">-0.0734157</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.7794470</td>
<td align="right">0.2844023</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b3006</td>
<td align="right">-0.0492718</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.1942494</td>
<td align="right">0.6307592</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b3006</td>
<td align="right">0.0241438</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.5851976</td>
<td align="right">0.9365340</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b3071</td>
<td align="right">0.0040897</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.0991252</td>
<td align="right">0.9996510</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b3071</td>
<td align="right">-0.0429313</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.0405682</td>
<td align="right">0.7256529</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b3071</td>
<td align="right">-0.0193434</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.4688439</td>
<td align="right">0.9658539</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b3071</td>
<td align="right">-0.0470210</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.1396935</td>
<td align="right">0.6650999</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b3071</td>
<td align="right">-0.0234330</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.5679692</td>
<td align="right">0.9415516</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b3071</td>
<td align="right">0.0235880</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.5717243</td>
<td align="right">0.9404782</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b3164</td>
<td align="right">-0.0537606</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.3030474</td>
<td align="right">0.5613809</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b3164</td>
<td align="right">-0.0799754</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.9384424</td>
<td align="right">0.2132895</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b3164</td>
<td align="right">-0.0457022</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.1077288</td>
<td align="right">0.6849361</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b3164</td>
<td align="right">-0.0262149</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.6353950</td>
<td align="right">0.9205603</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b3164</td>
<td align="right">0.0080584</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.1953185</td>
<td align="right">0.9973580</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b3164</td>
<td align="right">0.0342732</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.8307136</td>
<td align="right">0.8398705</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b3417</td>
<td align="right">0.4373322</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">10.6000477</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b3417</td>
<td align="right">-0.0535246</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.2973292</td>
<td align="right">0.5650324</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b3417</td>
<td align="right">0.1645875</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">3.9892676</td>
<td align="right">0.0004398</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b3417</td>
<td align="right">-0.4908568</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-11.8973769</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b3417</td>
<td align="right">-0.2727447</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-6.6107801</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b3417</td>
<td align="right">0.2181121</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">5.2865968</td>
<td align="right">0.0000011</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b3560</td>
<td align="right">0.0269104</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.6522544</td>
<td align="right">0.9147443</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b3560</td>
<td align="right">0.0099484</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.2411284</td>
<td align="right">0.9950662</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b3560</td>
<td align="right">-0.0794341</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.9253221</td>
<td align="right">0.2186619</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b3560</td>
<td align="right">-0.0169621</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.4111260</td>
<td align="right">0.9765379</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b3560</td>
<td align="right">-0.1063445</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.5775765</td>
<td align="right">0.0499797</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b3560</td>
<td align="right">-0.0893825</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.1664505</td>
<td align="right">0.1340593</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b3590</td>
<td align="right">-0.1095292</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.6547653</td>
<td align="right">0.0406734</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b3590</td>
<td align="right">0.0323026</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.7829501</td>
<td align="right">0.8621946</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b3590</td>
<td align="right">-0.1829942</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-4.4354106</td>
<td align="right">0.0000662</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b3590</td>
<td align="right">0.1418318</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">3.4377154</td>
<td align="right">0.0035370</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b3590</td>
<td align="right">-0.0734651</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.7806453</td>
<td align="right">0.2838181</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b3590</td>
<td align="right">-0.2152969</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-5.2183607</td>
<td align="right">0.0000016</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b3602</td>
<td align="right">0.0494614</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.1988440</td>
<td align="right">0.6278459</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b3602</td>
<td align="right">0.0123959</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.3004524</td>
<td align="right">0.9905696</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b3602</td>
<td align="right">-0.0593239</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.4378918</td>
<td align="right">0.4762439</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b3602</td>
<td align="right">-0.0370655</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.8983917</td>
<td align="right">0.8056747</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b3602</td>
<td align="right">-0.1087853</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.6367358</td>
<td align="right">0.0427034</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b3602</td>
<td align="right">-0.0717198</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.7383441</td>
<td align="right">0.3048686</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b3686</td>
<td align="right">0.0817559</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.9815983</td>
<td align="right">0.1962522</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b3686</td>
<td align="right">-0.0151520</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.3672531</td>
<td align="right">0.9830620</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b3686</td>
<td align="right">-0.0040805</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.0989030</td>
<td align="right">0.9996533</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b3686</td>
<td align="right">-0.0969079</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.3488514</td>
<td align="right">0.0885753</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b3686</td>
<td align="right">-0.0858364</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.0805013</td>
<td align="right">0.1608397</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b3686</td>
<td align="right">0.0110715</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.2683501</td>
<td align="right">0.9932348</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b3725</td>
<td align="right">-0.0313370</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.7595450</td>
<td align="right">0.8725477</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b3725</td>
<td align="right">0.0218425</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.5294180</td>
<td align="right">0.9519214</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b3725</td>
<td align="right">-0.1111272</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.6934985</td>
<td align="right">0.0365895</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b3725</td>
<td align="right">0.0531795</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.2889630</td>
<td align="right">0.5703769</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b3725</td>
<td align="right">-0.0797902</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.9339535</td>
<td align="right">0.2151175</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b3725</td>
<td align="right">-0.1329697</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.2229165</td>
<td align="right">0.0073558</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b4000</td>
<td align="right">-0.3492345</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-8.4647378</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b4000</td>
<td align="right">-0.0342983</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.8313208</td>
<td align="right">0.8395766</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b4000</td>
<td align="right">-0.1494884</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.6232964</td>
<td align="right">0.0018116</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b4000</td>
<td align="right">0.3149362</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">7.6334169</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b4000</td>
<td align="right">0.1997461</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">4.8414414</td>
<td align="right">0.0000101</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b4000</td>
<td align="right">-0.1151901</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.7919756</td>
<td align="right">0.0277576</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b4043</td>
<td align="right">0.4403134</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">10.6723076</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b4043</td>
<td align="right">-0.1193461</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.8927091</td>
<td align="right">0.0207036</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b4043</td>
<td align="right">-0.0269305</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.6527411</td>
<td align="right">0.9145731</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b4043</td>
<td align="right">-0.5596596</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-13.5650167</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b4043</td>
<td align="right">-0.4672439</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-11.3250487</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b4043</td>
<td align="right">0.0924156</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">2.2399680</td>
<td align="right">0.1139524</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b4172</td>
<td align="right">-0.1567681</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-3.7997409</td>
<td align="right">0.0009300</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b4172</td>
<td align="right">0.0696573</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.6883528</td>
<td align="right">0.3308512</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b4172</td>
<td align="right">-0.1993868</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-4.8327329</td>
<td align="right">0.0000106</td>
<td align="left">*</td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b4172</td>
<td align="right">0.2264254</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">5.4880937</td>
<td align="right">0.0000004</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b4172</td>
<td align="right">-0.0426187</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.0329920</td>
<td align="right">0.7301481</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b4172</td>
<td align="right">-0.2690441</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-6.5210857</td>
<td align="right">0.0000000</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b4203</td>
<td align="right">0.0678270</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.6439891</td>
<td align="right">0.3548625</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b4203</td>
<td align="right">-0.0549467</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.3317967</td>
<td align="right">0.5430476</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b4203</td>
<td align="right">0.0186766</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.4526838</td>
<td align="right">0.9690944</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b4203</td>
<td align="right">-0.1227737</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.9757858</td>
<td align="right">0.0161272</td>
<td align="left">*</td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b4203</td>
<td align="right">-0.0491504</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.1913053</td>
<td align="right">0.6326246</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b4203</td>
<td align="right">0.0736233</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.7844805</td>
<td align="right">0.2819532</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b4243</td>
<td align="right">0.0339361</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.8225425</td>
<td align="right">0.8438001</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b4243</td>
<td align="right">-0.0427536</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.0362610</td>
<td align="right">0.7282112</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b4243</td>
<td align="right">-0.0030235</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.0732825</td>
<td align="right">0.9998587</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b4243</td>
<td align="right">-0.0766897</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.8588035</td>
<td align="right">0.2472801</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b4243</td>
<td align="right">-0.0369596</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.8958249</td>
<td align="right">0.8070229</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b4243</td>
<td align="right">0.0397301</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.9629785</td>
<td align="right">0.7705456</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2cip</td>
<td align="left">b4373</td>
<td align="right">0.0396436</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">0.9608811</td>
<td align="right">0.7717209</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2 - inspi2hypoxic</td>
<td align="left">b4373</td>
<td align="right">-0.0498139</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.2073872</td>
<td align="right">0.6224221</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2 - inspi2lowmg</td>
<td align="left">b4373</td>
<td align="right">-0.0017911</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-0.0434130</td>
<td align="right">0.9999706</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2cip - inspi2hypoxic</td>
<td align="left">b4373</td>
<td align="right">-0.0894575</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-2.1682683</td>
<td align="right">0.1335316</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">inspi2cip - inspi2lowmg</td>
<td align="left">b4373</td>
<td align="right">-0.0414347</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">-1.0042941</td>
<td align="right">0.7469674</td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left">inspi2hypoxic - inspi2lowmg</td>
<td align="left">b4373</td>
<td align="right">0.0480227</td>
<td align="right">0.0412576</td>
<td align="right">516</td>
<td align="right">1.1639742</td>
<td align="right">0.6498806</td>
<td align="left"></td>
</tr>
</tbody>
</table>

``` r
contr2 <- contr2 %>% filter(sign == "*") %>% 
  arrange(contrast) %>% kable()

fix(contr2) # save as txt file
```

# G X E plot

``` r
rank_df1 = read.csv("fitness_ranked_by_inspi2_dpsq.csv", header = T)

plot_rank = rank_df1 %>%
  pivot_longer(cols = starts_with("In"), names_to = "Env", 
                        values_to = "mean_fitness")

plot_rank$Env = factor(plot_rank$Env, levels = c("InSPI2", "InSPI2_hypoxic", 
                                                   "InSPI2_cip", "InSPI2_lowmg"), labels = c("InSPI2", "InSPI2 Hypoxic","InSPI2 Ciprofloxacin", "InSPI2 Low Mg"))

lp = ggplot(plot_rank, aes(x = Rank, y = mean_fitness, color=Env)) +
  scale_color_paletteer_d("jcolors::pal5") +
  geom_line(size = 0.6) +
  xlab("Rank order of fitness of transferred genes\n(ranked by w in InSPI2)") +
  ylab("Mean Relative Fitness (w)") +
  geom_hline(yintercept = 1, lty = 2) +
  theme(legend.title = element_blank()) +
  scale_y_continuous(limits = c(0.0, 1.4), breaks = seq(0.0, 1.4, by = 0.2)) +
  theme(text = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white",
                                        size = 1, linetype = "solid"))

ggsave("fitness_ranked_by_inspi2.png", lp, height = 6, width = 8)
lp
```

![](HGT_barriers_files/figure-markdown_github/unnamed-chunk-18-1.png)
