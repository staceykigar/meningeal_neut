# Author: Stacey L. Kigar
# 20230923

# set-up ------------------------------------------------------------------

# load packages
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggplot2)
library(rcompanion)
library(report)
library(showtext)

#import data
setwd("/data/")

all <- read_csv("C57_data.csv")

# make batch column
all$batch <- paste0(all$study, "_", all$cohort)


# Assess independence of main predictors ----------------------------------
# i know i want to control for batch in my linear modeling, so compare
# group distribution within batch
table(all$group, all$batch)

#       Exo14NP_1 Exo14NP_2 NP10_1 NP10_2 NP10_3 NP10_4 NP13_1 NP3_1 NP3_3 NP3_4
# CSD         3         3      2      2      3      2      4     4     2     1
# HC          3         2      2      2      2      2      4     2     2     2

# all batches are balanced but numbers are small

# Wrangle data ------------------------------------------------------------

# make things factors:
all$group <- factor(all$group, 
                    levels = c("HC", "CSD"))
all$batch <- factor(all$batch)
all$study <- factor(all$study)
all$cohort <- factor(all$cohort)

# rename columns 
all  %<>%  rename(Group = group,
                  Blood_mono = bMO.iv,
                  Blood_NPs = bNP.iv,
                  Blood_T = bT.iv,
                  Blood_B = bB.iv)

# make 'unique ID column':
all$SampleID <- paste0(all$study, "_", all$sample)


# binarize USM
all %<>% mutate(USM_binary = if_else(USM > 0, 1, 0))


# Transform  --------------------------------------------------------------

cols <- all %>% select(mT.nv, mB.iv, Blood_T, Blood_mono, 
                       mNP.iv, mT.iv,  Blood_NPs) %>% colnames()

# store as new data frame because original column data will be overwritten
trans <- all %>% mutate_at(cols, sqrt) 


# Linear mixed modeling: OF  -----------------------------------------------

# set viewing window for looking at residual graphs:
par(mfrow=c(2,2))

# simplify dataframe by dropping levels for NA OF data:
final <- trans %>% drop_na(OF_cross) %>% droplevels()


#' `blood - neutrophils`
#' run model:
model = lm(OF_cross ~ Blood_NPs + batch, 
           data = final, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
report(model)


#' `meninges - iv+ neutrophils`
#' run model:
model = lm(OF_cross ~ mNP.iv + batch, 
           data = final, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)

#' `meninges - iv- neutrophils`
#' run model:
model = lm(OF_cross ~ mNP.nv + batch, 
           data = final, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


# Mixed logistic regression - USM   -----------------------------------------

# simplify dataframe by dropping levels for NA OF data:
final <- trans %>% drop_na(USM_binary) %>% droplevels()


#' `blood - neutrophils`
#' run model:
model = glm(USM_binary ~ Blood_NPs + batch, data = final, 
            family = binomial, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


#' `meninges - iv+ neutrophils`
#' run model:
model = glm(USM_binary ~ mNP.iv + batch, data = final, 
            family = binomial, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


#' `meninges - iv- neutrophils`
#' run model:
model = glm(USM_binary ~ mNP.nv + batch, data = final, 
            family = binomial, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model)
# get plain english summary:
report(model)


# linear model for blood vs meningeal neuts -------------------------------

# simplify dataframe by dropping levels for NA OF data:
final <- trans %>% drop_na(mNP.nv) %>% droplevels()

# generate model:
model <- final %>%  
  lm(formula = mNP.nv ~ Blood_NPs * Group + batch, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


# Close split plotting window
dev.off()

# interaction plot --------------------------------------------------------
library(interactions)

#' `meningeal neut (nv) vs blood neuts`

# generate model: 
mNPs_vbNPs_study <- final %>%  
  lm(formula = mNP.nv ~ Blood_NPs * Group + batch)


# plot and save
interactions::interact_plot(mNPs_vbNPs_study, pred = Blood_NPs, 
                            modx = Group, plot.points = T,
                            line.thickness = 2,
                            point.size = 12, partial.residuals = T) +
  theme_classic() + 
  xlab("\n%blood neutrophils") + 
  ylab("%(iv-) meningeal neutrophils\n") +
  theme(text = element_text(size = 36, 
                            family = "Arial", colour = "black")) +
  scale_color_manual(values = c('#636463', '#FD4441'))
# ggsave("f1_bNPs_vmNPsnv_lm.pdf")


# get parameters of line:
interactions::sim_slopes(mNPs_vbNPs_study, pred = Blood_NPs, 
           modx = Group)

# Slope of Blood_NPs when Group = CSD: 
#   
#   Est.   S.E.   t val.      p
# ------- ------ -------- ------
#   -1.56   1.13    -1.38   0.18
# 
# Slope of Blood_NPs when Group = HC: 
#   
#   Est.   S.E.   t val.      p
# ------- ------ -------- ------
#   -0.78   0.86    -0.91   0.37

# dot & whiskers graph set up --------------------------------------------------------------
library(stringr)
library(dotwhisker)


# select output directory for graphs
setwd("/results/")

# import Arial font
font_add("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_auto()

# make R stop doing scientific notation: 
options(scipen=999)


# Dot & whiskers plot - OF ----------------------------------------------------

final <- trans %>% drop_na(OF_cross) %>% droplevels()
colnames(final)

# save variables to be graphed:
ordered_vars <- c("mNP.nv",
                  "mNP.iv",
                  "Blood_NPs")

#' create tidy df for linear model - `OF cross`
m <- list()
m[[1]] <- lm(OF_cross ~ mNP.nv + batch, data = final, 
             na.action = na.exclude)

of_df <- m[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1")

# loop for remaining variables
for (i in 2:3) {
  m[[i]] <- lm(str_replace("OF_cross ~ meow + batch", "meow", 
                           ordered_vars[i]), 
               data = final, na.action = na.exclude)
  of_df <- rbind(of_df, 
                 m[[i]] %>% 
                   broom::tidy() %>% 
                   by_2sd(final) %>% 
                   mutate(model = paste("Model", i)))
}

# don't care about the cohort or intercept terms, filter:
of_df_mod <- of_df %>% filter(!grepl('batch*', term)) %>% 
  filter(!grepl('Intercept*', term))

#' make the graph - `OF cross`
dwplot(of_df_mod, 
       vline = geom_vline(
         xintercept = 0,
         colour = "black",
         linetype = 5,
         linewidth = 2
       ),
       vars_order = ordered_vars,
       dot_args = list(size = 5),
       whisker_args = list(size = 2)) %>% 
  relabel_predictors(
    c(mNP.nv = "%meningeal \nneutrophils (iv-)",
      mNP.iv = "%meningeal \nneutrophils (iv+)",
      Blood_NPs = "%blood \nneutrophils" )
  ) +
  theme_pubr() + xlab("Coefficient estimate (OF)") +
  theme(text = element_text(size = 28, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441", "#FD4441", "#FD4441")) +
  annotate(geom = "text", label = "*", x = -9.84, y = 1.25, size = 12) +
  xlim(-30,5)
# ggsave("f3_WT_OF_coefficients.pdf")

# multiple comparisons correction:
p.adjust(of_df_mod$p.value, method="BH")


# add p.adj values to tibble:
of_df_mod$p.adj <- p.adjust(of_df_mod$p.value, method="BH")

# save tibble for future table creation:
# of_df_mod %>% write_csv("WT_OF_neut_mlm_stats.csv")


# Dot & whiskers plot - USM ------------------------------------------------
final <- trans %>% drop_na(USM_binary) %>% droplevels()
colnames(final)

#' Create tidy df of coefficient estimates - `USM`
m2 <- list()

m2[[1]] <- glm(USM_binary ~ mNP.nv + batch, data = final, 
               family = binomial, na.action = na.exclude)

usm_df <- m2[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1") 

# loop for remaining variables
for (i in 2:3) {
  m2[[i]] <- glm(str_replace("USM_binary ~ meow + batch", "meow", 
                             ordered_vars[i]), 
                 data = final, family = binomial, na.action = na.exclude)
  usm_df <- rbind(usm_df, 
                  m2[[i]] %>% 
                    broom::tidy() %>% 
                    by_2sd(final) %>% 
                    mutate(model = paste("Model", i))) 
  
}

# don't care about the cohort or intercept terms, filter:
usm_df_mod <- usm_df %>% filter(!grepl('batch*', term)) %>% 
  filter(!grepl('Intercept*', term))

#' make the graph - `USM`
dwplot(usm_df_mod, 
       vline = geom_vline(
         xintercept = 0,
         colour = "black",
         linetype = 5,
         linewidth = 2
       ),
       vars_order = ordered_vars,
       dot_args = list(size = 5),
       whisker_args = list(size = 2)) %>% 
  relabel_predictors(
    c(mNP.nv = "%meningeal \nneutrophils (iv-)",
      mNP.iv = "%meningeal \nneutrophils (iv+)",
      Blood_NPs = "%blood \nneutrophils"
    )
  ) +
  theme_pubr() + xlab("Coefficient estimate (USM)") +
  theme(text = element_text(size = 28, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441", "#FD4441", "#FD4441")) +
  annotate(geom = "text", label = "*", x = -6.43, y = 1.25, size = 12) +
  annotate(geom = "text", label = "*", x = -5.49, y = 2.1, size = 12) +
  annotate(geom = "text", label = "*", x = -5.72, y = 2.97, size = 12) +
  xlim(-14.5,1)
# ggsave("f3_WT_USM_coefficients.pdf")

# multiple comparisons correction:
p.adjust(usm_df_mod$p.value, method="BH")


# add p.adj values to tibble:
usm_df_mod$p.adj <- p.adjust(usm_df_mod$p.value, method="BH")

# save tibble for future table creation:
# usm_df_mod %>% write_csv("WT_USM_neut_glm_stats.csv")



