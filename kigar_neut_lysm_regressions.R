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

# import data for LysM mice:
all <- read_csv("LysM_data.csv")

# Assess independence of main predictors ----------------------------------
# want to control for cohort in my linear modeling, so compare
# group distribution within Cohort
table(all$Group, all$Cohort)
#     1 2 4 5 6
# HC  4 2 4 2 2
# CSD 5 0 0 6 1

# save 'all' as a new df that can be modified:
meta <- all

# cohorts 2 & 4 are confounded (no CSD in either cohort); drop
meta  %<>%  filter(Cohort != "2")
meta  %<>%  filter(Cohort != "4")

# Convert things to factors:
all$Group <- factor(all$Group, 
                    levels = c("HC", "CSD"))
all$Cohort <- factor(all$Cohort)

meta$Group <- factor(meta$Group, 
                     levels = c("HC", "CSD"))
meta$Cohort <- factor(meta$Cohort)

# get before and after counts:
all %>% dplyr::group_by(Group) %>% summarise(count = n()) 
# Group    count
# 1 HC       14
# 2 CSD      12

meta %>% dplyr::group_by(Group) %>% summarise(count = n()) 
# Group    count
# 1 HC        8
# 2 CSD      12

# Transform data ----------------------------------------------------------

# create USM_binary column
meta %<>% mutate(USM_binary = 
                         if_else(USM_pref > 0, 1, 0))

# IHC neutrophil counts fail normality; transform the values
# first select the columns to transform:
cols <- meta %>% select(22:25) %>% colnames()
# store as new dataframe because original column data will be overwritten
meta.trans <- meta %>% mutate_at("Men_total", log) %>% 
  mutate_at(cols, log10)

# one formal outlier in 3xmeninges columns, remove:
scrub <- meta.trans %>% filter(Sample_ID != "4")
# CSD very high


# Linear mixed modeling: OF (meninges) ----------------------------------------
# change settings for 2x2 graphing window
par(mfrow=c(2,2))

#' `meninges total`
# run model:
model = lm(OF_cross_10m ~ Men_total + Cohort, 
           data = scrub, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


#' `meninges - parenchyma (>10um from blood vessel)`
# run model:
model = lm(OF_cross_10m ~ Men_parenchymal + Cohort, 
           data = scrub, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


#' `meninges - abluminal (<=10um from blood vessel)`
# run model:
model = lm(OF_cross_10m ~ Men_abluminal + Cohort, 
           data = scrub, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model)
# get plain english summary:
report(model)


#' `meninges - vascular`
# run model:
model = lm(OF_cross_10m ~ Men_vascular + Cohort, 
           data = scrub, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)

# Mixed logistic regression - USM (meninges) -------------------------------

#' `meninges total`
# run model:
glm.model <- glm(USM_binary ~ Men_total + Cohort, 
                 data = scrub, 
                 family = "binomial",
                 na.action = na.exclude)
# check output
summary(glm.model)
# check model assumptions
plot(glm.model) 
# get plain english summary:
report(glm.model)


#' `meninges - parenchyma (>10um from blood vessel)`
# run model:
glm.model <- glm(USM_binary ~ Men_parenchymal + Cohort, 
                 data = scrub, 
                 family = "binomial",
                 na.action = na.exclude)
# check output
summary(glm.model)
# check model assumptions
plot(glm.model) 
# get plain english summary:
report(glm.model)


#' `meninges - abluminal (<=10um from blood vessel)`
# run model:
glm.model <- glm(USM_binary ~ Men_abluminal + Cohort, 
                 data = scrub, 
                 family = "binomial",
                 na.action = na.exclude)
# check output
summary(glm.model)
# check model assumptions
plot(glm.model) 
# get plain english summary:
report(glm.model)

#' `meninges - vascular`
#' run model:
glm.model <- glm(USM_binary ~ Men_vascular + Cohort, 
                 data = scrub, 
                 family = "binomial",
                 na.action = na.exclude)
# check output
summary(glm.model)
# check model assumptions
plot(glm.model) 
# get plain english summary:
report(glm.model)

# Linear mixed modeling - OF (blood) -----------------------------------------

#' `blood - monocytes`
#' run model:
model = lm(OF_cross_10m ~ Blood_mono + Cohort, 
           data = meta, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model)
# get plain english summary:
report(model)


#' `blood - neutrophils`
#' run model:
model = lm(OF_cross_10m ~ Blood_NPs + Cohort, 
           data = meta, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


#' `blood - T cells`
#' run model:
model = lm(OF_cross_10m ~ Blood_T + Cohort, 
           data = meta, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


#' `blood - B cells`
#' run model:
model = lm(OF_cross_10m ~ Blood_B + Cohort, 
           data = meta, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model)
# get plain english summary:
report(model)

# Logistic regression - USM (blood)-----------------------------------------------

#' `blood - monocytes`
#' run model:
glm.model <- glm(USM_binary ~ Blood_mono  + Cohort, 
                 data = meta, 
                 family = "binomial",
                 na.action = na.exclude)
# check output
summary(glm.model)
# check model assumptions
plot(glm.model) 
# get plain english summary:
report(glm.model)


#' `blood - neutrophils`
#' run model:
glm.model <- glm(USM_binary ~ Blood_NPs  + Cohort, 
                 data = meta, 
                 family = "binomial",
                 na.action = na.exclude)
# check output
summary(glm.model)
# check model assumptions
plot(glm.model)
# get plain english summary:
report(glm.model)


#' `blood - T cells`
#' run model:
glm.model <- glm(USM_binary ~ Blood_T  + Cohort, 
                 data = meta, 
                 family = "binomial",
                 na.action = na.exclude)
# check output
summary(glm.model)
# check model assumptions
plot(glm.model) 
# get plain english summary:
report(glm.model)


#' `blood - B cells`
#' run model:
glm.model <- glm(USM_binary ~ Blood_B  + Cohort, 
                 data = meta, 
                 family = "binomial",
                 na.action = na.exclude)
# check output
summary(glm.model)
# check model assumptions
plot(glm.model) 
# get plain english summary:
report(glm.model)

# resets 2x2 graphing window
dev.off()


# dot & whiskers graph set up-----------------------------------------------------------------

library(stringr)
library(dotwhisker)


# import Arial font
font_add("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_auto()


# select output directory for graphs
setwd("/output")



# Dot & whiskers plot - OF & meninges ----------------------------------------
final <- scrub
colnames(final)

# save variables to be graphed:
ordered_vars <- c("Men_nonvascular",
                  "Men_vascular",
                  "Blood_NPs")

#' create tidy df for linear model - `OF cross`
m <- list()
m[[1]] <- lm(OF_cross_10m ~ Men_nonvascular + Cohort, data = final, 
             na.action = na.exclude)

of_df <- m[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1")

# loop for remaining variables
for (i in 2:3) {
  m[[i]] <- lm(str_replace("OF_cross_10m ~ meow + Cohort", "meow", 
                           ordered_vars[i]), 
               data = final, na.action = na.exclude)
  of_df <- rbind(of_df, 
                 m[[i]] %>% 
                   broom::tidy() %>% 
                   by_2sd(final) %>% 
                   mutate(model = paste("Model", i)))
}

# don't care about the cohort or intercept terms, filter:
of_df_mod <- of_df %>% filter(!grepl('Cohort*', term)) %>% 
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
    c(
      Men_parenchymal = "meningeal\nneutrophils (>10\u03bcm)",
      Men_abluminal = "meningeal\nneutrophils (\u226410\u03bcm)",
      Men_vascular = "meningeal\nneutrophils (iv+)"
    )
  ) +
  theme_pubr() + xlab("Coefficient estimate (OF)") +
  theme(text = element_text(size = 28, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441", "#FD4441", "#FD4441")) +
  annotate(geom = "text", label = "*", x = -18.3, y = 3, size = 12) +
  annotate(geom = "text", label = "*", x = -12, y = 2.15, size = 12) +
  xlim(-30,5)
# ggsave("sf5_OF_coefficients.pdf")

# multiple comparisons correction:
p.adjust(of_df_mod$p.value, method="BH")

# add p.adj values to tibble:
of_df_mod$p.adj <- p.adjust(of_df_mod$p.value, method="BH")

# save tibble for future table creation:
# of_df_mod %>% write_csv("OFvMen_neut_mlm_stats.csv")


# Dot & whiskers plot - USM & meninges ----------------------------------------
# dataframe and ordered_vars are the same as above section

#' Create tidy df of coefficient estimates - `USM`
m2 <- list()

m2[[1]] <- glm(USM_binary ~ Men_total + Cohort, data = final, 
               family = binomial, na.action = na.exclude)

usm_df <- m2[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1") 

# loop for remaining variables
for (i in 2:4) {
  m2[[i]] <- glm(str_replace("USM_binary ~ meow + Cohort", "meow", 
                             ordered_vars[i]), 
                 data = final, family = binomial, na.action = na.exclude)
  usm_df <- rbind(usm_df, 
                  m2[[i]] %>% 
                    broom::tidy() %>% 
                    by_2sd(final) %>% 
                    mutate(model = paste("Model", i))) 
  
}

# don't care about the cohort or intercept terms, filter:
usm_df_mod <- usm_df %>% filter(!grepl('Cohort*', term)) %>% 
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
    c(
      Men_parenchymal = "meningeal\nneutrophils (>10\u03bcm)",
      Men_abluminal = "meningeal\nneutrophils (\u226410\u03bcm)",
      Men_vascular = "meningeal\nneutrophils (iv+)"
    )
  ) +
  theme_pubr() + xlab("Coefficient estimate (USM)") +
  theme(text = element_text(size = 28, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441", "#FD4441", "#FD4441")) +
  xlim(-14.5,1)
# ggsave("sf5_USM_coefficients.pdf")

# multiple comparisons correction:
p.adjust(usm_df_mod$p.value, method="BH")


# add p.adj values to tibble:
usm_df_mod$p.adj <- p.adjust(usm_df_mod$p.value, method="BH")

# save tibble for future table creation:
# usm_df_mod %>% write_csv("USMvMen_neut_glm_stats.csv")


# Dot & whiskers plot - OF & blood ----------------------------------------
final <- meta
colnames(final)

# save variables to be graphed:
ordered_vars <- c("Blood_NPs",
                  "Blood_mono",
                  "Blood_B",
                  "Blood_T")

#' create tidy df for linear model - `OF cross`
m <- list()
m[[1]] <- lm(OF_cross_10m ~ Blood_NPs + Cohort, data = final, 
             na.action = na.exclude)

of_df <- m[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1")

# loop for remaining variables
for (i in 2:4) {
  m[[i]] <- lm(str_replace("OF_cross_10m ~ meow + Cohort", "meow", 
                           ordered_vars[i]), 
               data = final, na.action = na.exclude)
  of_df <- rbind(of_df, 
                 m[[i]] %>% 
                   broom::tidy() %>% 
                   by_2sd(final) %>% 
                   mutate(model = paste("Model", i)))
}

# don't care about the cohort or intercept terms, filter:
of_df_mod <- of_df %>% filter(!grepl('Cohort*', term)) %>% 
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
    c(
      Blood_NPs = "%blood neutrophils",
      Blood_mono = "%blood monocytes",
      Blood_B = "%blood B cells",
      Blood_T = "%blood T cells"
    )
  ) +
  theme_pubr() + xlab("Coefficient estimate (OF)") +
  theme(text = element_text(size = 28, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441", "#366CFC", "#B663E0", "#39BF6A"))
# ggsave("sf5_OFblood_coefficients.pdf")

# multiple comparisons correction:
p.adjust(of_df_mod$p.value, method="BH")


# add p.adj values to tibble:
of_df_mod$p.adj <- p.adjust(of_df_mod$p.value, method="BH")

# save tibble for future table creation:
# of_df_mod %>% write_csv("OFvblood_mlm_stats.csv")


# Dot & whiskers plot - USM & blood ----------------------------------------
# dataframe and ordered_vars are the same as above section

#' Create tidy df of coefficient estimates - `USM`
m2 <- list()

m2[[1]] <- glm(USM_binary ~ Blood_NPs + Cohort, data = final, 
               family = binomial, na.action = na.exclude)

usm_df <- m2[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1") 

# loop for remaining variables
for (i in 2:4) {
  m2[[i]] <- glm(str_replace("USM_binary ~ meow + Cohort", "meow", 
                             ordered_vars[i]), 
                 data = final, family = binomial, na.action = na.exclude)
  usm_df <- rbind(usm_df, 
                  m2[[i]] %>% 
                    broom::tidy() %>% 
                    by_2sd(final) %>% 
                    mutate(model = paste("Model", i))) 
  
}

# don't care about the cohort or intercept terms, filter:
usm_df_mod <- usm_df %>% filter(!grepl('Cohort*', term)) %>% 
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
    c(
      Blood_NPs = "%blood neutrophils",
      Blood_mono = "%blood monocytes",
      Blood_B = "%blood B cells",
      Blood_T = "%blood T cells"
    )
  ) +
  theme_pubr() + xlab("Coefficient estimate (USM)") +
  theme(text = element_text(size = 28, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441", "#366CFC", "#B663E0", "#39BF6A"))
# ggsave("sf5_USMblood_coefficients.pdf")

# multiple comparisons correction:
p.adjust(usm_df_mod$p.value, method="BH")


# add p.adj values to tibble:
usm_df_mod$p.adj <- p.adjust(usm_df_mod$p.value, method="BH")

# save tibble for future table creation:
# usm_df_mod %>% write_csv("USMvblood_glm_stats.csv")

