# Author: Stacey L. Kigar
# 20250317

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
setwd("~/~r_projects/Meninges_neuts/clean/")

all <- read_csv("2025_nat_comm_wildtype.csv")

# make batch column
all$batch <- paste0(all$experiment, "_", all$cohort)


# Assess independence of main predictors ----------------------------------
# i know i want to control for batch in my linear modeling, so compare
# group distribution within batch
table(all$group, all$batch)

#      Exo14NP_1 Exo14NP_2 NP-1_1 NP-1_2 NP-1_4 NP-13_1 NP-13_2 NP-3_1 NP-3_3 NP-3_4 NP-4_1 NP-4_3
# CSD         3         3      3      3      0       4       0      4      2      1      2      2
# HC          3         2      0      0      8       4       4      2      2      2      2      2
# 
# NP-4_4 NP-4_5 NP-4_6 NP-4_7 NP-4_8 NP-5_3 NP-5_4 NP-8_1 NP10_1 NP10_2 NP10_3 NP10_4
# CSD      2      2      1      1      1      0      2      2      2      2      3      2
# HC       2      2      2      0      2      2      1      2      2      2      2      2


# drop batches that are confounded for group
all  %<>%  filter(batch != "NP-1_1")
all  %<>%  filter(batch != "NP-1_2")
all  %<>%  filter(batch != "NP-1_4")
all  %<>%  filter(batch != "NP-13_2")
all  %<>%  filter(batch != "NP-4_7")
all  %<>%  filter(batch != "NP-5_3")


# remove groups potentially onfounded by injeciion stress:
all  %<>%  filter(batch != "NP-8_1")

# Wrangle data ------------------------------------------------------------

# make things factors:
all$group <- factor(all$group, 
                    levels = c("HC", "CSD"))
all$batch <- factor(all$batch)
all$study <- factor(all$experiment)
all$cohort <- factor(all$cohort)


# make 'unique ID column':
all$SampleID <- paste0(all$experiment, "_", all$sample)

colnames(all)

# Transform  --------------------------------------------------------------

cols <- all %>% select(perNP_m_nv, perMO_m_nv, perB_m_nv, perT_m_nv, 
                       perDC_m_nv, perNP_m_iv,  perMO_m_iv, perB_m_iv,
                       perT_m_iv, perDC_m_iv, perNP_b, perMO_b, perB_b,
                       perT_b, perDC_b) %>% colnames()

# store as new data frame because original column data will be overwritten
trans <- all %>% mutate_at(cols, sqrt) 


# Linear mixed modeling: OF  -----------------------------------------------

# set viewing window for looking at residual graphs:
par(mfrow=c(2,2))

# simplify dataframe by dropping levels for NA OF data:
final <- trans %>% drop_na(OF_m_05) %>% droplevels()


#' `blood - neutrophils`
#' run model:
model = lm(OF_m_05 ~ perNP_b + batch, 
           data = final, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
report(model)


#' `meninges - iv+ neutrophils`
#' run model:
model = lm(OF_m_05 ~ perNP_m_iv + batch, 
           data = final, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)

#' `meninges - iv- neutrophils`
#' run model:
model = lm(OF_m_05 ~ perNP_m_nv + batch, 
           data = final, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


# Linear mixed modeling: LD  -----------------------------------------------

# set viewing window for looking at residual graphs:
par(mfrow=c(2,2))

# simplify dataframe by dropping levels for NA LD data:
final <- trans %>% drop_na(LD_cross) %>% droplevels()


#' `blood - neutrophils`
#' run model:
model = lm(LD_cross ~ perNP_b + batch, 
           data = final, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
report(model)


#' `meninges - iv+ neutrophils`
#' run model:
model = lm(LD_cross ~ perNP_m_iv + batch, 
           data = final, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)

#' `meninges - iv- neutrophils`
#' run model:
model = lm(LD_cross ~ perNP_m_nv + batch, 
           data = final, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


# Linear mixed modeling: SI  -----------------------------------------------

# set viewing window for looking at residual graphs:
par(mfrow=c(2,2))

# simplify dataframe by dropping levels for NA SI data:
final <- trans %>% drop_na(SI_approach_610) %>% droplevels()


#' `blood - neutrophils`
#' run model:
model = lm(SI_approach_610 ~ perNP_b + batch, 
           data = final, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
report(model)


#' `meninges - iv+ neutrophils`
#' run model:
model = lm(SI_approach_610 ~ perNP_m_iv + batch, 
           data = final, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)

#' `meninges - iv- neutrophils`
#' run model:
model = lm(SI_approach_610 ~ perNP_m_nv + batch, 
           data = final, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)

# Mixed logistic regression - USM   -----------------------------------------

# simplify dataframe by dropping levels for NA OF data:
final <- trans %>% drop_na(USM_mark) %>% droplevels()


#' `blood - neutrophils`
#' run model:
model = glm(USM_mark ~ perNP_b + batch, data = final, 
            family = binomial, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


#' `meninges - iv+ neutrophils`
#' run model:
model = glm(USM_mark ~ perNP_m_iv + batch, data = final, 
            family = binomial, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


#' `meninges - iv- neutrophils`
#' run model:
model = glm(USM_mark ~ perNP_m_nv + batch, data = final, 
            family = binomial, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model)
# get plain english summary:
report(model)


# linear model for blood vs meningeal neuts -------------------------------

# simplify dataframe by dropping levels for NA OF data:
final <- trans %>% drop_na(perNP_m_nv) %>% drop_na(perNP_b) %>% droplevels()

# generate model:
model <- final %>%  
  lm(formula = perNP_m_nv ~ perNP_b * group + batch, na.action = na.exclude)
# check output
summary(model)
# check model assumptions
plot(model) 
# get plain english summary:
report(model)


# Close split plotting window
dev.off()

# interaction plot --------------------------------------------------------
#' library(interactions)
#' 
#' setwd("~/~r_projects/Meninges_neuts/results/2025/")
#' 
#' #' `meningeal neut (nv) vs blood neuts`
#' 
#' # generate model: 
#' mNPs_vbNPs_study <- final %>%  
#'   lm(formula = perNP_m_nv ~ perNP_b * group + batch)
#' 
#' 
#' # plot and save
#' interactions::interact_plot(mNPs_vbNPs_study, pred = perNP_b, 
#'                             modx = group, plot.points = T,
#'                             line.thickness = 2,
#'                             point.size = 12, partial.residuals = T) +
#'   theme_classic() + 
#'   xlab("\n%blood neutrophils") + 
#'   ylab("%(iv-) meningeal neutrophils\n") +
#'   theme(text = element_text(size = 36, 
#'                             family = "Arial", colour = "black")) +
#'   scale_color_manual(values = c('#636463', '#FD4441'))
#' ggsave("f1_bNPs_vmNPsnv_lm.pdf")
#' 
#' 
#' # get parameters of line:
#' interactions::sim_slopes(mNPs_vbNPs_study, pred = perNP_b, 
#'            modx = group)
#' 
#' # Slope of Blood_NPs when Group = CSD: 
#' #   
#' #   Est.   S.E.   t val.      p
#' # ------- ------ -------- ------
#' #   -1.56   1.13    -1.38   0.18
#' # 
#' # Slope of Blood_NPs when Group = HC: 
#' #   
#' #   Est.   S.E.   t val.      p
#' # ------- ------ -------- ------
#' #   -0.78   0.86    -0.91   0.37

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

final <- trans %>% drop_na(OF_m_05) %>% droplevels()
colnames(final)
 
# save variables to be graphed:
ordered_vars <- c("perNP_m_nv",
                  "perNP_m_iv",
                  "perNP_b")

#' create tidy df for linear model - `OF cross`
m <- list()
m[[1]] <- lm(OF_m_05 ~ perNP_m_nv + batch, data = final, 
             na.action = na.exclude)

of_df <- m[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1")

# loop for remaining variables
for (i in 2:3) {
  m[[i]] <- lm(str_replace("OF_m_05 ~ meow + batch", "meow", 
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

# multiple comparisons correction:
p.adjust(of_df_mod$p.value, method="BH")

# add p.adj values to tibble:
of_df_mod$p.adj <- p.adjust(of_df_mod$p.value, method="BH")

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
    c(perNP_m_nv = "%meningeal \nneutrophils (iv-)",
      perNP_m_iv = "%meningeal \nneutrophils (iv+)",
      perNP_b = "%blood \nneutrophils" )
  ) +
  theme_pubr() + xlab("Coefficient estimate (OF)") +
  theme(text = element_text(size = 28, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441", "#FD4441", "#FD4441")) +
  annotate(geom = "text", label = "*", x = -3.113539, y = 2.09, size = 12) +
  annotate(geom = "text", label = "**", x = -2.803703, y = 1.24, size = 12) +
  xlim(-6,1)
ggsave("WT_OF_coefficients.pdf")

# save tibble for future table creation:
of_df_mod %>% write_csv("WT_OF_neut_mlm_stats.csv")


# Dot & whiskers plot - OF (all cells) ---------------------------------------

final <- trans %>% drop_na(OF_m_05) %>% droplevels()
colnames(final)

# save variables to be graphed:
ordered_vars <- cols

#' create tidy df for linear model - `OF cross`
m <- list()
m[[1]] <- lm(OF_m_05 ~ perNP_m_nv + batch, data = final, 
             na.action = na.exclude)

of_df <- m[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1")

# loop for remaining variables
for (i in 2:15) {
  m[[i]] <- lm(str_replace("OF_m_05 ~ meow + batch", "meow", 
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

# multiple comparisons correction:
p.adjust(of_df_mod$p.value, method="BH")

# add p.adj values to tibble:
of_df_mod$p.adj <- p.adjust(of_df_mod$p.value, method="BH")

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
    c(perMO_m_nv = "%men. monocytes (iv-)",
      perNP_m_nv = "%men. neutrophils (iv-)",
      perT_m_nv = "%men. T cells (iv-)",
      perB_m_nv = "%men. B cells (iv-)",
      perDC_m_nv = "%men. dendritic cells (iv-)",
      perMO_m_iv = "%men. monocytes (iv+)",
      perNP_m_iv = "%men. neutrophils (iv+)",
      perT_m_iv = "%men. T cells (iv+)",
      perB_m_iv = "%men. B cells (iv+)",
      perDC_m_iv = "%men. dendritic cells (iv+)",
      perMO_b = "%blood monocytes",
      perNP_b = "%blood neutrophils",
      perT_b = "%blood T cells",
      perB_b = "%blood B cells",
      perDC_b = "%blood dendritic cells")
  ) +
  theme_pubr() + xlab("Coefficient estimate (OF)") +
  theme(text = element_text(size = 12, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441","#366CFC","#B663E0", "#39BF6A", "#FD8F24",
                              "#FD4441","#366CFC","#B663E0", "#39BF6A", "#FD8F24",
                              "#FD4441","#366CFC","#B663E0", "#39BF6A", "#FD8F24")) + 
  annotate(geom = "text", label = "*", x = -2.8037027, y = 4.4, size = 10) 
ggsave("WT_OF_coefficients_allcells.pdf")

# save tibble for future table creation:
of_df_mod %>% write_csv("WT_OF_allcells_mlm_stats.csv")



# Dot & whiskers plot - LD ----------------------------------------------------

final <- trans %>% drop_na(LD_cross) %>% droplevels()
colnames(final)

# save variables to be graphed:
ordered_vars <- c("perNP_m_nv",
                  "perNP_m_iv",
                  "perNP_b")

#' create tidy df for linear model - `LD cross`
m <- list()
m[[1]] <- lm(LD_cross ~ perNP_m_nv + batch, data = final, 
             na.action = na.exclude)

ld_df <- m[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1")

# loop for remaining variables
for (i in 2:3) {
  m[[i]] <- lm(str_replace("LD_cross ~ meow + batch", "meow", 
                           ordered_vars[i]), 
               data = final, na.action = na.exclude)
  ld_df <- rbind(ld_df, 
                 m[[i]] %>% 
                   broom::tidy() %>% 
                   by_2sd(final) %>% 
                   mutate(model = paste("Model", i)))
}

# don't care about the cohort or intercept terms, filter:
ld_df_mod <- ld_df %>% filter(!grepl('batch*', term)) %>% 
  filter(!grepl('Intercept*', term))

# multiple comparisons correction:
p.adjust(ld_df_mod$p.value, method="BH")

# add p.adj values to tibble:
ld_df_mod$p.adj <- p.adjust(ld_df_mod$p.value, method="BH")

#' make the graph - `LD cross`
dwplot(ld_df_mod, 
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
    c(perNP_m_nv = "%meningeal \nneutrophils (iv-)",
      perNP_m_iv = "%meningeal \nneutrophils (iv+)",
      perNP_b = "%blood \nneutrophils" )
  ) +
  theme_pubr() + xlab("Coefficient estimate (LD)") +
  theme(text = element_text(size = 28, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441", "#FD4441", "#FD4441"))
ggsave("WT_LD_coefficients.pdf")

# save tibble for future table creation:
ld_df_mod %>% write_csv("WT_ld_neut_mlm_stats.csv")


# Dot & whiskers plot - LD (all cells) ----------------------------------------

final <- trans %>% drop_na(LD_cross) %>% droplevels()

# save variables to be graphed:
ordered_vars <- cols

#' create tidy df for linear model - `LD cross`
m <- list()
m[[1]] <- lm(LD_cross ~ perNP_m_nv + batch, data = final, 
             na.action = na.exclude)

ld_df <- m[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1")

# loop for remaining variables
for (i in 2:15) {
  m[[i]] <- lm(str_replace("LD_cross ~ meow + batch", "meow", 
                           ordered_vars[i]), 
               data = final, na.action = na.exclude)
  ld_df <- rbind(ld_df, 
                 m[[i]] %>% 
                   broom::tidy() %>% 
                   by_2sd(final) %>% 
                   mutate(model = paste("Model", i)))
}

# don't care about the cohort or intercept terms, filter:
ld_df_mod <- ld_df %>% filter(!grepl('batch*', term)) %>% 
  filter(!grepl('Intercept*', term))

# multiple comparisons correction:
p.adjust(ld_df_mod$p.value, method="BH")

# add p.adj values to tibble:
ld_df_mod$p.adj <- p.adjust(ld_df_mod$p.value, method="BH")

#' make the graph - `LD cross`
dwplot(ld_df_mod, 
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
    c(perMO_m_nv = "%men. monocytes (iv-)",
      perNP_m_nv = "%men. neutrophils (iv-)",
      perT_m_nv = "%men. T cells (iv-)",
      perB_m_nv = "%men. B cells (iv-)",
      perDC_m_nv = "%men. dendritic cells (iv-)",
      perMO_m_iv = "%men. monocytes (iv+)",
      perNP_m_iv = "%men. neutrophils (iv+)",
      perT_m_iv = "%men. T cells (iv+)",
      perB_m_iv = "%men. B cells (iv+)",
      perDC_m_iv = "%men. dendritic cells (iv+)",
      perMO_b = "%blood monocytes",
      perNP_b = "%blood neutrophils",
      perT_b = "%blood T cells",
      perB_b = "%blood B cells",
      perDC_b = "%blood dendritic cells")
  ) +
  theme_pubr() + xlab("Coefficient estimate (LD)") +
  theme(text = element_text(size = 12, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") + 
  scale_color_manual(values=c("#FD4441","#366CFC","#B663E0", "#39BF6A", "#FD8F24", 
                              "#FD4441","#366CFC","#B663E0", "#39BF6A", "#FD8F24",
                              "#FD4441","#366CFC","#B663E0", "#39BF6A", "#FD8F24"))
ggsave("WT_LD_coefficients_allcells.pdf")

# save tibble for future table creation:
ld_df_mod %>% write_csv("WT_ld_allcells_mlm_stats.csv")




# Dot & whiskers plot - SI ----------------------------------------------------

final <- trans %>% drop_na(SI_approach_610) %>% droplevels()
colnames(final)

# save variables to be graphed:
ordered_vars <- c("perNP_m_nv",
                  "perNP_m_iv",
                  "perNP_b")

#' create tidy df for linear model - `social approach`
m <- list()
m[[1]] <- lm(SI_approach_610 ~ perNP_m_nv + batch, data = final, 
             na.action = na.exclude)

si_df <- m[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1")

# loop for remaining variables
for (i in 2:3) {
  m[[i]] <- lm(str_replace("SI_approach_610 ~ meow + batch", "meow", 
                           ordered_vars[i]), 
               data = final, na.action = na.exclude)
  si_df <- rbind(si_df, 
                 m[[i]] %>% 
                   broom::tidy() %>% 
                   by_2sd(final) %>% 
                   mutate(model = paste("Model", i)))
}

# don't care about the cohort or intercept terms, filter:
si_df_mod <- si_df %>% filter(!grepl('batch*', term)) %>% 
  filter(!grepl('Intercept*', term))

# multiple comparisons correction:
p.adjust(si_df_mod$p.value, method="BH")

# add p.adj values to tibble:
si_df_mod$p.adj <- p.adjust(si_df_mod$p.value, method="BH")

#' make the graph - `social approach`
dwplot(si_df_mod, 
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
    c(perNP_m_nv = "%meningeal \nneutrophils (iv-)",
      perNP_m_iv = "%meningeal \nneutrophils (iv+)",
      perNP_b = "%blood \nneutrophils" )
  ) +
  theme_pubr() + xlab("Coefficient estimate (SI)") +
  theme(text = element_text(size = 28, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441", "#FD4441", "#FD4441"))
ggsave("WT_SI_coefficients.pdf")

# save tibble for future table creation:
si_df_mod %>% write_csv("WT_SI_neut_mlm_stats.csv")


# Dot & whiskers plot - SI (all cell types) ------------------------------------

final <- trans %>% drop_na(SI_approach_610) %>% droplevels()
colnames(final)

# save variables to be graphed:
ordered_vars <- cols

#' create tidy df for linear model - `social approach`
m <- list()
m[[1]] <- lm(SI_approach_610 ~ perNP_m_nv + batch, data = final, 
             na.action = na.exclude)

si_df <- m[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1")

# loop for remaining variables
for (i in 2:15) {
  m[[i]] <- lm(str_replace("SI_approach_610 ~ meow + batch", "meow", 
                           ordered_vars[i]), 
               data = final, na.action = na.exclude)
  si_df <- rbind(si_df, 
                 m[[i]] %>% 
                   broom::tidy() %>% 
                   by_2sd(final) %>% 
                   mutate(model = paste("Model", i)))
}

# don't care about the cohort or intercept terms, filter:
si_df_mod <- si_df %>% filter(!grepl('batch*', term)) %>% 
  filter(!grepl('Intercept*', term))

# multiple comparisons correction:
p.adjust(si_df_mod$p.value, method="BH")

# add p.adj values to tibble:
si_df_mod$p.adj <- p.adjust(si_df_mod$p.value, method="BH")

#' make the graph - `social approach`
dwplot(si_df_mod, 
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
    c(perMO_m_nv = "%men. monocytes (iv-)",
      perNP_m_nv = "%men. neutrophils (iv-)",
      perT_m_nv = "%men. T cells (iv-)",
      perB_m_nv = "%men. B cells (iv-)",
      perDC_m_nv = "%men. dendritic cells (iv-)",
      perMO_m_iv = "%men. monocytes (iv+)",
      perNP_m_iv = "%men. neutrophils (iv+)",
      perT_m_iv = "%men. T cells (iv+)",
      perB_m_iv = "%men. B cells (iv+)",
      perDC_m_iv = "%men. dendritic cells (iv+)",
      perMO_b = "%blood monocytes",
      perNP_b = "%blood neutrophils",
      perT_b = "%blood T cells",
      perB_b = "%blood B cells",
      perDC_b = "%blood dendritic cells")
  ) +
  theme_pubr() + xlab("Coefficient estimate (SI)") +
  theme(text = element_text(size = 12, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441","#366CFC","#B663E0", "#39BF6A", "#FD8F24",
                              "#FD4441","#366CFC","#B663E0", "#39BF6A", "#FD8F24",
                              "#FD4441","#366CFC","#B663E0", "#39BF6A", "#FD8F24"))
ggsave("WT_SI_coefficients_allcells.pdf")

# save tibble for future table creation:
si_df_mod %>% write_csv("WT_SI_allcells_mlm_stats.csv")


# Dot & whiskers plot - USM ------------------------------------------------
final <- trans %>% drop_na(USM_mark) %>% droplevels()
colnames(final)

# save variables to be graphed:
ordered_vars <- c("perNP_m_nv",
                  "perNP_m_iv",
                  "perNP_b")

#' Create tidy df of coefficient estimates - `USM`
m2 <- list()

m2[[1]] <- glm(USM_mark ~ perNP_m_nv + batch, data = final, 
               family = binomial, na.action = na.exclude)

usm_df <- m2[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1") 

# loop for remaining variables
for (i in 2:3) {
  m2[[i]] <- glm(str_replace("USM_mark ~ meow + batch", "meow", 
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

# multiple comparisons correction:
p.adjust(usm_df_mod$p.value, method="BH")

# add p.adj values to tibble:
usm_df_mod$p.adj <- p.adjust(usm_df_mod$p.value, method="BH")

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
    c(perNP_m_nv = "%meningeal \nneutrophils (iv-)",
      perNP_m_iv = "%meningeal \nneutrophils (iv+)",
      perNP_b = "%blood \nneutrophils"
    )
  ) +
  theme_pubr() + xlab("Coefficient estimate (USM)") +
  theme(text = element_text(size = 28, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441", "#FD4441", "#FD4441")) +
  annotate(geom = "text", label = "**", x = -4.877790, y = 1.2, size = 12) +
  annotate(geom = "text", label = "*", x = -4.641654, y = 2.09, size = 12) +
  annotate(geom = "text", label = "*", x = -4.565118, y = 2.95, size = 12) 
ggsave("WT_USM_coefficients.pdf")

# save tibble for future table creation:
usm_df_mod %>% write_csv("WT_USM_neut_glm_stats.csv")


# Dot & whiskers plot - USM (all cell types)------------------------------------------
final <- trans %>% drop_na(USM_mark) %>% droplevels()
colnames(final)

# save variables to be graphed:
ordered_vars <- cols

#' Create tidy df of coefficient estimates - `USM`
m2 <- list()

m2[[1]] <- glm(USM_mark ~ perNP_m_nv + batch, data = final, 
               family = binomial, na.action = na.exclude)

usm_df <- m2[[1]] %>% 
  broom::tidy() %>% 
  by_2sd(final) %>% 
  mutate(model = "Model 1") 

# loop for remaining variables
for (i in 2:15) {
  m2[[i]] <- glm(str_replace("USM_mark ~ meow + batch", "meow", 
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

# multiple comparisons correction:
p.adjust(usm_df_mod$p.value, method="BH")

# add p.adj values to tibble:
usm_df_mod$p.adj <- p.adjust(usm_df_mod$p.value, method="BH")

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
    c(perMO_m_nv = "%men. monocytes (iv-)",
      perNP_m_nv = "%men. neutrophils (iv-)",
      perT_m_nv = "%men. T cells (iv-)",
      perB_m_nv = "%men. B cells (iv-)",
      perDC_m_nv = "%men. dendritic cells (iv-)",
      perMO_m_iv = "%men. monocytes (iv+)",
      perNP_m_iv = "%men. neutrophils (iv+)",
      perT_m_iv = "%men. T cells (iv+)",
      perB_m_iv = "%men. B cells (iv+)",
      perDC_m_iv = "%men. dendritic cells (iv+)",
      perMO_b = "%blood monocytes",
      perNP_b = "%blood neutrophils",
      perT_b = "%blood T cells",
      perB_b = "%blood B cells",
      perDC_b = "%blood dendritic cells")
  ) +
  theme_pubr() + xlab("Coefficient estimate (USM)") +
  theme(text = element_text(size = 12, family = "Arial", colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values=c("#FD4441","#366CFC","#B663E0", "#39BF6A", "#FD8F24",
                              "#FD4441","#366CFC","#B663E0", "#39BF6A", "#FD8F24",
                              "#FD4441","#366CFC","#B663E0", "#39BF6A", "#FD8F24")) +
  annotate(geom = "text", label = "*", x = -4.8777895, y = 4.4, size = 10) 
ggsave("WT_USM_coefficients_allcells.pdf")

# save tibble for future table creation:
usm_df_mod %>% write_csv("WT_USM_allcells_glm_stats.csv")


# power calculations ------------------------------------------------------

library(pwr)

### OF  vs blood neuts ####
final <- trans %>% drop_na(OF_m_05) %>% droplevels() %>% 
  drop_na(perNP_b) %>% droplevels()

table(final$group)
# HC CSD 
# 14  15 


model <- lm(OF_m_05 ~ perNP_b + batch, data = final, 
             na.action = na.exclude)
summary(model)

# u = parameters in model = 1 continuous IDV, 1 factor IDV with 5 levels = 6
# n = # of samples (40) - 11 NAs = 29
# v = n - u - 1 = 29 - 6 - 1 = 22

# calculate Cohen's f2 from betas: f^2 = beta^2 / (1 - beta^2)
(-0.6325 )^2 / (1 - (-0.6325 )^2)

pwr.f2.test(u = 6, v = 22, f2 = 0.6668229, sig.level = 0.05, power = NULL)
# 83% power to detect an effect

# how many animals if 1 IDV and 80% power?
pwr.f2.test(u = 1, v = NULL, f2 = 0.6668229, sig.level = 0.05, power = .8)
# v = 11.97925 
# calculate for n = v + u + 1 = 12 + 1 + 1 = 14

### OF  vs iv- meningeal neuts ####
final <- trans %>% drop_na(OF_m_05) %>% droplevels() %>% 
  drop_na(perNP_m_nv) %>% droplevels()

table(final$group)
# HC CSD 
# 19  21 


model <- lm(OF_m_05 ~ perNP_m_nv + batch, data = final, 
            na.action = na.exclude)
summary(model)
plot(model)

# u = parameters in model = 1 continuous IDV, 1 factor IDV with 7 levels = 8
# n = # of samples = 40
# v = n - u - 1 = 40 - 8 - 1 = 31

# calculate Cohen's f2 from betas: f^2 = beta^2 / (1 - beta^2)
(-0.9730 )^2 / (1 - (-0.9730 )^2)

pwr.f2.test(u = 8, v = 31, f2 = 17.77194, sig.level = 0.05, power = NULL)
# 100% power to detect an effect

# do not reject H0

### SI  vs blood neuts ####
final <- trans %>% drop_na(SI_approach_610) %>% droplevels() %>% 
  drop_na(perNP_b) %>% droplevels()

table(final$group)
# HC CSD 
# 17  19 


model <- lm(SI_approach_610 ~ perNP_b + batch, data = final, 
            na.action = na.exclude)
summary(model)
plot(model)

# u = parameters in model = 1 continuous IDV, 1 factor IDV with 8 levels = 9
# n = # of samples (40) - 11 NAs = 9
# v = n - u - 1 = 36 - 9 - 1 = 26

# calculate Cohen's f2 from betas: f^2 = beta^2 / (1 - beta^2)
(-0.72418  )^2 / (1 - (-0.72418 )^2)

pwr.f2.test(u = 9, v = 26, f2 = 1.102769, sig.level = 0.05, power = NULL)
# 98%% power to detect an effect

# do not reject H0

### SI  vs blood neuts ####
final <- trans %>% drop_na(SI_approach_610) %>% droplevels() %>% 
  drop_na(perNP_m_iv) %>% droplevels()

table(final$group)
# HC CSD 
# 22  26 


model <- lm(SI_approach_610 ~ perNP_m_iv + batch, data = final, 
            na.action = na.exclude)
summary(model)
plot(model)

# u = parameters in model = 1 continuous IDV, 1 factor IDV with 8 levels = 9
# n = # of samples (40) - 11 NAs = 9
# v = n - u - 1 = 36 - 9 - 1 = 26

# calculate Cohen's f2 from betas: f^2 = beta^2 / (1 - beta^2)
((-8.7209)^2) / (1 - ((-8.7209)^2)) #-1.013324
# apparently - values mean there's something wrong with the model.

pwr.f2.test(u = 9, v = 26, f2 = -1.013324, sig.level = 0.05, power = NULL)
# 98%% power to detect an effect

# do not reject H0


# test absolute count data for relationships ------------------------------


abs.cols <- trans %>% select(NP_m_nv, MO_m_nv, B_m_nv, T_m_nv, DC_m_nv,
                             NP_m_iv, MO_m_iv, B_m_iv, T_m_iv, DC_m_iv,
                             NP_ul_b, MO_ul_b, B_ul_b, T_ul_b, DC_ul_b) %>% colnames()

# update df
trans <- trans %>% mutate_at(abs.cols, sqrt) 


#### OF vs blood neut ####
final <- trans %>% drop_na(OF_m_05) %>% droplevels() %>% 
  drop_na(NP_ul_b) %>% droplevels()

model <- lm(OF_m_05 ~ NP_ul_b + batch, data = final, na.action = na.exclude)
summary(model)
# 0.0579

#### OF vs iv- meningeal neuts ####
final <- trans %>% drop_na(OF_m_05) %>% droplevels() %>% 
  drop_na(NP_m_nv) %>% droplevels()
# i am left with one batch, so drop from model

model <- lm(OF_m_05 ~ NP_m_nv, data = final, na.action = na.exclude)
summary(model)
# 0.785
