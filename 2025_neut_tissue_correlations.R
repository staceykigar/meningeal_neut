# followed tutorial from: 
# https://ryjohnson09.netlify.app/post/how-to-make-a-heatmap-in-r/

# Set-up ------------------------------------------------------------------

# load libraries
library(tidyverse)
library(magrittr)
library(gplots) # makes pretty correlation plots
library(viridis) #nice color palette
library ("Hmisc") # calculates significance within correlation matrix
library(corrplot) # makes correlation plots


# import data
setwd("~/~r_projects/Meninges_neuts/clean/")
df <- read_csv("np10_neuts.csv")


# set working directory for graph image output:
setwd("~/~r_projects/Meninges_neuts/results/2025")


# prep for neutrophil graph -----------------------------------------------

# create neut specific df:
df1 <- df %>% select(!contains("MO")) %>% select(!contains("L"))

# rename columns so pretty for graphing:
df1 %<>% rename(`meninges (iv-)` = men_nvNP,
               `meninges (iv+)` = men_ivNP, 
               blood = bl_NP,
               spleen = sp_NP,
               skull = sk_NP,
               tibia = tib_NP)

# create correlation matrix ---------------------------------------------

# generate correlation matrix with p values
# exclude descriptive variables
M1 <- rcorr (as.matrix (df1[,-(1:6)] )) 


# graph data -----------------------------------------------------------


# set file name
png(filename = "np10_neutrophil_corr_inferno.png", 
    width = 7, height = 5, units = "in", res = 300)

# make graph
heatmap.2(M1$r, # contains correlation values from correlation matrix
          margins = c(10,10),
          density.info = "none",
          trace = "none",
          dendrogram = "none",
          colsep=1:nrow(M1$r),
          rowsep=1:nrow(M1$r),
          sepcolor = "black",
          col = viridis::viridis_pal(option = "B"),
          breaks = c(seq(-0.03,0.4,length=1),
                     seq(0.44, 1, length=14))) # alters color scale to capture the range of correlation values present in the monocyte matrix (for consistency between graphs)

# turn off file save
dev.off()


# prep for monocyte graph -------------------------------------------------

# create mono specific df:
df2 <- df %>% select(!contains("NP")) %>% select(!contains("L"))

# want to match the order of clustering from the neut graph, so 
# change column order:
df2 %<>% select(1:6, tib_MO, sk_MO, men_nvMO,
                sp_MO, bl_MO, men_ivMO )


df2 %<>% select(1:6, men_ivMO, bl_MO,  sp_MO,
                men_nvMO, sk_MO, tib_MO)

# rename columns so pretty for graphing:
df2 %<>% rename(`meninges (iv-)` = men_nvMO,
                `meninges (iv+)` = men_ivMO, 
                blood = bl_MO,
                spleen = sp_MO,
                skull = sk_MO,
                tibia = tib_MO)



# generate correlation matrix: monocytes --------------------------------------

# generate correlation matrix with p values
# exclude descriptive variables
M2 <- rcorr (as.matrix (df2[,-(1:6)] )) 



# graph data: monocytes -----------------------------------------------------


# set file name
png(filename = "np10_monocyte_corr_inferno.png", 
    width = 7, height = 5, units = "in", res = 300)

# make graph
heatmap.2(M2$r, # contains correlation values from correlation matrix
          margins = c(10,10),
          density.info = "none",
          trace = "none",
          dendrogram = "none",
          colsep=1:nrow(M2$r),
          rowsep=1:nrow(M2$r),
          sepcolor = "black",
          col = viridis::viridis_pal(option = "B"),
          symbreaks = F, # prevents default centering on 0
          symkey = F, # prevents centering the color key on 0
          breaks = c(seq(-0.03,0.4,length=1),
                     seq(0.44, 1, length=14)),
          Rowv = F, # prevents reordering of the rows 
          Colv = F, # prevents reordering of the columns
          revC = T) # flips the data on the diagonal to match neut data

# turn off file save
dev.off()


# create reference heatmap for p values --------------------------------

#' corrplot doesn't offer all the same visualization options, but it does
#' permit significance assessment. generate another heatmap as reference so
#' significance can be added to the nicer graph elsewhere

#' Multiple comparisons testing, 15 unique tests
#' Bonferroni correction: 0.05/15 of unique tests
0.05/15


# generate significance plot as reference for prettier plot:
corrplot(M1$r, order = "hclust", tl.col = "black", tl.srt = 45, tl.cex = 1,
         p.mat = M1$P, sig.level = 0.003333333, insig = "label_sig") 

# generate significance plot as reference for prettier plot:
corrplot(M2$r, order = "hclust", tl.col = "black", tl.srt = 45, tl.cex = 1,
         p.mat = M2$P, sig.level = 0.003333333, insig = "label_sig") 


# prep for lymphocyte graph -------------------------------------------------

# create mono specific df:
df3 <- df %>% select(!contains("NP")) %>% select(!contains("MO"))

# want to match the order of clustering from the neut graph, so 
# change column order:
df3 %<>% select(1:6, tib_L, sk_L, men_nvL,
                sp_L, bl_L, men_ivL )


df3 %<>% select(1:6, men_ivL, bl_L,  sp_L,
                men_nvL, sk_L, tib_L)

# rename columns so pretty for graphing:
df3 %<>% rename(`meninges (iv-)` = men_nvL,
                `meninges (iv+)` = men_ivL, 
                blood = bl_L,
                spleen = sp_L,
                skull = sk_L,
                tibia = tib_L)



# generate correlation matrix: lymphocytes --------------------------------------

# generate correlation matrix with p values
# exclude descriptive variables
M3 <- rcorr (as.matrix (df3[,-(1:6)] )) 



# graph data: lymphocytes -----------------------------------------------------


# set file name
png(filename = "np10_lymphocyte_corr_inferno.png", 
    width = 7, height = 5, units = "in", res = 300)

# make graph
heatmap.2(M3$r, # contains correlation values from correlation matrix
          margins = c(10,10),
          density.info = "none",
          trace = "none",
          dendrogram = "none",
          colsep=1:nrow(M3$r),
          rowsep=1:nrow(M3$r),
          sepcolor = "black",
          col = viridis::viridis_pal(option = "B"),
          symbreaks = F, # prevents default centering on 0
          symkey = F, # prevents centering the color key on 0
          breaks = c(seq(-0.03,0.4,length=1),
                     seq(0.44, 1, length=14)),
          Rowv = F, # prevents reordering of the rows 
          Colv = F, # prevents reordering of the columns
          revC = T) # flips the data on the diagonal to match neut data

# turn off file save
dev.off()
