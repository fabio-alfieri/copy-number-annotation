# rm(list=ls())
gc(full=T)

library(tidyverse)
library(ggplot2)
library(IRanges) 
library(stringr)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggiraph)
library(readr)
library(htmlwidgets)
library(htmltools)
library(ggnewscale)

shap.list <- readRDS("dev/Data/shap_Mid-length_AmplDel.rds")
output.annotation <- readRDS("dev/Data/output_annotation.rds") # NEW
centromere_table <- read.table("dev/Data/centomere.tsv", header = T)
load("dev/Data/All_levels_backbonetables.RData")
out_annot_list <- readRDS("dev/Data/output_annotation.rds")

lis_names <- c("ampl", "del")
paths_vec <- c("dev/Data/pred_ampl.rds", "dev/Data/pred_del.rds")

pred_list <- lapply(X = paths_vec, FUN = readRDS)
names(pred_list) <- lis_names

# toplot.plot_before <- output.annotation$ampl$toplot # updated
# clusters_explained <- output.annotation$ampl$aggregated # updated
# rm(output.annotation)
