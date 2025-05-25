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
pred_ampl <- readRDS("dev/Data/pred_ampl.rds")
pred_del <- readRDS("dev/Data/pred_del.rds")

toplot.plot_before <- output.annotation$ampl$toplot # updated
# clusters_explained <- output.annotation$ampl$aggregated # updated
# rm(output.annotation)
