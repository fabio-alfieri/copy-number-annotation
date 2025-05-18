rm(list=ls())
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
toplot.plot <- output.annotation$ampl$toplot # updated
clusters_explained <- output.annotation$ampl$aggregated # updated
centromere_table <- read.table("dev/Data/centomere.tsv", header = T)
load("dev/Data/All_levels_backbonetables.RData")

# rm(output.annotation)
