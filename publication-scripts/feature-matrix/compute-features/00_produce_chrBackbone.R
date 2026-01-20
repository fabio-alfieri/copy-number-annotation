gc(full=T)

packages <- c("parallel")

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

wd <- 'path/to/GitHub/copy-number-annotation/'
setwd(wd)

chr_info_path <- "data/centromeres_and_chromosomes/chr_info_h19.txt"
backbone_outpath <- 'data/centromeres_and_chromosomes/chr_backbone.rds'

chr_info <-
  read.table(chr_info_path, header = TRUE)

fixed_bin_length <- 100000

chr_backbone_1mbp <- mclapply(1:22, mc.cores = 4, function(chr){
  n_bins <-
    as.integer(chr_info[chr_info$Chromosome == paste0("chr", chr),]$Length /
                 fixed_bin_length)
  length_bin <-
    as.integer(chr_info[chr_info$Chromosome == paste0("chr", chr),]$Length /
                 n_bins)
  
  end <- length_bin
  start_bin <- 1
  end_bin <- length_bin
  
  chr_backbone <- data.frame()
  for(bin in 1:n_bins){
    chr_backbone <- rbind(chr_backbone,
                          cbind(chr,
                                bin,
                                start_bin,
                                end_bin))
    start_bin <- start_bin+length_bin
    end_bin <- end_bin+length_bin
  }
  return(chr_backbone)
})


add.col <- function(df, segment_length) {
  df <- cbind(df,
              resize = rep(
                1:round(dim(df)[1] / segment_length + 0.5),
                each = round(dim(df)[1] / round(
                  dim(df)[1] / segment_length + 0.5
                ) + 0.5)
              )[1:nrow(df)])
  return(df)
}
merge.bins <- function(df) {
  df2 <- data.frame()
  for (i in levels(factor(df$resize))) {
    start <- df[df$resize == i, ]$start_bin[1]
    end <-
      df[df$resize == i, ]$end_bin[length(df[df$resize == i, ]$start_bin)]
    
    df2 <- rbind.data.frame(
      df2,
      c(
        i,
        start,
        end
      ),
      stringsAsFactors = FALSE
    )
  }
  
  colnames(df2) <- c(
    "bin",
    "bin_start",
    "bin_end"
  )
  
  return(df2)
}

chr_backbone <- list()
chr_backbone[[paste0(1,'Mbp')]] <- chr_backbone_1mbp

for(chr in 1:22){
  bin_gene <- chr_backbone[['1Mbp']][[chr]]
  
  for(segment_length in 2:50){
    bin_gene_sg <- add.col(bin_gene, segment_length)
    bin_gene_merged <- merge.bins(bin_gene_sg)
    chr_backbone[[paste0(segment_length,'Mbp')]][[chr]] <- cbind(chr = chr,bin_gene_merged)
    
  }
}

saveRDS(chr_backbone, file = backbone_outpath)
