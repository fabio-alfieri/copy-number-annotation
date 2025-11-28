library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

extract_gene_names_from_windows <- function(regions, genes, outname, write = FALSE){
  
  grang <- GRanges(regions)
  
  overlaps <- subsetByOverlaps(genes, grang)
  
  gene_names <- annotate::getSYMBOL(x = overlaps$gene_id, data='org.Hs.eg')
  
  names(gene_names) <- NULL
  gene_names <- sort(gene_names)
  
  if(write){
    
    df <- as.data.frame(gene_names)
    write.table(x = df, file = outname, 
                col.names = F, 
                row.names = F, 
                sep = "\t", 
                quote = F)
    
  }
  
  return(gene_names)
  
}


regions
names(regions) <- c("Positive")

genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

gene_lists <- list()

for (idx in seq_along(along.with = region_paths)) {
  
  for(alt in alternatives){
    
    gene_list <- extract_gene_names_from_windows(path, paste0(outname, "_", alt), genes, FALSE, what = alt)
    gene_lists[[paste0(path, "_", alt)]] <- gene_list
  
  }
}

# this snippet takes as input a named list of gene names (a list of vectors)
# and performs GO enrichment analysis over a set of ontologies 
# (Biological Process, Cellular Component, Molecular Function)

pvalueCutoff <- 0.01 # p-value
qvalueCutoff <- 0.1  # q-value, NOTA: q-value secondo clusterProfiler non è la correzione del p-value, ma il minimo p-value adjusted per avere risultati significativi
showCategory <- 15 # top terms to show

go_res <- lapply(X = seq_along(gene_lists), FUN = function(idx){
  
  list_name <- names(gene_lists[idx]) # take sublist name
  list <- gene_lists[[idx]] # take sublist
  
  print(paste0("Processing ", list_name))
  
  gene_df <- bitr(list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # take list of names and convert to entrez ID
  entrez_ids <- gene_df$ENTREZID
  
  ontologies <- c("BP","CC","MF") # define ontologies
  
  to.return <- lapply(X = seq_along(ontologies), FUN = function(idx_ont){
    
    ont <- ontologies[idx_ont] # extract ontology
    print(paste0("Processing ", ont))
    
    go_results <- enrichGO(gene             = entrez_ids, # input genes
                           OrgDb            = org.Hs.eg.db, # database
                           keyType          = "ENTREZID",
                           ont              = ont,
                           pAdjustMethod    = "BH", # p-value correction
                           pvalueCutoff     = pvalueCutoff,
                           qvalueCutoff     = qvalueCutoff)
    
    if (nrow(go_results) > 0) { # if at least one term is enriched
      go_results <- pairwise_termsim(go_results)
      dotplot <- dotplot(go_results, showCategory = showCategory) 
      barplot <- barplot(go_results, showCategory = showCategory)
      emapplot <- emapplot(go_results)
      
      return(list(dotplot = dotplot, #return plots and output
                  barplot = barplot, 
                  emapplot = emapplot, 
                  result = go_results))
    
      } else {
      print(paste0("No enrichment results for ", list_name, " - ", ont))
      return(NULL)
    }
    
  })
  names(to.return) <- ontologies
  return(to.return)
})
names(go_res) <- names(gene_lists)

# does the same but with KEGG
kegg_res <- lapply(X = seq_along(gene_lists), FUN = function(idx){
  
  list_name <- names(gene_lists[idx])
  list <- gene_lists[[idx]]
  
  print(paste0("Processing ", list_name))
  
  gene_df <- bitr(list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- gene_df$ENTREZID
  
  kegg_results <- enrichKEGG(gene          = entrez_ids,
                             organism      = 'hsa',
                             pAdjustMethod = "BH", 
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.1)
  

  if (nrow(kegg_results) > 0) {
    
    kegg_results <- pairwise_termsim(kegg_results)
    dotplot <- dotplot(kegg_results, showCategory = 10)
    barplot <- barplot(kegg_results, showCategory = 10)

    return(list(dotplot = dotplot, 
                barplot = barplot, 
                result = kegg_results))
    
  } else {
    print(paste0("No enrichment results for KEGG - ", list_name))
    return(NULL)
  }
})
names(kegg_res) <- names(gene_lists)


cancer_gene_list <- read.table("/home/ieo7429/Desktop/THESIS_GAB/tables/ONCO_TSG_TABLE_hg19_0.1Mbp.tsv", header = T, sep = "\t", fill = TRUE)
