rm(list=ls())
gc(full=T)

packages <- c(
  "stringr", "tidyr", "ggplot2", "fastcluster", "tidyverse"
)

installed <- rownames(installed.packages())

for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

model.outputs_ext <- list.files(path = "/home/ieo7429/Scrivania/results_regressor_gab/InteractomeINSIDER/CSP", full.names = TRUE)
chr_nums <- 1:22
patterns <- c('Mid-length')

for(chr_num_ext in chr_nums){
  
  for(chr_num_int in chr_nums){
    
    model.outputs <- model.outputs_ext[grepl(pattern = paste0("CSP", "_", chr_num_ext, "_", chr_num_int), x = model.outputs_ext)]
    
    if (length(model.outputs) == 0){
      print("skipped")  
      next
    }
    
    for(i in patterns){
      
      load(model.outputs[grep(i, model.outputs)][grep('ampl', model.outputs[grep(i, model.outputs)])])
      model_ampl <- Output.regressor; model_ampl <- model_ampl$ampl_score[[1]]
      
      rm(Output.regressor)
      
      ampl_pred_list <- lapply(X = seq_along(along.with = model_ampl), FUN = function(x){
        model_output <- model_ampl[[x]]
        full_preds <- rbind(model_output$Train.labels, model_output$Test.labels)
      })
      
      mean_preds_ampl <- Reduce(f = rbind, x = ampl_pred_list) %>%
        group_by(bin, Type) %>%
        summarise(
          ampl_score = mean(ampl_score),
          prediction = mean(prediction)) %>% 
        ungroup()
      
      # --- extract chromosome from bin and label train/test ---
      mean_preds_ampl <- mean_preds_ampl %>%
        mutate(chr = sub("_.*", "", bin))
      
      chr_counts <- mean_preds_ampl %>% count(chr)
      
      if (nrow(chr_counts) == 1){
        mean_preds_ampl <- mean_preds_ampl %>%
          mutate(set = "Train + Test (same chr)")
      } else {
        train_chr <- chr_counts$chr[which.max(chr_counts$n)]
        mean_preds_ampl <- mean_preds_ampl %>%
          mutate(set = ifelse(chr == train_chr, "Train", "Test"))
      }
      
      # --- plot ---
      p <- ggplot(mean_preds_ampl, aes(x = ampl_score, y = prediction, color = chr, shape = set)) +
        geom_point(alpha = 0.6) +
        labs(
          title = paste0("CSP_", chr_num_ext, "_", chr_num_int, " - ", i),
          x = "Ampl score", y = "Prediction",
          color = "Chromosome", shape = "Set"
        ) +
        theme_minimal()
      
      #pred_ampl_path <- paste0("/home/ieo7429/Scrivania/results_regressor_gab/predictions/CSP/", i, "_", chr_num_ext, "_", chr_num_int, "-pred_ampl.rds")
      #saveRDS(object = mean_preds_ampl, file = pred_ampl_path)
      
      #plot_path <- paste0("/home/ieo7429/Scrivania/results_regressor_gab/predictions/CSP/", i, "_", chr_num_ext, "_", chr_num_int, "-plot_ampl.png")
      #ggsave(filename = plot_path, plot = p, width = 7, height = 5)
    }
  }
}