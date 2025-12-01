# Required libraries
library(tidyr)
library(ggplot2)

# 1. Prepare the data
shap_df <- models.shap.df[[model]]
shap_df <- shap_df[,colnames(shap_df) != "BIAS"]

X_train <- models.X.train[[model]]
X_train$labels <- paste(X_train$bin,X_train$Type,sep = "-")
X_train <- X_train[,!colnames(X_train) %in% c("bin","Type","ampl_score","del_score")]

X_test <- models.X.test[[model]]
X_test$labels <- paste(X_test$bin,X_test$Type,sep = "-")
X_test <- X_test[,!colnames(X_test) %in% c("bin","Type","ampl_score","del_score")]

X <- rbind(X_train, X_test)

# Assuming X_test is the original feature data frame
# shap_long <- shap_df %>%
#   pivot_longer(
#     cols = -labels,          # Exclude the label column from pivoting
#     names_to = "Feature",    # Name for the reshaped column of feature names
#     values_to = "SHAP_value" # Name for the reshaped column of SHAP values
#   )
shap_long <- shap_df
colnames(shap_long) <- c('labels','Feature','SHAP_value')


# Add feature values for coloring (optional, but improves plot understanding)
feature_values <- as.data.frame(X) %>%
  pivot_longer(
    cols = -labels,
    names_to = "Feature",
    values_to = "Feature_value"
  )

# Merge SHAP values with feature values
shap_long <- cbind(shap_long, Feature_value = feature_values$Feature_value)


# Sort by importance
shap_importance <- models.shap.values[[model]]
shap_importance <- shap_importance[shap_importance$Feature != "BIAS",]
shap_importance <- shap_importance[shap_importance$Mean_Abs_SHAP > 0,]

# In the test data, some features' real values are NA, remove them
shap_long_rmna <- shap_long[complete.cases(shap_long),]
shap_long_rmna <- shap_long_rmna[shap_long_rmna$Feature %in% shap_importance$Feature,]

for(fea in as.character(unique(shap_long_rmna$Feature))){
  m <- shap_long_rmna[shap_long_rmna$Feature == fea,]
  pearson_corr <- cor.test(m$SHAP_value, m$Feature_value, method = "pearson")
  spearman_corr <- cor.test(m$SHAP_value, m$Feature_value, method = "spearman", exact = F) ### edit here added exact = F
  cor.plot <- ggplot(m, aes(x = SHAP_value, y = Feature_value)) +
    geom_point(color = background,size = 0.3) +
    geom_smooth(method = "lm", color = lines, se = FALSE, linewidth = 0.80) +
    annotate("text",
             x = min(m[,"SHAP_value"]),
             y = max(m[,"Feature_value"]),
             label = paste0(paste0("Pearson's R: ", round(pearson_corr$estimate, 3)),"\n",
                            paste0("Pearson's p-value: ", round(pearson_corr$p.value, 4)),"\n",
                            paste0("Spearman's rho: ", round(spearman_corr$estimate, 3)), "\n",
                            paste0("Spearman's p-value: ", round(spearman_corr$p.value, 4))),
             hjust = 0, vjust = 1, size = 4, color = "black") +
    labs(title = model, subtitle = fea) +
    theme_classic() +
    theme(
      text = element_text(size = 12, family = "Helvetica"),
      plot.title = element_text(size = 13, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(
        color = "black",  # Border color
        fill = NA,        # Transparent fill inside the border
        linewidth = 0.75          # Border thickness
      )
    )
  direction.plots[[model]][[fea]] <- cor.plot
  direction.m <- rbind(direction.m,c(model, fea, pearson_corr$estimate, pearson_corr$p.value, spearman_corr$estimate, spearman_corr$p.value))
}
