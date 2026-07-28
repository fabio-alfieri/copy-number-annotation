df <- read.table("results_regressor_gab/avg_performances_complete_chromosome_specific.tsv", header = T, sep = "\t", row.names = 1)

df_chr <- df %>%
  tibble::rownames_to_column("feature") %>%
  dplyr::filter(stringr::str_detect(feature, "Length_Mid-length_ampl")) %>%
  dplyr::mutate(
    chromosome = stringr::str_extract(feature, "(?<=CS_)\\d+"),
    chromosome = factor(chromosome, levels = as.character(1:22))
  ) %>%
  dplyr::select(chromosome, avg_r2) %>%
  dplyr::filter(!is.na(avg_r2))

ggplot(df_chr, aes(x = chromosome, y = avg_r2)) +
  geom_col(fill = "steelblue") +
  theme_bw() +
  labs(
    x = "Chromosome",
    y = expression(R^2),
    title = "Mid-length Amplifications"
  )
