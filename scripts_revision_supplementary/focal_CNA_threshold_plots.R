# ============================================================
# Show that 3 Mb is stricter than GISTIC 50% chromosome-arm cutoff
# ============================================================

library(data.table)
library(ggplot2)
library(rtracklayer)

genome_build <- "hg19"
small_scale_cutoff_bp <- 3e6

# -------------------------
# Download cytoband table from UCSC
# -------------------------

session <- browserSession("UCSC")
genome(session) <- genome_build

query <- ucscTableQuery(session, table = "cytoBandIdeo")
cyto <- as.data.table(getTable(query))

cyto <- cyto[chrom %in% paste0("chr", c(1:22, "X"))]

# -------------------------
# Compute chromosome arm lengths
# -------------------------

chrom_sizes <- cyto[
  ,
  .(chrom_end = max(chromEnd)),
  by = chrom
]

centromeres <- cyto[
  gieStain == "acen",
  .(
    centromere_start = min(chromStart),
    centromere_end = max(chromEnd)
  ),
  by = chrom
]

arms <- merge(chrom_sizes, centromeres, by = "chrom")

arm_lengths <- rbindlist(list(
  arms[
    ,
    .(
      chrom,
      arm = "p",
      arm_start = 0,
      arm_end = centromere_start,
      arm_length_bp = centromere_start
    )
  ],
  arms[
    ,
    .(
      chrom,
      arm = "q",
      arm_start = centromere_end,
      arm_end = chrom_end,
      arm_length_bp = chrom_end - centromere_end
    )
  ]
))

arm_lengths[
  ,
  chromosome_arm := paste0(chrom, arm)
]

arm_lengths[
  ,
  gistic_50pct_arm_bp := 0.5 * arm_length_bp
]

arm_lengths[
  ,
  small_scale_cutoff_bp := small_scale_cutoff_bp
]

arm_lengths[
  ,
  gistic_50pct_arm_Mb := gistic_50pct_arm_bp / 1e6
]

arm_lengths[
  ,
  small_scale_cutoff_Mb := small_scale_cutoff_bp / 1e6
]

arm_lengths[
  ,
  cutoff_as_pct_of_arm := 100 * small_scale_cutoff_bp / arm_length_bp
]

arm_lengths[
  ,
  fold_stricter_than_gistic := gistic_50pct_arm_bp / small_scale_cutoff_bp
]

arm_lengths[
  ,
  threeMb_is_below_gistic_threshold := small_scale_cutoff_bp < gistic_50pct_arm_bp
]

arm_lengths <- arm_lengths[order(gistic_50pct_arm_Mb)]

arm_lengths[
  ,
  gistic_label := paste0(round(gistic_50pct_arm_Mb, 1), " Mb")
]

# Plots
(p1 <- ggplot(
  arm_lengths,
  aes(
    x = reorder(chromosome_arm, gistic_50pct_arm_Mb),
    y = gistic_50pct_arm_Mb
  )
) +
  geom_col(fill = 'lightblue') +
  geom_text(
    aes(label = gistic_label),
    hjust = -0.15,
    size = 3
  ) +
  geom_hline(
    yintercept = 3,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18))
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.margin = margin(5.5, 40, 5.5, 5.5)
  ) +
  labs(
    title = "Focal CNA definition: GISTIC 50% arm vs. 3Mbp",
    subtitle = paste0("Dashed line indicates 3 Mb."),
    x = "Chromosome arm",
    y = "50% of chromosome arm length (Mb)"
  ))

ggsave(
  "../plot_3Mb_vs_GISTIC_50pct_chromosome_arm_threshold.pdf",
  p1,
  width = 7,
  height = 7
)

library(patchwork)

# ------------------------------------------------------------
# Common chromosome-arm order
# ------------------------------------------------------------

arm_order <- arm_lengths[
  order(gistic_50pct_arm_Mb)
]$chromosome_arm

arm_lengths[
  ,
  chromosome_arm := factor(chromosome_arm, levels = arm_order)
]

arm_lengths[
  ,
  gistic_label := paste0(round(gistic_50pct_arm_Mb, 1), " Mb")
]

arm_lengths[
  ,
  pct_label := paste0(round(cutoff_as_pct_of_arm, 1), "%")
]

p1 <- ggplot(
  arm_lengths,
  aes(x = chromosome_arm, y = gistic_50pct_arm_Mb)
) +
  geom_col(fill = 'lightblue') +
  geom_text(
    aes(label = gistic_label),
    hjust = -0.15,
    size = 3
  ) +
  geom_hline(
    yintercept = 3,
    linetype = "dashed",
    linewidth = 0.8,
    col = 'red'
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.22))
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  ) +
  labs(
    title = "GISTIC <50% chromosome-arm",
    subtitle = "Dashed line = 3 Mb cutoff",
    x = "Chromosome arm",
    y = "50% of chromosome arm length (Mb)"
  )


p2 <- ggplot(
  arm_lengths,
  aes(x = chromosome_arm, y = cutoff_as_pct_of_arm)
) +
  geom_col(fill = 'lightblue') +
  geom_text(
    aes(label = pct_label),
    hjust = -0.15,
    size = 3
  ) +
  geom_hline(
    yintercept = 50,
    linetype = "dashed",
    linewidth = 0.8,
    col = 'red'
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.22))
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  ) +
  labs(
    title = "3 Mb as fraction of chromosome arm",
    subtitle = "Dashed line = 50% arm boundary",
    x = NULL,
    y = "3 Mb as % of chromosome arm"
  )


combined_plot <- p1 + p2 +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title = "A 3 Mb cutoff is substantially more stringent than the GISTIC definition",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )
combined_plot

ggsave(
  "../combined_3Mb_vs_GISTIC_arm_threshold.pdf",
  combined_plot,
  width = 14,
  height = 8
)
