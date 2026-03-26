## -------------------------------------------------------------------------
## Figure 6B-D: ITS taxa composition barplots
##
## Purpose:
##   Load the cleaned phyloseq object
##   and generate stacked barplots for genus-level relative abundance.
##   
## -------------------------------------------------------------------------

library(tidyverse)
library(phyloseq)
library(RColorBrewer)

physeq <- readRDS("./data/ps.ITS.RDS")

group_levels <- c(
    "NoAbxNoCA", "NoCA", "NoSCFA",
    "AceticAcid", "PropionicAcid", "ButyricAcid", "MixedAcid"
)

day_list <- phyloseq::sample_data(physeq) %>%
    data.frame() %>%
    dplyr::pull(Date) %>%
    unique()

for (day in day_list) {
    
    physeq_temp <- phyloseq::subset_samples(physeq, Date == day)
    physeq_temp <- phyloseq::subset_samples(physeq_temp, !is.na(phyloseq::sample_data(physeq_temp)[["Group"]]))
    physeq_temp <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq_temp) > 0, physeq_temp)
    
    phyloseq::sample_data(physeq_temp)[["Group"]] <- factor(
        phyloseq::sample_data(physeq_temp)[["Group"]],
        levels = group_levels
    )
    
    unclassified_label <- "Unclassified"
    other_label <- "Other (low abundance)"
    
    physeq_rank <- phyloseq::tax_glom(physeq_temp, taxrank = "Genus", NArm = FALSE)
    physeq_rel <- phyloseq::transform_sample_counts(physeq_rank, function(x) x / sum(x))
    
    df <- phyloseq::psmelt(physeq_rel) %>%
        dplyr::mutate(
            raw_tax = as.character(.data[["Genus"]]),
            Taxa = dplyr::case_when(
                is.na(raw_tax) | raw_tax == "" | raw_tax == "NA" ~ unclassified_label,
                grepl("^[A-Za-z]__$", raw_tax) ~ unclassified_label,
                TRUE ~ raw_tax
            )
        )
    
    # keep the top 10 taxa
    top_taxa <- df %>%
        dplyr::group_by(Taxa) %>%
        dplyr::summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
        dplyr::slice_max(mean_abund, n = 10) %>%
        dplyr::pull(Taxa)
    
    df2 <- df %>%
        dplyr::mutate(
            Taxa = dplyr::case_when(
                Taxa == unclassified_label ~ unclassified_label,
                Taxa %in% top_taxa ~ Taxa,
                TRUE ~ other_label
            )
        ) %>%
        dplyr::group_by(Sample, Group, Taxa) %>%
        dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")
    
    taxa_levels <- unique(as.character(df2$Taxa))
    normal_taxa <- setdiff(taxa_levels, c(unclassified_label, other_label))
    
    # set taxa color
    pal_name <- "Spectral"
    max_n <- RColorBrewer::brewer.pal.info[pal_name, "maxcolors"]
    base_pal <- RColorBrewer::brewer.pal(min(max_n, max(3, length(normal_taxa))), pal_name)
    
    pal_vals <- if (length(normal_taxa) <= length(base_pal)) {
        base_pal[seq_along(normal_taxa)]
    } else {
        grDevices::colorRampPalette(base_pal)(length(normal_taxa))
    }
    
    normal_cols <- setNames(rev(pal_vals), normal_taxa)
    special_cols <- setNames(c("#7F7F7F", "#D9D9D9"), c(unclassified_label, other_label))
    taxa_color <- c(normal_cols, special_cols)
    
    df2$Taxa <- factor(df2$Taxa, levels = c(normal_taxa, unclassified_label, other_label))
    
    g <- ggplot2::ggplot(df2, ggplot2::aes(x = Sample, y = Abundance, fill = Taxa)) +
        ggplot2::geom_bar(position = "stack", stat = "identity") +
        ggplot2::scale_fill_manual(values = taxa_color, drop = FALSE) +
        facet_wrap(~Group, scales = "free_x", ncol = length(group_levels) ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.9, hjust = 0.7),
            axis.ticks.x = ggplot2::element_blank(),
            panel.spacing = grid::unit(0, "lines"),
            strip.background = ggplot2::element_rect(color = "grey30", fill = "grey90")
        ) +
        ggplot2::labs(x = "", y = "Relative abundance", fill = "Genus") +
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1, byrow = FALSE))
    g
    
    plot_filename <- sprintf("ITS_profile.taxa_Genus.stacked_barplot.%s.pdf", day)
    ggsave(plot = g, 
           filename = plot_filename, 
           width = 15,
           height = 7,
           dpi = 300)
}






