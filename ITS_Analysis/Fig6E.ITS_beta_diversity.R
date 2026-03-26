## -------------------------------------------------------------------------
## Figure 6E: ITS beta-diversity plots 
##
## Purpose:
##   Load the cleaned phyloseq object 
##   and generate day-stratified beta-diversity (Bray-Curtis) PCoA plots
##   and PERMANOVA
##   summaries for treatment-group comparisons.
## -------------------------------------------------------------------------

library(tidyverse)
library(phyloseq)
library(vegan)
library(ggalt)

physeq <- readRDS("./data/ps.ITS.RDS")

group_levels <- c(
    "NoAbxNoCA", "NoCA", "NoSCFA",
    "AceticAcid", "PropionicAcid", "ButyricAcid", "MixedAcid"
)

group_colors = c("#7D8EA1", "#45B06A", "#0152A1", "#84D2F6", "#9D1B28", "#FFA733", "#FFCE00")


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
    
    color_map <- stats::setNames(group_colors, group_levels)
    
    ps.dist <- phyloseq::distance(physeq_temp, "bray")
    pcoa.ape <- ape::pcoa(ps.dist, correction = "cailliez")
    
    pcoa_plot_data <- data.frame(
        pco1 = pcoa.ape[["vectors.cor"]][, 1],
        pco2 = pcoa.ape[["vectors.cor"]][, 2]
    )
    proportion_of_variance <- pcoa.ape[["values"]][["Rel_corr_eig"]]
    
    df <- data.frame(
        pco1 = pcoa_plot_data$pco1,
        pco2 = pcoa_plot_data$pco2,
        Group = phyloseq::sample_data(physeq_temp)[["Group"]]
    )
    
    xlabel <- sprintf("PCo1 (Proportion of Variance: %0.2f%%)", proportion_of_variance[1] * 100)
    ylabel <- sprintf("PCo2 (Proportion of Variance: %0.2f%%)", proportion_of_variance[2] * 100)
    
    metadata <- data.frame(phyloseq::sample_data(physeq_temp))
    category <- metadata[["Group"]]
    
    disper <- vegan::betadisper(ps.dist, category)
    print(disper)
    
    permutest.result <- vegan::permutest(disper, permutations = how(nperm = 9999))
    permutest.result
    
    if (permutest.result$tab$`Pr(>F)`[1] < 0.05){
        warning(
                "    Dispersion of data is significantly different among groups. 
    Interpret PERMANOVA results with caution. ")
    }
    
    formula_text <- as.formula("ps.dist ~ Group")
    permanova.result <- vegan::adonis2(
        formula = formula_text,
        data = metadata,
        permutations = 9999,
        method = "bray"
    )
    
    g <- ggplot2::ggplot(df, ggplot2::aes(x = pco1, y = pco2)) +
        ggplot2::stat_ellipse(
            geom = "polygon",
            ggplot2::aes(fill = Group),
            alpha = 0.2,
            linetype = 2,
            show.legend = FALSE,
            level = 0.95
        ) +
        ggplot2::geom_point(ggplot2::aes(fill = Group), colour = "black", pch = 21, size = 5) +
        ggplot2::scale_fill_manual(values = group_colors) +
        ggplot2::theme_bw() +
        ggplot2::xlab(xlabel) +
        ggplot2::ylab(ylabel) 
    
    plot_filename <- sprintf("beta_diversity_bray.pcoa.%s.pdf", day)
    ggsave(g, filename = plot_filename, width = 8, height = 7, device = "pdf", dpi = 300)
    
}
