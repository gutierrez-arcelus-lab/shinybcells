library(tidyverse)
library(vroom)
library(Seurat)
library(hdf5r)

usethis::use_package("shiny")
usethis::use_package("pkgload")
usethis::use_package("dplyr")
usethis::use_package("tidyr")
usethis::use_package("readr")
usethis::use_package("stringr")
usethis::use_package("purrr")
usethis::use_package("tibble")
usethis::use_package("vroom")
usethis::use_package("ggplot2")
usethis::use_package("patchwork")
usethis::use_package("dplyr")
usethis::use_package("ggbeeswarm")
usethis::use_package("ensembldb")
usethis::use_package("hdf5r")
usethis::use_package("locuszoomr")
usethis::use_package("GenomicRanges")
usethis::use_package("rtracklayer")
usethis::use_package("leafviz")



###############################################################################
# fig colors
plot_colors <- 
    read_tsv("./figure_colors_v2.txt", 
	     col_names = c("stim", "timep", "col")) |>
    unite("lab", c(stim, timep), sep = "_") |>
    deframe()

atac_colors <- plot_colors[c("Unstim_0", "Unstim_24", "IL-4c_24", "TLR7c_24", "BCRc_24", "DN2c_24")]
atac_colors[c("TLR7c_24", "BCRc_24", "DN2c_24")] <- plot_colors[c("TLR7c_24", "BCRc_48", "DN2c_48")] 

cluster_colors <- 
    c("C0" = "#A8CDE2", "C1" = "#3B83B9", "C2" = "#E3362C", "C3" = "#F9B56F", 
      "C4" = "#FC9230", "C5" = "#DDA086", "C6" = "#9F7BB8", "C7" = "#987898", 
      "C8" = "#F1E78D", "C9" = "#B05D2F", "C10" = "#83BF98", "C11" = "#6ABD5D", 
      "C12" = "#6F8544", "C13" = "#F4817F")

usethis::use_data(plot_colors)
usethis::use_data(atac_colors)
usethis::use_data(cluster_colors)

###############################################################################
# Bulk RNA-seq
dat <- read_rds("./deseq_normalized_counts.rds")

bulk_genes <- unique(dat$gene_label)

usethis::use_data(bulk_genes)

gene_exp <- dat |>
    pivot_wider(names_from = gene_label, values_from = norm_counts) |>
    unite("condition", c(stim, timep), sep = "_", remove = FALSE) |>
    mutate_at(vars(stim, timep), ~fct_inorder(as.character(.)))

usethis::use_data(gene_exp)

###############################################################################
# Single-cell data
sc_data <- SeuratDisk::LoadH5Seurat("../../lupus/citeseq/data/bcells.h5Seurat")

genes_expressed <- 
    GetAssayData(object = sc_data, assay = "RNA", slot = "data") |>
    apply(1, function(x) sum(x > 0)) |>
    {function(x) which(x >= 10)}() |>
    names()

sc_data_sub <- subset(sc_data, features = genes_expressed)

expr_mat <- 
    as.matrix(GetAssayData(sc_data_sub, assay = "RNA", slot = "data"))

chunk_dims <- c(1, ncol(expr_mat)) 

h5file <- H5File$new("../inst/extdata/bcells_expressed.h5", mode = "w")

dset <- 
    h5file$create_dataset(
			  name = "expr_data",
			  dtype = h5types$H5T_NATIVE_FLOAT,  
			  dims = dim(expr_mat),
			  chunk_dims = chunk_dims
)

for (i in 1:nrow(expr_mat)) 
  dset[i, ] <- expr_mat[i, , drop = FALSE]

h5file$close_all()

sc_genes <- rownames(expr_mat)
sc_cells <- colnames(expr_mat)

usethis::use_data(sc_genes)
usethis::use_data(sc_cells)

sc_meta <- 
    sc_data_sub@meta.data |>
    as_tibble(rownames = "barcode") |>
    select(barcode, hto = dmm_hto_call, cluster = seurat_clusters) |>
    mutate(cluster = factor(cluster, levels = paste0("C", 0:13)))

usethis::use_data(sc_meta, overwrite = TRUE)

umap_df <- 
    Embeddings(sc_data_sub, "umap") |>
    as_tibble(rownames = "barcode") |>
    left_join(sc_meta, join_by(barcode)) |>
    mutate(hto = recode(hto, 
			"Unstim 0h" = "Unstim_0",
			"IL4 24h" = "IL-4c_24",
			"IL4 72h" = "IL-4c_72",
			"TLR7 24h" = "TLR7c_24",
			"TLR7 72h" = "TLR7c_72",
			"BCR 24h" = "BCRc_24",
			"BCR 72h" = "BCRc_72",
			"DN2 24h" = "DN2c_24",
			"DN2 72h" = "DN2c_72"))

usethis::use_data(umap_df, overwrite = TRUE)

cluster_labels <-
    umap_df |>
    group_by(cluster) |>
    summarise_at(vars(UMAP_1, UMAP_2), mean) |>
    ungroup()

sc_clusters_plot <- 
    ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(fill = cluster), 
               size = 1.5, 
               shape = 21, 
               stroke = .05, 
               color = "black") +
    geom_label(data = cluster_labels, 
               aes(x = UMAP_1, y = UMAP_2, label = cluster),
               label.padding = unit(0.1, "lines"),
               size = 12, size.unit = "pt", alpha = .5, fontface = "bold") +
    scale_fill_manual(values = cluster_colors) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          panel.grid = element_blank()) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(fill = guide_legend(title = "Stim:",
                               title.position = "top",
                               override.aes = list(size = 3)))
                
sc_hto_plot <- 
    ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(fill = hto), 
               size = 1.5, 
               shape = 21, 
               stroke = .05, 
               color = "black") +
    scale_fill_manual(values = plot_colors) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          panel.grid = element_blank()) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(fill = guide_legend(title = "Stim:",
                               title.position = "top",
                               override.aes = list(size = 3))) 

sc_var_plots <- list("HTO" = sc_hto_plot, "Clusters" = sc_clusters_plot)

usethis::use_data(cluster_labels)
usethis::use_data(sc_clusters_plot)
usethis::use_data(sc_hto_plot)
usethis::use_data(sc_var_plots)

###############################################################################
# ATAC-seq
ah <- AnnotationHub::AnnotationHub()
ens_data <- ah[["AH98047"]]

db_file <- dbconn(ens_data)@dbname 
file.copy(db_file, "../inst/extdata/Homo_sapiens.GRCh38.ensdb.sqlite")

atac_genes <- keys(ens_data, keytype = "SYMBOL")
atac_genes <- atac_genes[! atac_genes == ""]

usethis::use_data(atac_genes)


###############################################################################
# Splicing
splicing_contrasts <- 
    c("Unstim 0h vs. TLR7c 24h" = "unstday0vs.TLR7",
      "Unstim 0h vs. BCRc 24h" = "unstday0vs.BCR",
      "Unstim 0h vs. DN2c 24h" = "unstday0vs.DN2",
      "TLR7c 24h vs. BCRc 24h" = "TLR7vs.BCR",
      "TLR7c 24h vs. DN2c 24h" = "TLR7vs.DN2",
      "BCRc 24h vs. DN2c 24h" = "BCRvs.DN2")

usethis::use_data(splicing_contrasts, overwrite = TRUE)

