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
#usethis::use_package("AnnotationHub")
usethis::use_package("hdf5r")
usethis::use_package("locuszoomr")
usethis::use_package("GenomicRanges")
usethis::use_package("rtracklayer")
usethis::use_package("leafviz")
usethis::use_package("shinycssloaders")


###############################################################################
# fig colors
plot_colors <- 
    read_tsv("../../bcellactivation/paper_plots/figure_colors.txt", 
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

###############################################################################
# Bulk RNA-seq
dat <- read_rds("../../bcellactivation/bcell_lowinput/data/deseq_normalized_counts.rds")

bulk_genes <- unique(dat$gene_label)

gene_exp <- 
    dat |>
    pivot_wider(names_from = gene_label, values_from = norm_counts) |>
    unite("condition", c(stim, timep), sep = "_", remove = FALSE) |>
    mutate_at(vars(stim, timep), ~fct_inorder(as.character(.)))


###############################################################################
# Single-cell data
sc_data <- 
    SeuratDisk::LoadH5Seurat("../../bcellactivation//citeseq/data/bcells.h5Seurat")

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

sc_meta <- 
    sc_data_sub@meta.data |>
    tibble::as_tibble(rownames = "barcode") |>
    dplyr::select(barcode, hto = dmm_hto_call, cluster = seurat_clusters) |>
    dplyr::mutate(cluster = factor(cluster, levels = paste0("C", 0:13)))

umap_df <- 
    Embeddings(sc_data_sub, "umap") |>
    tibble::as_tibble(rownames = "barcode") |>
    dplyr::left_join(sc_meta, dplyr::join_by(barcode)) |>
    dplyr::mutate(hto = recode(hto, 
			"Unstim 0h" = "Unstim_0",
			"IL4 24h" = "IL-4c_24",
			"IL4 72h" = "IL-4c_72",
			"TLR7 24h" = "TLR7c_24",
			"TLR7 72h" = "TLR7c_72",
			"BCR 24h" = "BCRc_24",
			"BCR 72h" = "BCRc_72",
			"DN2 24h" = "DN2c_24",
			"DN2 72h" = "DN2c_72"))

cluster_labels <-
    umap_df |>
    dplyr::group_by(cluster) |>
    dplyr::summarise_at(dplyr::vars(UMAP_1, UMAP_2), mean) |>
    dplyr::ungroup()

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

sc_var_plots <- list("hto" = sc_hto_plot, "cluster" = sc_clusters_plot)


###############################################################################
# ATAC-seq

# Use a new R version just to run AnnotationHub, otherwise keep R 4.1
#ah <- AnnotationHub::AnnotationHub()
#
#AnnotationHub::query(ah, c("EnsDb", "Homo sapiens"))
#
#ens_data <- ah[["AH119325"]]
#
#db_file <- dbconn(ens_data)@dbname 
#file.copy(db_file, "../inst/extdata/Homo_sapiens.GRCh38.ensdb.sqlite", overwrite = TRUE)
#

ens_data <- ensembldb::EnsDb("../inst/extdata/Homo_sapiens.GRCh38.ensdb.sqlite")
atac_genes <- ensembldb::keys(ens_data, keytype = "SYMBOL")
atac_genes <- atac_genes[! atac_genes == ""]

# Download the MANE Select list
mane_data <- 
    "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz" |>
    readr::read_tsv()

mane_tx <- 
    mane_data |>
    dplyr::select(symbol, MANE_status, Ensembl_nuc) |>
    dplyr::mutate_at(vars(Ensembl_nuc), ~str_remove(., "\\.\\d+$")) |>
    pull(Ensembl_nuc)


gtf_file <- "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v49.primary_assembly.annotation.gtf.gz"

gtf <- vroom::vroom(gtf_file, comment = "#", col_names = FALSE)

gtf_tx <- 
    dplyr::filter(gtf, X3 == "transcript") |>
    tibble::rowid_to_column() |>
    dplyr::mutate(gene_id = stringr::str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
		  transcript_id = stringr::str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
		  is_mane = grepl('tag "MANE_Select"|tag "MANE_Plus_Clinical"', X9)) |>
    dplyr::select(rowid, gene_id, transcript_id, is_mane)

# Build the blacklist
exc_tx <- 
    gtf_tx |>
    dplyr::group_by(gene_id) |>
    dplyr::filter(any(is_mane == TRUE)) |>
    dplyr::ungroup() |>
    dplyr::filter(is_mane == FALSE) |>
    dplyr::pull(transcript_id)

# Free up memory
rm(gtf, gtf_tx)
gc()

# Load full GTF via rtracklayer
gtf_rtrack <- rtracklayer::import(gtf_file)

# THE FIXES FOR LOCUSZOOMR
# Force Ensembl chromosome names (chr1 -> 1) to match GWAS data
GenomeInfoDb::seqlevelsStyle(gtf_rtrack) <- "Ensembl"

# Map GENCODE biotypes to Ensembl biotypes so gg_genetracks can color them
gtf_rtrack$gene_biotype <- gtf_rtrack$gene_type
gtf_rtrack$transcript_biotype <- gtf_rtrack$transcript_type

# Drop blacklisted transcripts AND drop top-level 'gene' rows
# Dropping 'gene' rows forces ensembldb to recalculate boundaries, fixing the long line bug
exclude_mask <- (gtf_rtrack$transcript_id %in% exc_tx) | (gtf_rtrack$type == "gene")
gtf_filtered <- gtf_rtrack[!exclude_mask]

# Export and Build
rtracklayer::export(gtf_filtered, "gencode_v49_filtered.gtf", format = "gtf")

ensembldb::ensDbFromGtf(
    gtf = "gencode_v49_filtered.gtf",
    outfile = "GencodeDb.sqlite",
    organism = "Homo_sapiens",
    genomeVersion = "GRCh38",
    version = "49"
)

file.copy("GencodeDb.sqlite", "../inst/extdata/Homo_sapiens.GRCh38.ensdb2.sqlite", overwrite = TRUE)

unlink("gencode_v49_filtered.gtf")
unlink("GencodeDb.sqlite")


###############################################################################
# Splicing
splicing_contrasts <- 
    c("Unstim 0h vs. TLR7c 24h" = "unstday0vs.TLR7",
      "Unstim 0h vs. BCRc 24h" = "unstday0vs.BCR",
      "Unstim 0h vs. DN2c 24h" = "unstday0vs.DN2",
      "TLR7c 24h vs. BCRc 24h" = "TLR7vs.BCR",
      "TLR7c 24h vs. DN2c 24h" = "TLR7vs.DN2",
      "BCRc 24h vs. DN2c 24h" = "BCRvs.DN2")

usethis::use_data(plot_colors, 
		  atac_colors, 
		  cluster_colors, 
		  gene_exp,
		  bulk_genes,
		  sc_genes,
		  sc_cells,
		  sc_meta,
		  umap_df,
		  cluster_labels,
		  sc_clusters_plot,
		  sc_hto_plot,
		  sc_var_plots,
		  atac_genes,
		  mane_tx,
		  splicing_contrasts,
		  internal = TRUE, overwrite = TRUE)
