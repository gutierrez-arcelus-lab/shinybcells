format_timecourse <- function(gene_i) {
    
    d <- 
        gene_exp |> 
        dplyr::select(1:4, norm_counts = {{ gene_i }})
    
    d_tmp <- dplyr::filter(d, condition != "Unstim_0")
    
    dplyr::filter(d, condition == "Unstim_0") |> 
        dplyr::select(-stim) |>
        tidyr::expand_grid(dplyr::distinct(d_tmp, stim)) |>
        dplyr::bind_rows(d_tmp)
}  

get_sc_data <- function(gene_i) {
    
    require(hdf5r)
    
    idx <- which(sc_genes == gene_i)
    h5file <- H5File$new(system.file("extdata", "bcells_expressed.h5", package = "shinybcells"))
    h5data <- h5file[["expr_data"]]
    gene_data <- h5data$read(args = list(idx, quote(expr=)))
    gene_df <- tibble::tibble( "{gene_i}" := gene_data)
    h5file$close_all()
    return(gene_df)
}

get_ensdb <- function() {
    dbfile <- system.file("extdata", "Homo_sapiens.GRCh38.ensdb.sqlite", package = "shinybcells")
    if (!nzchar(dbfile)) stop("EnsDb file not found.")
    EnsDb(dbfile)
}