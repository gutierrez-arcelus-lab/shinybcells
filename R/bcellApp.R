library(shiny)
library(ggplot2)
library(patchwork)
library(ensembldb)


bcellApp <- function(...) {
    ui <- navbarPage(
        "B cell activation",
        tabPanel("Time course RNA-seq",
                 selectizeInput("generna", "Select gene:", choices = NULL),
                 plotOutput("plotrna", width = "100%", height = "200px")
        ),
        tabPanel("Bulk ATAC-seq",
                 radioButtons("atac_plot_mode", "Plot by:", choices = c("Gene" = "gene", "Genomic Position" = "position"),
                              selected = "gene"),
                 conditionalPanel(
                     condition = "input.atac_plot_mode == 'gene'",
                     selectizeInput("atac_gene", "Select gene:", choices = NULL)
                 ),
                 conditionalPanel(
                     condition = "input.atac_plot_mode == 'position'",
                     textInput("atac_coords", "Enter GRCh38 Genomic Coordinates (e.g. chr1:10000)")
                 ),
                 numericInput("atac_window", "Window Size", value = 10e3, min = 1e3, max = 500e3),
                 actionButton("makeatacplot", "Click to plot!"),
                 plotOutput("plotatac", width = "100%", height = "600px")
        ),
        tabPanel("Single-cell RNA-seq", 
                 fluidRow(
                     column(6, selectizeInput("varsc", "Select variable:", choices = c("HTO" = "hto", "Clusters" = "cluster"))),
                     column(6, selectizeInput("genesc", "Select gene:", choices = NULL))
                 ),
                 fluidRow(
                     column(6, plotOutput("plotscvars", width = "100%", height = "600px")),
                     column(6, plotOutput("plotsc", width = "100%", height = "600px"))
                 ),
                 titlePanel("Marker genes"),
                 fluidRow(
                     column(4, selectizeInput("markergenes", "Select genes:", multiple = TRUE, choices = NULL)),
                     column(8, plotOutput("bubbleplot", inline = TRUE))
                 )
        ),
        tabPanel("Splicing",
                 HTML('<p style="text-align:right;"><a href="https://github.com/jackhump/leafviz" target="_blank">Powered by LeafViz</a></p>'),
                 selectizeInput("contrast", "Select contrast:", choices = names(splicing_contrasts)),
                 fluidRow(
                     column(6,
                            div(id = "clusterTable",
                                h4(id = "title","Differential splicing events (clusters)"),
                                hr(),
                                div(DT::DTOutput("all_clusters"))
                            )
                     ),
                     column(6,
                            div(id="clusterView",
                                h4(id="title","Splicing event visualization"),
                                hr(),
                                div(plotOutput("select_cluster_plot", width = "100%")),
                                DT::DTOutput("cluster_view") 
                            )
                     )
                 )
        )
    )
    
    server <- function(input, output, session) {
    
        #ensdb <- get_ensdb()
        ah <- AnnotationHub::AnnotationHub()
        ensdb <- ah[["AH98047"]]
        
        bigwigs <- 
            list.files(system.file("extdata", package = "shinybcells"),
                       pattern = "*_filtered.bigWig")
        
        names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", bigwigs)
        
        ## Bulk RNA
        updateSelectizeInput(session, 
                             "generna", 
                             selected = "CD19",
                             choices = bulk_genes, 
                             server = TRUE)
        
        datarna <- reactive({
            req(input$generna)
            format_timecourse(input$generna)
        })
        
        output$plotrna <- 
            renderPlot(res = 96, {
                req(input$generna)
                ggplot(data = datarna()) +
                    ggbeeswarm::geom_quasirandom(aes(x = timep, y = norm_counts, fill = condition),
                                                 size = 2.5, stroke = .25, shape = 21,
                                                 method = "smiley", width = .2) +
                    scale_fill_manual(values = plot_colors) +
                    facet_grid(cols = vars(stim), space = "free", scale = "free_x") +
                    theme_bw() +
                    guides(fill = "none") +
                    labs(x = "hours", y = "Norm. counts")
            })
        
        # Bulk ATAC
        updateSelectizeInput(session,
                             "atac_gene",
                             selected = "CD19",
                             choices = atac_genes,
                             server = TRUE)
        
        atac_plot <-
            eventReactive(input$makeatacplot, {
                
                if (input$atac_plot_mode == "gene") {
                    
                    loc <-
                        locuszoomr::locus(gene = input$atac_gene,
                                          flank = input$atac_window,
                                          ens_db = ensdb)
                    
                } else if (input$atac_plot_mode == "position") {
                    
                    
                    coords_split <- 
                        stringr::str_split(input$atac_coords, ":") |>
                        unlist()
                    
                    atac_chrom <- coords_split[[1]]
                    atac_pos <- readr::parse_number(coords_split[[2]])
                    
                    region <- c(atac_pos - input$atac_window, 
                                atac_pos + input$atac_window)
                    
                    loc <-
                        locuszoomr::locus(seqname = atac_chrom,
                                          xrange = region,
                                          ens_db = ensdb)
                }
                
                atac_peaks <-
                    interv <-
                    GenomicRanges::GRanges(paste0("chr", loc$seqname),
                                           IRanges(loc$xrange[1], loc$xrange[2]))
                
                atac_ranges <- 
                    system.file("extdata", bigwigs, package = "shinybcells") |>
                    setNames(names(bigwigs)) |>
                    purrr::map(~rtracklayer::import(., which = interv))
                
                atac_covered <-
                    atac_ranges |>
                    purrr::map_dfr(as.data.frame, .id = "stim") |>
                    dplyr::select(stim, start, end, score)
                
                atac_gaps <-
                    atac_ranges |>
                    purrr::map_dfr(~ranges(.) |> 
                                       keepSeqlevels(loc$seqname) |> 
                                       gaps(start = loc$xrange[1], end = loc$xrange[2]) |> 
                                       as.data.frame(),
                                   .id = "stim") |>
                    dplyr::mutate(score = 0) |>
                    dplyr::select(stim, start, end, score)
                
                atac_peaks <-
                    dplyr::bind_rows(atac_covered, atac_gaps) |>
                    dplyr::mutate(stim = stringr::str_replace(stim, "unst", "Unstim"),
                                  stim = stringr::str_replace(stim, "IL4", "IL-4c"),
                                  stim = stringr::str_replace(stim, "TLR7", "TLR7c"),
                                  stim = stringr::str_replace(stim, "BCR", "BCRc"),
                                  stim = stringr::str_replace(stim, "DN2", "DN2c"),
                                  stim = factor(stim, levels = names(atac_colors))) |>
                    dplyr::arrange(stim, start) |>
                    tidyr::pivot_longer(start:end, names_to = "dummy", values_to = "pos")
                
                plot_atac <-
                    ggplot(atac_peaks) +
                    geom_ribbon(aes(x = pos, ymin = 0, ymax = score, color = stim, fill = stim),
                                linewidth = .5, outline.type = "full", alpha = .5) +
                    scale_x_continuous(limits = loc$xrange,
                                       labels = function(x) round(x/1e6L, 2),
                                       expand = c(0, 0)) +
                    scale_color_manual(values = atac_colors) +
                    scale_fill_manual(values = atac_colors) +
                    facet_wrap(~stim, ncol = 1, strip.position = "right") +
                    theme_minimal() +
                    theme(
                        axis.text = element_text(size = 12),
                        legend.position = "none",
                        strip.text.y.right = element_text(angle = 0, size = 12),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor.y = element_blank(),
                        plot.background = element_rect(color = "white", fill = "white")) +
                    labs(x = NULL)
                
                gene_tracks <-
                    locuszoomr::gg_genetracks(loc, cex.text = 1) +
                    scale_x_continuous(limits = loc$xrange/1e6,
                                       labels = function(x) round(x, 2),
                                       expand = c(0, 0)) +
                    theme_minimal() +
                    theme(
                        axis.text = element_text(size = 12),
                        axis.title = element_text(size = 12),
                        plot.background = element_rect(color = "white", fill = "white"))
                
                plot_atac / gene_tracks + plot_layout(heights = c(1, .3))
            })
        
        output$plotatac <-
            renderPlot({
                input$makeatacplot
                atac_plot()
            })
        
        #Single-cell RNA-seq
        output$plotscvars <- renderPlot(sc_var_plots[[input$varsc]])
        
        updateSelectizeInput(session, 
                             "genesc", 
                             selected = "CD19",
                             choices = sc_genes, 
                             server = TRUE)
        
        sc_data_gene <- 
            reactive({
                req(input$genesc)
                get_sc_data(input$genesc) |>
                    tibble::add_column(barcode = sc_cells, .before = 1) |>
                    dplyr::select(barcode, value = 2) |>
                    dplyr::right_join(umap_df, dplyr::join_by(barcode)) |>
                    dplyr::arrange(value) |> 
                    dplyr::mutate(barcode = forcats::fct_inorder(barcode))
            })
        
        output$plotsc <- 
            renderPlot(
                ggplot(data = sc_data_gene(), 
                       aes(x = UMAP_1, y = UMAP_2)) +
                    geom_point(aes(color = value), size = 1.5, stroke = 0) +
                    scale_color_gradientn(colors = c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84",
                                                     "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000")) +
                    theme_minimal() +
                    theme(axis.text = element_text(size = 12),
                          axis.title = element_text(size = 12),
                          legend.text = element_text(size = 12),
                          legend.title = element_text(size = 12),
                          panel.grid = element_blank()) +
                    guides(color = guide_colorbar(barwidth = .5, barheight = 8))
            )
        
        ### Bubble plot
        updateSelectizeInput(session, 
                             "markergenes", 
                             choices = sc_genes, 
                             server = TRUE)
        
        markers_plot_data <-
            reactive({
                req(input$markergenes, input$varsc)
                
                n_genes <- length(input$markergenes)
                grouping_var <- input$varsc
                var <- rlang::sym(grouping_var)
                
                if (n_genes == 1) {
                    
                    dat <- 
                        get_sc_data(input$markergenes) |>
                        tibble::add_column(barcode = sc_cells, .before = 1) |>
                        dplyr::left_join(sc_meta, dplyr::join_by(barcode)) |>
                        tidyr::pivot_longer(-c(barcode, hto, cluster), names_to = "gene") |>
                        dplyr::group_by(!!var, gene) |>
                        dplyr::summarise(avg.exp = mean(expm1(value)),
                                         pct.exp = mean(value > 0) * 100) |>
                        dplyr::ungroup()
                    
                } else if (n_genes > 1) {
                    
                    dat <- 
                        purrr::map_dfc(marker_genes, get_sc_data) |>
                        tibble::add_column(barcode = sc_cells, .before = 1) |>
                        dplyr::left_join(sc_meta, dplyr::join_by(barcode)) |>
                        dplyr::select({{ grouping_var }}, dplyr::all_of(marker_genes)) |>
                        dplyr::group_by(!!var) |>
                        dplyr::summarise_all(list(avg.exp = ~mean(expm1(.)),
                                                  pct.exp = ~mean(. > 0) * 100)) |>
                        dplyr::ungroup() |>
                        tidyr::pivot_longer(-1, 
                                            names_to = c("gene", ".value"),
                                            names_pattern = c("(.+)_(avg.exp|pct.exp)"))
                    
                    ident <- 
                        dat |>
                        dplyr::select({{ grouping_var }}, gene, avg.exp) |>
                        tidyr::pivot_wider(names_from = {{ grouping_var }}, values_from = avg.exp) |>
                        tibble::column_to_rownames("gene") |>
                        dist() |>
                        hclust()
                    
                    dat <- dat |>
                        dplyr::mutate(gene = factor(gene, levels = ident$labels[ident$order]))
                }
                
                dat |>
                    dplyr::group_by(gene) |>
                    dplyr::mutate(avg.exp.scaled = as.numeric(scale(log1p(avg.exp))),
                                  avg.exp.scaled = dplyr::case_when(avg.exp.scaled > 2.5 ~ 2.5,
                                                                    avg.exp.scaled < -2.5 ~ -2.5,
                                                                    .default = avg.exp.scaled)) |>
                    dplyr::ungroup()
                
               
                
               
            })
        
        plot_height <- reactive({
            req(input$markergenes)
            min_len <- 300
            custom_len <- 50 * length(input$markergenes)
            
            return(max(min_len, custom_len))
        })
        
        output$bubbleplot <- 
            renderPlot(res = 96, width = 800, height = function() plot_height(),
                       {
                           req(input$markergenes)
                           ggplot(markers_plot_data(), aes(x = .data[[input$varsc]], y = gene)) +
                               geom_point(aes(size = pct.exp, fill = avg.exp.scaled),
                                          stroke = 0.2, shape = 21) +
                               scale_radius(range = c(0, 6),
                                            breaks = c(0, .25, .5, .75, 1) * 100,
                                            limits = c(0, 100)) +
                               scale_fill_gradient2(low = "Light Sky Blue", 
                                                    mid = "lightyellow", 
                                                    high = "Dark Red",
                                                    midpoint = 0) +
                               theme_minimal() +
                               theme(
                                   axis.text.x = element_text(size = 12),
                                   axis.text.y = element_text(size = 12, face = 'italic'),
                                   panel.grid = element_line(linewidth = .25, color = "grey90"),
                                   legend.title = element_text(size = 12),
                                   legend.key.spacing.y = unit(-.5, "lines")) +
                               guides(fill = guide_colorbar(order = 1, position = "right",
                                                            barwidth = .5, barheight = 5)) +
                               labs(x = NULL, y = NULL, fill = "Scaled\nExpression", size = "%\nExpressed")
                       }
            )
        
        ## Splicing
        get_splicing_data <- 
            reactive({
                req(input$contrast)
                splicing_contrast <- splicing_contrasts[[input$contrast]]
                load(system.file("extdata", paste0(splicing_contrast, ".Rdata"), package = "shinybcells"))
                
                return(list("clusters" = clusters, 
                            "exons_table" = exons_table, 
                            "meta" = meta,
                            "cluster_ids" = cluster_ids,
                            "counts" = counts,
                            "introns" = introns))
            })
        
        splicing_data <- reactive({function() get_splicing_data()}())
        
        output$all_clusters <- 
            DT::renderDT({
                DT::datatable(splicing_data()$clusters[, c("gene", "coord", "N", "FDR", "annotation")],
                              escape = FALSE,
                              rownames = FALSE,
                              colnames = c("Genomic location" = "coord", "Gene" = "gene", "N" = "N", "Annotation" = "annotation", "q" = "FDR"),
                              selection = "single",
                              caption = "Click on a row to plot the corresponding visualization. N: number of introns within a cluster. q: Benjamini–Hochberg q-value.",
                              fillContainer = FALSE,
                              options = list(pageLength = 15,
                                             columnDefs = list(list(className = 'dt-center', targets = 0:4)))
                )
            })
        
        values <- reactiveValues(default = 1) 
        
        observeEvent(input$contrast, {values$default <- 1})
        
        observeEvent(input$all_clusters_rows_selected, {
            values$default <- input$all_clusters_rows_selected
        })
        
        plot_cluster_data <- 
            eventReactive(values$default, {
                sel <- values$default 
                gene  <- splicing_data()$clusters[ sel, ]$gene
                gene <- gsub("<.*?>", "", gene) # strip out html italic tags
                width <- leafviz::getGeneLength(splicing_data()$exons_table, gene)
                clusterID <- splicing_data()$clusters[ sel, ]$clusterID
                coord <- splicing_data()$clusters[ sel, ]$coord
                return(list(gene = gene, width = width, cluster = clusterID, coord = coord) )
            })
        
        output$select_cluster_plot <- 
            renderPlot(width = "auto", height = "auto", res = 90, {
                plotTitle <- c(plot_cluster_data()$gene, as.character(plot_cluster_data()$cluster) )
                leafviz::make_cluster_plot(plot_cluster_data()$cluster,
                                           main_title = plotTitle,
                                           meta = splicing_data()$meta,
                                           cluster_ids = splicing_data()$cluster_ids,
                                           exons_table = splicing_data()$exons_table,
                                           counts = splicing_data()$counts,
                                           introns = splicing_data()$introns)
            })
        
        output$cluster_view <- DT::renderDT({
            clu <- plot_cluster_data()$cluster
            if(!is.null(clu)){
                if(length(splicing_data()$introns)){
                    DT::datatable(leafviz::filter_intron_table(splicing_data()$introns, clu, toSave=FALSE),
                                  rownames = TRUE,
                                  options <- list(searching = FALSE, paging = FALSE, info = FALSE)
                    )
                }
            } else {
                print("no cluster selected!")
            }
            
        }) 
        
        
    }
    
    shinyApp(ui, server, ...)
}