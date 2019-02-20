#' @title groupedViolins
#' 
#' @description Produce grouped violin plots.  Suitable for 
#' plotting cluster or cell type specific marker genes.
#'
#' @param object Seurat object to plot
#' @param genes List of genes to plot
#' @param marker_list A table of cluster markers as produced by FindMarkers()
#'  or FindAllMarkers or a two column data.frame with the columns "cluster" and 
#'  "gene".  If provided, supercedes the genes argument.  Note: if a marker is 
#'  present in multiple clusters, only the first (as determined by arrange()) is used. 
#'  Default: NULL
#' @param group_by Grouping variable from ident or meta data to use.
#' @param plot_title What it says on the tin.  Default: NULL
#' @param assay_use Assay to plot.  Default: "RNA"
#' @param flip Should genes be displayed along the y-axis? Default: TRUE
#'
#' @importFrom Seurat FetchData
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#' @importFrom plyr mapvalues
#' @import ggplot2
#'
#' @return
#' @export
#'
#' @examples
groupedViolins <- function(object, 
                         genes = NULL, 
                         marker_list = NULL,
                         group_by,
                         plot_title = NULL,
                         assay_use = "RNA",
                         flip = TRUE){
  
  # if there is a simpler/better way to handle arguments 
  # being either quoted or unquoted, I'd love to hear about it.
  if (is.null(genes) & is.null(marker_list)){
    stop("No genes or marker_list provided.  What am I supposed to plot?")
  }
  try(
    if (is.character(group_by)) {
      group_by <- as.name(substitute(group_by))
    }, silent = TRUE
  )

  if (!is.null(marker_list)){
    marker_list %<>% select(cluster, gene) %>%
      arrange(gene) %>%
      group_by(gene) %>%
      slice(1)
    genes <- marker_list %>% pull(gene)
  }
  
  group_by <- enquo(group_by)
  type_expr <- FetchData(object, 
                         vars = c(unique(genes), quo_name(group_by)))
  type_expr %<>% 
    rownames_to_column('name') %>% 
    gather(key = "gene",
           value = "value",
           -name,
           -!!group_by)

  if (!is.null(marker_list)){
    type_expr$cluster_marker <- mapvalues(x = type_expr[["gene"]],
                                          from = marker_list[["gene"]],
                                          to = marker_list[["cluster"]])
    # fix the ordering of the clusters in the event they are all numbers
    # so that 0 < 1 < 2 ... < n + 1 instead or 0 < 1 < 10 < 2
    if(all(grep(pattern = "[::digit::]+", x = type_expr[["cluster_marker"]]))){
    type_expr[["cluster_marker"]] %<>% as.integer() %>% as.factor()
    }
  }
  
  type_violins <- type_expr %>%
    ggplot() +
    theme(legend.position = "none",
          strip.text = element_text(size = 9),
          axis.text = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(.05,
                               "lines"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(face = "italic")) +
    labs(x = "", y = "")
  
  # if (!is.null(marker_list)){
  #   type_violins <- type_violins + 
  # } else {
  #   type_violins <- type_violins +
  # }

  if (isTRUE(flip)){
    type_violins <- type_violins + coord_flip()
  }

  if (!is.null(marker_list)){
    type_violins <- type_violins +
      geom_violin(aes(x = !!group_by,
                      y = value, 
                      fill = cluster_marker),
                  scale = "width",
                  trim = TRUE) +
      facet_grid(~cluster_marker+gene,
                 scales = "free")
  } else {
    type_violins <- type_violins +
      geom_violin(aes(x = !!group_by,
                      y = value, 
                      fill = gene),
                  scale = "width",
                  trim = TRUE) +
      facet_grid(~gene, scales = "free")
  }

  if(!is.null(plot_title)){
    type_violins <- type_violins + labs(title = plot_title)
  }
  type_violins
}


marker_by_cluster_violins <- function(object, 
                                      markers,
                                      ...){
  lapply(unique(markers$cluster), function(i){
    j <- markers %>% filter(cluster == i) %>% pull(gene)
    lineage_violins(object = object, genes = j, ...) +
      labs(title = glue::glue("Cluster {i}"))
  })
}
