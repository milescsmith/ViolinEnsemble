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
#' @param group_var Grouping variable from ident or meta data to use. If NULL,
#' the current active.ident is used.  Default: NULL
#' @param plot_title What it says on the tin.  Default: NULL
#' @param assay_use Assay to plot.  Default: "RNA"
#' @param slot_use Slot to draw data from.  Default: "data"
#' @param flip Should genes be displayed along the y-axis? Default: TRUE
#'
#' @importFrom Seurat GetAssayData
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr group_by arrange select slice
#' @importFrom stringr str_detect
#' @importFrom magrittr %>% %<>%
#' @import ggplot2
#'
#' @return
#' @export
#'
#' @examples
groupedViolins <- function(object,
                         genes = NULL,
                         marker_list = NULL,
                         group_var = NULL,
                         plot_title = NULL,
                         assay_use = "RNA",
                         slot_use = "data",
                         flip = TRUE){

  # if there is a simpler/better way to handle arguments
  # being either quoted or unquoted, I'd love to hear about it.
  if (is.null(genes) & is.null(marker_list)){
    stop("No genes or marker_list provided.  What am I supposed to plot?")
  }
  try(
    if (is.null(group_var)){
      group_var <- "ident"
    }, silent = TRUE)
  try(
    if (is.character(group_var)) {
      group_var <- as.name(substitute(group_var))
    }, silent = TRUE
  )

  if (!is.null(marker_list)){
    marker_list %<>% select(cluster, gene) %>%
      mutate_if(is.factor, as.character) %>%
      arrange(gene) %>%
      group_by(gene) %>%
      slice(1)
    genes <- marker_list %>% pull(gene) %>% as.character()
  }

  group_var <- enquo(group_var)
  
  if (assay_use %in% names(object)){
    DefaultAssay(object) <- assay_use
  } else {
    stop("That assay does not appear to be present in your data object")
  }
  
  type_expr <- FetchData(object,
                         vars = c(unique(genes), quo_name(group_var)),
                         slot = slot_use)
  type_expr %<>%
    rownames_to_column('name') %>%
    gather(key = "gene",
           value = "value",
           -name,
           -!!group_var)

  if (!is.null(marker_list)){
    type_expr %<>% inner_join(marker_list, by = "gene")
    # fix the ordering of the clusters in the event they are all numbers
    # so that 0 < 1 < 2 ... < n + 1 instead or 0 < 1 < 10 < 2
    if(all(str_detect(pattern = "[0-9]+",
                      string = type_expr[["cluster"]]))){
      type_expr[["cluster"]] %<>% as.integer() %>% as.factor()
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

  if (isTRUE(flip)){
    type_violins <- type_violins + coord_flip()
  }

  if (!is.null(marker_list)){
    type_violins <- type_violins +
      geom_violin(aes(x = !!group_var,
                      y = value,
                      fill = cluster),
                  scale = "width",
                  trim = TRUE) +
      facet_grid(~cluster+gene,
                 scales = "free")
  } else {
    type_violins <- type_violins +
      geom_violin(aes(x = !!group_var,
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
