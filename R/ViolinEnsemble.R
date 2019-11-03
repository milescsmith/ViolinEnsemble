#' @title ViolinEnsemble
#'
#' @description Produce grouped violin plots.  Suitable for
#' plotting cluster or cell type specific marker features.
#'
#' @param object Seurat object to plot
#' @param features List of features to plot
#' @param marker_list A table of cluster markers as produced by FindMarkers()
#' or FindAllMarkers or a two column data.frame with the columns "cluster" and
#' "feature".  If provided, supercedes the features argument.  Note: if a marker is
#' present in multiple clusters, only the first (as determined by arrange()) is used.
#' Default: NULL
#' @param grouping_var Grouping variable from ident or meta data to use. If NULL,
#' the current active.ident is used.  Default: NULL
#' @param cluster_x Arrange the x-axis variables using hierarchical clustering. Default: TRUE
#' @param cluster_y Arrange the y-axis variables using hierarchical clustering. Default: FALSE
#' @param show_points Display data points in addition to violin using geom_jitter? Default: FALSE
#' @param alpha Alpha value to use with geom_jitter. Default: 1
#' @param pt_size Point size to use with geom_jitter. Default: 0.25
#' @param plot_title What it says on the tin.  Default: NULL
#' @param assay_use Assay to plot.  Default: "RNA"
#' @param slot_use Slot to draw data from.  Default: "data"
#' @param sort If TRUE, attempt to sort the groups alpha-numerically. Default: FALSE
#' @param flip Should features be displayed along the y-axis? Default: TRUE
#'
#' @importFrom Seurat GetAssayData FetchData DefaultAssay<-
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr group_by arrange select slice mutate_if inner_join pull
#' @importFrom stringr str_detect
#' @importFrom magrittr %>% %<>%
#' @importFrom gtools mixedsort
#' @importFrom stats reorder
#' @import ggplot2
#'
#' @return
#' @export
#'
#' @examples
ViolinEnsemble <- function(object,
                           features = NULL,
                           marker_list = NULL,
                           grouping_var = NULL,
                           cluster_x = FALSE,
                           cluster_y = FALSE,
                           show_points = FALSE,
                           alpha = 0.5,
                           pt_size = 0.25,
                           plot_title = NULL,
                           assay_use = "RNA",
                           slot_use = "data",
                           sort = FALSE,
                           flip = TRUE){

  # if there is a simpler/better way to handle arguments
  # being either quoted or unquoted, I'd love to hear about it.
  if (is.null(features) & is.null(marker_list)){
    stop("No features or marker_list provided.  What am I supposed to plot?")
  }
  try(
    if (is.null(grouping_var)){
      grouping_var <- "ident"
    }, silent = TRUE)
  try(
    if (is.character(grouping_var)) {
      grouping_var <- as.name(substitute(grouping_var))
    }, silent = TRUE
  )

  if (!is.null(marker_list)){
    marker_list %<>% select(cluster, feature) %>%
      mutate_if(is.factor, as.character) %>%
      arrange(feature) %>%
      group_by(feature) %>%
      slice(1)
    features <- marker_list %>% pull(feature) %>% as.character()
  }

  grouping_var <- enquo(grouping_var)

  if (assay_use %in% names(object)){
    DefaultAssay(object) <- assay_use
  } else {
    stop("That assay does not appear to be present in your data object")
  }

  type_expr <- FetchData(object,
                         vars = c(unique(features), quo_name(grouping_var)),
                         slot = slot_use)
  type_expr %<>%
    rownames_to_column('name') %>%
    gather(key = "feature",
           value = "value",
           -name,
           -!!grouping_var)

  if (!is.null(marker_list)){
    type_expr %<>% inner_join(marker_list, by = "feature")
    # fix the ordering of the clusters in the event they are all numbers
    # so that 0 < 1 < 2 ... < n + 1 instead or 0 < 1 < 10 < 2
    if(all(str_detect(pattern = "[0-9]+",
                      string = type_expr[["cluster"]]))){
      type_expr[["cluster"]] %<>% as.integer() %>% as.factor()
    }
  }

  if (isTRUE(sort)){
    try(
      if(is.factor(type_expr[[quo_name(grouping_var)]])){
        type_expr[[quo_name(grouping_var)]] %<>% reorder()
        type_expr %<>% arrange(!!grouping_var)
      })
  }

  if (isTRUE(cluster_x)) {

    feature_dendro <- type_expr %>%
      group_by(ident, feature) %>%
      summarize(avg = mean(value)) %>%
      pivot_wider(names_from = feature,
                  values_from = avg) %>%
      t() %>%
      dist() %>%
      hclust() %>%
      as.dendrogram()

      type_expr$feature <- factor(type_expr$feature,
                                 levels = labels(feature_dendro),
                                 ordered = TRUE)
  }

  if (isTRUE(cluster_y)) {

    id_dendro <- type_expr %>%
      group_by(ident, feature) %>%
      summarize(avg = mean(value)) %>%
      pivot_wider(names_from = feature,
                  values_from = avg) %>%
      dist() %>%
      hclust() %>%
      as.dendrogram()

    type_expr$ident <- factor(type_expr$ident,
                                 levels = labels(id_dendro),
                                 ordered = TRUE)
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
      geom_violin(aes(x = !!grouping_var,
                      y = value,
                      fill = cluster),
                  scale = "width",
                  trim = TRUE) +
      facet_grid(~cluster+feature,
                 scales = "free")
  } else {
    type_violins <- type_violins +
      geom_violin(aes(x = !!grouping_var,
                      y = value,
                      fill = feature),
                  scale = "width",
                  trim = TRUE) +
      facet_grid(~feature, scales = "free")
  }

  if (show_points){
    type_violins <- type_violins + geom_jitter(aes(x = !!grouping_var,
                                                   y = value),
                                               alpha = alpha,
                                               size = pt_size
                                               )
  }

  if(!is.null(plot_title)){
    type_violins <- type_violins + labs(title = plot_title)
  }
  type_violins
}
