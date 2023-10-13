
## Plotting ----------

#' Plot pairwise enrichment heatmap
#'
#' @inheritParams pairwise_score
#'
#' @importFrom ggplot2 ggplot aes geom_tile labs theme_minimal scale_fill_viridis_c geom_text coord_fixed theme scale_y_discrete element_text
#' @return ggplot object
#' @export
#'
#' @examples
#' gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
#' gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")
#'
#' gr_b1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")
#' gr_b2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 25), strand = "+")
#'
#' gsize <- 30
#'
#' score_func <- purrr::partial(fc_enrichment, genome_size = gsize, ignore.strand = TRUE)
#'
#' plot_pairwise_score(
#'   list("A" = gr_a1, "B" = gr_a2),
#'   list("X" = gr_b1, "Y" = gr_b2),
#'   score_func
#' )
plot_pairwise_score <- function(grlist_1, grlist_2, score_func) {
  # seqInfo merge issues this annoying warning everytime two GRanges have
  # seqnames not in common, but this is a perfectly valid situation, because not
  # every GRanges object will encompass whole genome, so it does not really add
  # value and tends to fill console with warnings
  values <- suppressWarnings(
    pairwise_score(grlist_1, grlist_2, score_func)
  )

  ggplot(values, aes(x=!!quote(gr1), y=!!quote(gr2), fill=!!quote(score), label = round(!!quote(score), 2))) +
    geom_tile(color="white", linewidth = 1) +
    geom_text(color = "#aaaaaa") +
    theme_minimal() +
    scale_fill_viridis_c(option = "D") +
    coord_fixed() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
    labs(x= "",
         y = "",
         title = "Pairwise enrichment heatmap",
         caption = .package_caption(list(), list(), TRUE)) +
    scale_y_discrete(limits=rev)
}

#' Plot triangular pairwise enrichment heatmap for a single GRanges list
#'
#' @importFrom ggplot2 ggplot aes geom_tile labs theme_minimal
#' @importFrom ggplot2 scale_fill_viridis_c geom_text coord_fixed
#' @importFrom ggplot2 theme scale_y_discrete element_text theme element_blank
#' @inheritParams combinations_score
#' @return A ggplot object
#' @export
#'
#' @examples
#' gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
#' gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")
#'
#' gr_a3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")
#'
#' score_func <- purrr::partial(fc_enrichment, genome_size = 30, ignore.strand = TRUE)
#'
#' plot_combinations_score(
#'   list("A" = gr_a1, "B" = gr_a2, "C" = gr_a3),
#'   score_func
#' )
plot_combinations_score <- function(grlist, score_func) {
  # seqInfo merge issues this annoying warning everytime two GRanges have no
  # seqnames in common, but this is a perfectly valid situation, because not
  # every GRanges object will encompass whole genome, so it does not really add
  values <- suppressWarnings(
    combinations_score(grlist, score_func)
  )

  ggplot(values, aes(
      x = !!quote(gr1),
      y = !!quote(gr2),
      fill = !!quote(score),
      label = round(!!quote(score), 2))) +
    geom_tile(color="white", linewidth = 1) +
    geom_text(color="#aaaaaa") +
    theme_minimal() +
    scale_fill_viridis_c(option = "D", limits = c(0, max(values$score))) +
    coord_fixed() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
    labs(x= "",
         y = "",
         title = "Pairwise scored heatmap",
         caption = .package_caption(list(), list(), TRUE)) +
    scale_y_discrete(limits=rev) +
    theme(panel.grid = element_blank())
}

#' Plot a Euler diagram with intersections between GRanges objects
#'
#' Further parameters can be passed to the plot function (see plot.euler)
#'
#' @param grlist List of GRanges objects
#' @param names Names to give the sets
#' @param ignore.strand Whether to ignore strand for intersections
#' @param fills Color fills for each set. Eulerr interpolates the intersections
#' @importFrom eulerr euler
#' @return A euler.diagram object that can be further processed, or plotted via print, ggsave
#' @export
#'
#' @examples
#' gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
#' gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
#' gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 18), strand = "+")
#' plot_euler(list(gr_1, gr_2, gr_3), names = c("a", "b", "c"))
plot_euler <- function(grlist, names = NULL, ignore.strand = TRUE, fills = NULL) {
  v_int <- calculate_venn_intersections(
    grlist,
    names = names,
    ignore.strand = ignore.strand
  )

  plot_specs <- eulerr::euler(v_int, input = "union")

  if (is.null(fills)) {
    plot(plot_specs, quantities = TRUE)
  } else {
    plot(plot_specs, fills = fills, quantities = TRUE)
  }
}

## Helpers ----------

#' Make a string out of a named list.
#'
#' @param named_list A named list
#' @return A string
.key_value_string <- function(named_list) {
  paste(names(named_list), lapply(named_list, .format_value), sep = ":", collapse = ", ")
}


.format_value <- function(v) {
  if(is.numeric(v)) {
    v <- sprintf("%.4f", v)
    v <- sub("\\.?0+$", "", v)
  }
  v
}

#' Make a string out of a named list. Split into lines if too wide.
#'
#' @param named_list A named list
#'
#' @return A string
.limited_size_caption_line <- function(named_list) {
  size_limit <- 3
  chunks <- split(named_list, ceiling(seq_along(named_list)/size_limit))
  paste(vapply(chunks, .key_value_string, character(1)), collapse="\n")
}


#' Make a string to put as caption in verbose mode. Includes system date.
#'
#' @param params Named list with relevant parameters and their values
#' @param outcome Named values with relevant outcomes and their values
#' @param verbose Logical. If TRUE all information is printed. If FALSE returns
#'   a blank caption.
#' @importFrom utils packageVersion
#' @return A caption string
.package_caption <- function(params, outcome, verbose) {
  if (verbose) {
    verbose_params <- .limited_size_caption_line(params)
    verbose_crop <- .limited_size_caption_line(outcome)

    date <- format(Sys.time(), "%a %b %d %X %Y")
    pkg_version <- paste("bedscout v.", packageVersion("bedscout"))
    date <- paste(date, pkg_version, sep = ' - ')

    paste(verbose_params, verbose_crop, date, sep = "\n\n")
  }
}

