#' Threshold Indicator Taxa Analysis
#'
#' Uses indicator species scores across binary partitions of a sample set to
#' detect congruence in taxon-specific changes of abundance and occurrence
#' frequency along an environmental gradient as evidence of an ecological
#' community threshold.
#'
#' Relevant references include Baker and King (2010)
#' <doi:10.1111/j.2041-210X.2009.00007.x>, King and Baker (2010)
#' <doi:10.1899/09-144.1>, and Baker and King (2013) <doi:10.1899/12-142.1>.
#'
#' @import parallel ggplot2
#' @docType package
#' @name TITAN2
#' @importFrom graphics axis box legend mtext par plot points polygon segments
#'   symbols
#' @importFrom stats approxfun median quantile runif sd reorder
#' @importFrom utils read.table write.table
#' @importFrom glue glue
#' @importFrom cowplot plot_grid
#' @importFrom dplyr %>% data_frame as_data_frame mutate filter case_when slice
#'   bind_rows rename arrange n slice
#' @importFrom purrr pluck map map_dfr
#' @importFrom ggridges geom_density_ridges
#' @importFrom cli cli_alert_info cli_alert col_red col_green col_black
#'   cli_alert_success cli_warn cli_progress_bar cli_progress_update
#'   cli_progress_done cli_abort cli_alert_warning
#' @aliases TITAN2 package-TITAN2
NULL
