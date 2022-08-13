#' Get open age interval of a life table
#'
#' @param lt `data.table()` representing a single life table.
#'
#' @return Start of open age interval
#'
#' @export
get_open_age <- function(lt) {

  lt[is.infinite(age_end), age_start]

}
