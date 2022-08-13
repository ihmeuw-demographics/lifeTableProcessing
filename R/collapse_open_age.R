#' Collapse open age of life table
#'
#' If the given open age is less that the current open age of the life table,
#' collapse deaths and population counts for higher age groups into the new
#' open age.
#'
#' @param lt `data.table()` representing a single life table to (possibly)
#'   collapse.
#' @param open_age New open age interval start.
#' @param catch_lower_input `logical()` determining if an error should be thrown
#'   when the input life table's open age is lower than the given open age.
#'
#' @return Life table with new open age interval, or original life table if the
#'   original open age is lower than the given open age.
#' @export
#'
#' @examples
collapse_open_age <- function(lt,
                              open_age,
                              catch_lower_input = FALSE) {

  const_lt_cols <- setdiff(
    colnames(lt),
    c("age_start", "age_end", "deaths", "population")
  )

  original_open_age <- get_open_age(lt)

  if (original_open_age <= open_age) {

    msg <- paste(
      "Original open age of", original_open_age,
      "is less than or equal to given open age of", open_age
    )

    if (!catch_lower_input) {
      message(msg)
      return(data.table::copy(lt))
    } else {
      stop(msg)
    }

  }

  lt_open <- lt[
    age_start >= open_age,
    .(
      age_start = open_age,
      age_end = Inf,
      deaths = sum(deaths),
      population = sum(population)
    ),
    by = const_lt_cols
  ]

  rbind(lt[age_start < open_age], lt_open)

}
