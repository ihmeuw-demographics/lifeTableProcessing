#' Prepare HMD data
#'
#' Compile "by statistic" raw data from the [Human Mortality database][https://mortality.org/]
#'
#' Currently there are functions to prep data in 5-year age groups and single
#' year time periods for:
#'
#' * deaths: [prep_hmd_deaths_5x1()]
#' * population: [prep_hmd_pop_5x1()]
#' * life tables: [prep_hmd_period_lt_5x1()]
#'
#' @param dir_input Path to directory containing raw location-specific files for
#'   all sexes.
#' @param dir_input_male Path to directory containing raw location-specific
#'   files for males.
#' @param dir_input_female Path to directory containing raw location-specific
#'   files for females.
#'
#' @return [data.table()] of HMD data for all locations.
#'
#' @name prep_hmd
NULL


#' @name prep_hmd
#' @export
prep_hmd_deaths_5x1 <- function(dir_input) {

  dt <- load_raw_hmd_files(dir_input, skip = 2, na = ".", col_types = "icnnn")
  setnames(dt, tolower(colnames(dt)))
  prep_hmd_age(dt, "age")
  dt[, total := NULL]

  dt_long <- data.table::melt(
    dt,
    measure.vars = c("female", "male"),
    variable.name = "sex",
    value.name = "deaths",
    variable.factor = FALSE
  )

  dt_long

}



#' @name prep_hmd
#' @export
prep_hmd_pop_5x1 <- function(dir_input) {

  dt <- load_raw_hmd_files(dir_input, skip = 2, na = ".", col_types = "ccnnn")
  setnames(dt, tolower(colnames(dt)))
  prep_hmd_age(dt, "age")
  dt[, total := NULL]

  dt[year %like% "-$", territory_change_status := "before"]
  dt[year %like% "\\+$", territory_change_status := "after"]
  dt[, year := as.integer(substr(year, 1, 4))]

  dt_long <- data.table::melt(
    dt,
    measure.vars = c("female", "male"),
    variable.name = "sex",
    value.name = "population",
    variable.factor = FALSE
  )

  dt_long

}

#' @name prep_hmd
#' @export
prep_hmd_period_lt_5x1 <- function(dir_input_male, dir_input_female) {

  dt <-
    list(
      male = dir_input_male,
      female = dir_input_female
    ) |>
    lapply(
      load_raw_hmd_files,
      skip = 2, na = ".", col_types = "icnnnnnnnn"
    ) |>
    rbindlist(idcol = "sex", use.names = TRUE)

  prep_hmd_age(dt, "Age")
  setnames(dt, "Year", "year")
  setnames(dt, "Lx", "nLx")
  setcolorder(dt, c("hmd_loc_id", "year", "sex", "age_start", "age_end"))

  dt

}

#' Load raw HMD files
#'
#' @param dir Path do directory of raw location-specific HMD files
#' @param ... Additional parameters passed to [readr::read_table()]
#'
#' @return [data.table()] compiled from all location files in `dir`.
load_raw_hmd_files <- function(dir, ...) {

  hmd_raw_files <- fs::dir_ls(dir)

  names(hmd_raw_files) <-
    hmd_raw_files |>
    names() |>
    basename() |>
    tools::file_path_sans_ext() |>
    sub("\\..*","", x = _)

  dt <- hmd_raw_files |>
    lapply(\(f) readr::read_table(file = f, ...)) |>
    data.table::rbindlist(idcol = "hmd_loc_id", use.names = TRUE) |>
    data.table::setDT()

  return(dt)

}

#' Prepare age information from HMD data
#'
#' @param dt [data.table()] of HMD data.
#' @param age_col name of column representing age
#'
#' @return Invisibly returns `dt` with `age_col` removed in place of `age_start` and
#'   `age_end`.
prep_hmd_age <- function(dt, age_col) {

  dt[, c("age_start", "age_end") := data.table::tstrsplit(get(age_col), "-", fixed = TRUE)]

  dt[, `:=`(
    age_start = as.numeric(gsub("\\+$", "", age_start)),
    age_end = as.numeric(age_end) + 1
  )]

  dt[age_start == 0, age_end := 1]
  dt[age_start == 110, age_end := Inf]
  dt[, (age_col) := NULL]

  invisible(dt)

}
