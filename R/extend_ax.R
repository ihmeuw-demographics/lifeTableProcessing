#' Extend ax values for ages 75-105 using HMD regression parameters
#'
#' @param empir_lt data.table with columns: ihme_loc_id, sex, year, age (numeric), mx, ax
#' @param hmd_ax_results data.table with variables par_mx, par_smx, par_con, sex, age
#'
#' @return returns empir_lt with the same variables, but modified values for ax for ages 75 and above
#'
#' @export
extend_ax <- function(empir_lt, hmd_ax_results) {

  empir_lt <- copy(empir_lt)

  empir_lt[
    hmd_ax_results,
    ax := i.par_mx * mx + i.par_smx * (mx ^ 2) + i.par_con,
    on = c("sex", "age_start", "age_end")
  ]

  return(empir_lt)

}
