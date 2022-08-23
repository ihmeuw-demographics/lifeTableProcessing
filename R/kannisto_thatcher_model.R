#' Kannisto-Thatcher model in probability of death space
#'
#' \deqn{
#' q_x = 1 - (\frac{\alpha * \exp(\beta x) + 1}{\alpha * \exp(\beta * (x+5)) + 1})^{1/\beta}
#' }
#'
#' @param x age group start
#' @param params [`numeric(2)`] vector of Kannisto-Thatcher model parameters
#'   \eqn{\alpha} (a) and \eqn{\beta} (b).
#'
#' @return Prediction of Kannisto-Thatcher model at age x
#'
#' @export
model_kt_qx <- function(x, params) {

  a <- params[1]
  b <- params[2]

  numerator <- a * exp(b * x) + 1
  denominator <- a * exp(b * (x + 5)) + 1

  1 - (numerator / denominator)^(1/b)

}


#' Calculate sum of squared-errors for qx Kannisto-Thatcher model
#'
#' \deqn{
#' \text{SSE} = \sum ((q_x - q_{x,\text{pred}})^2)
#' }
#'
#' @param params [`numeric(2)`] vector of Kannisto-Thatcher model parameters
#'   \eqn{\alpha} (a) and \eqn{\beta} (b).
#' @param data [data.table()] Life table with at columns `age_start` and `qx`.
#'
#' @return Sum of square errors of \eqn{q_x}.
#'
#' @export
sse_kt_qx <- function(params, data) {

  pred_qx <- model_kt_qx(data$age_start, params)
  squared_error <- (data$qx - pred_qx)^2

  sum(squared_error)

}


#' Calculate Kannisto-Thatcher extension parameters
#'
#' Calculate optimal parameters \eqn{\alpha} and \eqn{\beta} in the
#' Kannisto-Thatcher model for use in extending a single life table.
#'
#' @param lt [data.table()] Life table with columns `age_start`, `mx`, and `qx`.
#' @param age_start_range [`numeric(2)`] vector of lower and upper inclusive
#'   age group start bounds to calculate parameters from.
#' @param fn Optimization function passed to [`optim()`].
#'
#' @details
#'
#' ## Parameter selection
#'
#' Initial parameters are derived from the model:
#'
#' \deqn{\text{logit}(m_x) \sim \ln(\alpha) + \beta * \text{age}}
#'
#' Optimization is then done using [`optim()`] to minimize the function `fn`.
#' `fn` must have the method signature `fn(params, data)`, where:
#'
#'  * `params` is a named vector of the parameters \eqn{\alpha} (`a`) and
#'     \eqn{\beta} (`b`)
#'  * `data` is a life table similar to `lt`
#'
#' @return `numeric(2)` a named vector of optimal parameters `a` and `b`.
#'
#' @export
calculate_kt_params <- function(lt, age_start_range, fn) {

  data <- lt[age_start %between% age_start_range]

  fit_init <- glm(
    z ~ age_start,
    data = data[, .(age_start, z = demUtils::logit(mx))]
  )

  param_init <- c(
    a = unname(exp(coef(fit_init)[1])),
    b = unname(coef(fit_init)[2])
  )

  fit <- optim(
    par = param_init,
    fn = fn,
    data = data,
    method = "L-BFGS-B",
    lower = c(a = 1e-8, b = 1e-8),
    upper = c(a = 5, b = 5),
    control = list(maxit = 1e4)
  )

  fit$par

}


#' Predict probability of death with the Kannisto-Thatcher model
#'
#' @param params [`numeric(2)`] vector of Kannisto-Thatcher model parameters
#'   \eqn{\alpha} (a) and \eqn{\beta} (b).
#' @param base_age start of age group to begin making predictions from/
#' @param open_age_start start of pen age interval (\eqn{q_x = 1}).
#'
#' @return [`data.table()`] life table of \eqn{q_x} values from `base_age` to
#'  the age group preceding the open age interval.
#'
#' @export
predict_kt_qx <- function(params, base_age, open_age_start = 110) {

  stopifnot(base_age < (open_age_start - 5))

  new_ages <- seq(base_age, open_age_start - 5, 5)
  pred_qx <- model_kt_qx(new_ages, params)

  data.table(
    age_start = new_ages,
    age_end = new_ages + 5,
    qx = pred_qx
  )

}


#' Extend qx with the Kannisto-Thatcher model
#'
#' @param lt [`data.table()`] Life table with columns `age_start`, `mx`, `qx`, and
#'   `id_cols` identifying a single life table.
#' @param id_cols [`character()`] vector of column names identifying a single
#'   life table
#' @param age_start_range [`numeric(2)`] vector of lower and upper inclusive
#'   age group start bounds to calculate model parameters from.
#' @param open_age_start [`numeric(1)`] start of the open age interval to
#'   predict out to.
#'
#' @return [`data.table()`] life table(s) with \eqn{q_x} extended out to
#'   `open_age_start`.
#'
#' @export
extend_qx_kt <- function(lt, id_cols, age_start_range, open_age_start = 110) {

  lt <- prep_lt_extension(lt, id_cols, age_start_cutoff = 65)

  pred_func <- function(data) {
    params <- unlist(data$params)
    base_age <- data$base_age
    predict_kt_qx(params, base_age, open_age_start)
  }

  lt_params <- lt[
    j = .(
      params = lapply(
        list(.SD),
        calculate_kt_params,
        age_start_range = age_start_range,
        fn = sse_kt_qx
      ),
      base_age = max(age_start) + 5
    ),
    by = id_cols
  ]

  lt_predict <- lt_params[
    j = .(dt = lapply(list(.SD), pred_func)),
    by = id_cols
  ]

  lt_predict <- setDT(tidyr::unnest(lt_predict, dt))

  lt_extended <- rbind(lt, lt_predict, use.names = TRUE, fill = TRUE)
  setorderv(lt_extended, id_cols)

}


prep_lt_extension <- function(lt, id_cols, age_start_cutoff) {

  data <- copy(lt)

  # Remove any old age groups with qx > 1
  data[
    age_start >= age_start_cutoff & qx > 1,
    age_qx_over_1 := min(age_start),
    by = id_cols
  ]
  data[is.na(age_qx_over_1), age_qx_over_1 := 999]
  data[, min_age_qx_over_1 := min(age_qx_over_1, na.rm = T), by = id_cols]
  data <- data[age_start < min_age_qx_over_1]

  # Drop old age groups after qx starts to decrease
  data[
    age_start >= age_start_cutoff,
    qx_change := shift(qx, 1, type = "lead") - qx,
    by = id_cols
  ]
  data[
    age_start >= age_start_cutoff & qx_change < 0,
    min_age_qx_change_neg := min(age_start),
    by = id_cols
  ]

  change_neg <- data[
    age_start == min_age_qx_change_neg,
    c(..id_cols, "min_age_qx_change_neg")
  ]

  data[, min_age_qx_change_neg := NULL]

  if (nrow(change_neg) > 0) {
    data <- merge(data, change_neg, by = id_cols, all.x=T)
    data <- data[is.na(min_age_qx_change_neg) | age_start < min_age_qx_change_neg]
  }

  data[, -c("age_qx_over_1", "min_age_qx_over_1", "qx_change", "min_age_qx_change_neg")]

}
