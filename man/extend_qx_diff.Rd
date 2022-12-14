% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extend_qx_diff.R
\name{extend_qx_diff}
\alias{extend_qx_diff}
\title{Extend qx up to age 105 using HMD regression parameters}
\usage{
extend_qx_diff(empir_lt, hmd_qx_results, by_vars)
}
\arguments{
\item{empir_lt}{data.table with columns: ihme_loc_id, sex, year, age (numeric), mx, ax}

\item{hmd_qx_results}{data.table with variables sex, age, slope, intercept}

\item{by_vars}{character vector containing the variable names of all variables that uniquely identify the observations (except for age)}
}
\value{
returns empir_lt with the same variables, but modified qx extended to age 105
}
\description{
\deqn{
logit({}_{5}q_{x+5}) - logit({}_{5}q_x) = \beta_0 + \beta_1 * (logit({}_{5}q_x) - logit({}_{5}q{x-5})) + \text{location}_{RE}
}
}
\details{
Recursive prediction
}
