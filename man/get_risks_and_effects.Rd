% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_risks_and_effects}
\alias{get_risks_and_effects}
\title{Helper function to format risk and effect estimates}
\usage{
get_risks_and_effects(
  psi11n,
  psi00n,
  psi10n,
  psi01n,
  cov_matrix,
  cov_matrix_cv
)
}
\arguments{
\item{psi11n}{Estimate of E[Y(1, S(1))]}

\item{psi00n}{Estimate of E[Y(0, S(0))]}

\item{psi10n}{Estimate of E[Y(1, S(0))]}

\item{psi01n}{Estimate of E[Y(0, S(1))]}

\item{cov_matrix}{Covariance matrix}

\item{cov_matrix_cv}{Cross-validated covariance matrix}
}
\value{
A \code{list} with named objects \code{risk}, \code{eff}, \code{eff2}
}
\description{
Helper function to format risk and effect estimates
}
