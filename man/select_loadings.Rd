% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.R
\name{select.poisson_nmf_fit}
\alias{select.poisson_nmf_fit}
\alias{select}
\alias{select.multinom_topic_model_fit}
\alias{select_loadings}
\title{Extract or Re-order Data Rows in Poisson NMF or Multinomial Topic Model Fit}
\usage{
\method{select}{poisson_nmf_fit}(.data, loadings, ...)

\method{select}{multinom_topic_model_fit}(.data, loadings, ...)

select_loadings(.data, loadings, ...)
}
\arguments{
\item{.data}{Poisson NMF or Multinomial Topic Model fit; that is,
an object of class \dQuote{poisson_nmf_fit} or
\dQuote{multinom_topic_model_fit}, such as an output from
\code{\link{fit_poisson_nmf}} or \code{\link{fit_topic_model}}.}

\item{loadings}{Indices (names or numbers) giving data rows to
keep. If not specified, all rows are kept.}

\item{\dots}{Other arguments passed to the generic select function.}
}
\value{
A Poisson NMF or multinomial topic model fit containing the
  selected data rows only.
}
\description{
This function can be used to extract estimates for a
  subset of the count data, or to re-order the rows of the loadings
  matrix.
}
