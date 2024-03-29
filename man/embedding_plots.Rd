% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/embedding_plots.R
\name{embedding_plot_2d}
\alias{embedding_plot_2d}
\alias{embedding_plot_2d_ggplot_call}
\alias{pca_plot}
\alias{tsne_plot}
\alias{umap_plot}
\alias{pca_hexbin_plot}
\alias{pca_hexbin_plot_ggplot_call}
\title{PCA, t-SNE and UMAP Plots}
\usage{
embedding_plot_2d(
  fit,
  Y,
  fill = "loading",
  k,
  fill.label,
  ggplot_call = embedding_plot_2d_ggplot_call,
  plot_grid_call = function(plots) do.call(plot_grid, plots)
)

embedding_plot_2d_ggplot_call(
  Y,
  fill,
  fill.type = c("loading", "numeric", "factor", "none"),
  fill.label,
  font.size = 9
)

pca_plot(
  fit,
  Y,
  pcs = 1:2,
  n = 10000,
  fill = "loading",
  k,
  fill.label,
  ggplot_call = embedding_plot_2d_ggplot_call,
  plot_grid_call = function(plots) do.call(plot_grid, plots),
  ...
)

tsne_plot(
  fit,
  Y,
  n = 2000,
  fill = "loading",
  k,
  fill.label,
  ggplot_call = embedding_plot_2d_ggplot_call,
  plot_grid_call = function(plots) do.call(plot_grid, plots),
  ...
)

umap_plot(
  fit,
  Y,
  n = 2000,
  fill = "loading",
  k,
  fill.label,
  ggplot_call = embedding_plot_2d_ggplot_call,
  plot_grid_call = function(plots) do.call(plot_grid, plots),
  ...
)

pca_hexbin_plot(
  fit,
  Y,
  pcs = 1:2,
  bins = 40,
  breaks = c(0, 1, 10, 100, 1000, Inf),
  ggplot_call = pca_hexbin_plot_ggplot_call,
  ...
)

pca_hexbin_plot_ggplot_call(Y, bins, breaks, font.size = 9)
}
\arguments{
\item{fit}{An object of class \dQuote{poisson_nmf_fit} or
\dQuote{multinom_topic_model_fit}.}

\item{Y}{The n x 2 matrix containing the 2-d embedding, where n is
the number of rows in \code{fit$L}. If not provided, the embedding
will be computed automatically.}

\item{fill}{The quantity to map onto the fill colour of the points
in the PCA plot. Set \code{fill = "loading"} to vary the fill
colour according to the loadings (or topic proportions) of the
select topiced or topics. Alternatively, \code{fill} may be set to a
data vector with one entry per row of \code{fit$L}, in which case
the data are mapped to the fill colour of the points. When
\code{fill = "none"}, the fill colour is not varied.}

\item{k}{The dimensions or topics selected by number or name. When
\code{fill = "loading"}, one plot is created per selected dimension
or topic; when \code{fill = "loading"} and \code{k} is not
specified, all dimensions or topics are plotted.}

\item{fill.label}{The label used for the fill colour legend.}

\item{ggplot_call}{The function used to create the plot. Replace
\code{embedding_plot_2d_ggplot_call} or \code{pca_hexbin_plot_ggplot_call}
with your own function to customize the appearance of the plot.}

\item{plot_grid_call}{When \code{fill = "loading"} and multiple
topics (\code{k}) are selected, this is the function used to
arrange the plots into a grid using \code{\link[cowplot]{plot_grid}}.
It should be a function accepting a single argument, \code{plots},
a list of \code{ggplot} objects.}

\item{fill.type}{The type of variable mapped to fill colour. The
fill colour is not varied when \code{fill.type = "none"}.}

\item{font.size}{Font size used in plot.}

\item{pcs}{The two principal components (PCs) to be plotted,
specified by name or number.}

\item{n}{The maximum number of points to plot. If \code{n} is less
than the number of rows of \code{fit$L}, the rows are subsampled at
random. This argument is ignored if \code{Y} is provided.}

\item{\dots}{Additional arguments passed to
\code{\link{pca_from_topics}}, \code{\link{tsne_from_topics}} or
\code{\link{umap_from_topics}}. These additional arguments are only
used if \code{Y} is not provided.}

\item{bins}{Number of bins used to create hexagonal 2-d
histogram. Passed as the \dQuote{bins} argument to
\code{\link[ggplot2]{stat_bin_hex}}.}

\item{breaks}{To produce the hexagonal histogram, the counts are
subdivided into intervals based on \code{breaks}. Passed as the
\dQuote{breaks} argument to \code{\link{cut}}.}
}
\value{
A \code{ggplot} object.
}
\description{
Visualize the structure of the Poisson NMF loadings or
  the multinomial topic model topic proportions by projection onto
  a 2-d surface. \code{plot_hexbin_plot} is most useful for
  visualizing the PCs of a data set with thousands of samples or
  more.
}
\details{
This is a lightweight interface primarily intended to
expedite creation of plots for visualizing the loadings or topic
proportions; most of the heavy lifting is done by
\sQuote{ggplot2}. The 2-d embedding itself is computed by invoking
\code{\link{pca_from_topics}}, \code{\link{tsne_from_topics}} or
\code{\link{umap_from_topics}}. For more control over the plot's
appearance, the plot can be customized by modifying the
\code{ggplot_call} and \code{plot_grid_call} arguments.

An effective 2-d visualization may also require some fine-tunning
of the settings, such as the t-SNE \dQuote{perplexity}, or the
number of samples included in the plot. The PCA, UMAP, t-SNE
settings can be controlled by the additional arguments
(\dots). Alternatively, a 2-d embedding may be pre-computed, and
passed as argument \code{Y}.
}
\examples{
set.seed(1)
data(pbmc_facs)

# Get the Poisson NMF and multinomial topic models fitted to the
# PBMC data.
fit1 <- multinom2poisson(pbmc_facs$fit)
fit2 <- pbmc_facs$fit

# Plot the first two PCs of the loadings matrix (for the
# multinomial topic model, "fit2", the loadings are the topic
# proportions).
subpop <- pbmc_facs$samples$subpop
p1 <- pca_plot(fit1,k = 1)
p2 <- pca_plot(fit2)
p3 <- pca_plot(fit2,fill = "none")
p4 <- pca_plot(fit2,pcs = 3:4,fill = "none")
p5 <- pca_plot(fit2,fill = fit2$L[,1])
p6 <- pca_plot(fit2,fill = subpop)
p7 <- pca_hexbin_plot(fit1)
p8 <- pca_hexbin_plot(fit2)

\donttest{
# Plot the loadings using t-SNE.
p1 <- tsne_plot(fit1,k = 1)
p2 <- tsne_plot(fit2)
p3 <- tsne_plot(fit2,fill = subpop)

# Plot the loadings using UMAP.
p1 <- umap_plot(fit1,k = 1)
p2 <- umap_plot(fit2)
p3 <- umap_plot(fit2,fill = subpop)
}

}
\seealso{
\code{\link{pca_from_topics}},
  \code{\link{tsne_from_topics}},
  \code{\link{umap_from_topics}}
}
