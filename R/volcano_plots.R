#' @rdname volcano_plot
#'
#' @title Volcano Plots for Visualizing Results of Differential Expression Analysis
#'
#' @description Create a \dQuote{volcano} plot to visualize the
#'   results of a differential count analysis using a topic model. Here,
#'   the volcano plot is a scatterplot in which the posterior mean
#'   log-fold change (LFC), estimated by running the methods implemented
#'   in \code{\link{de_analysis}}, is plotted against the estimated
#'   z-score. Variations on this volcano plot may also be created, for
#'   example by showing f0 (the null-model estimates) instead of the
#'   z-scores. Use \code{volcano_plotly} to create an interactive
#'   volcano plot.
#'
#' @details Interactive volcano plots can be created using the
#'   \sQuote{plotly} package. The \dQuote{hover text} shows the label
#'   and detailed LFC statistics.
#' 
#' @param de An object of class \dQuote{topic_model_de_analysis},
#'   usually an output from \code{\link{de_analysis}}. It is better to
#'   run \code{de_analysis} with \code{shrink.method = "ash"} so that
#'   the points in the volcano plot can be coloured by their local false
#'   sign rate (lfsr).
#'
#' @param k The topic, selected by number or name.
#'
#' @param labels Character vector specifying how the points in the
#'   volcano plot are labeled. This should be a character vector with
#'   one entry per LFC estimate (row of \code{de$postmean}). When not
#'   specified, the row names of \code{de$postmean} are used. When
#'   available. labels are added to the plot using
#'   \code{\link[ggrepel]{geom_text_repel}}.
#'
#' @param y Which quantity to plot on along the vertical (y)
#'   axis. Choose \code{y = "z"} to show z-score magnitudes, or \code{y
#'   = "f0"} to show the estimated expression levels.
#'
#' @param do.label The function used to deetermine which LFC estimates
#'   to label. Replace \code{volcano_plot_do_label_default} with your
#'   own function to customize the labeling of points in the volcano
#'   plot.
#' 
#' @param ymin Y-axis values less than \code{ymin} are shown as
#'   \code{ymin}.
#' 
#' @param ymax Y-axis values greater than \code{ymax} are shown as
#' \code{ymax}. When \code{y = "z"}, setting \code{ymax} to a finite
#' value can improve the volcano plot when some z-scores are much
#' larger (in magnitude) than others.
#'
#' @param max.overlaps Argument passed to
#'   \code{\link[ggrepel]{geom_text_repel}}.
#'
#' @param plot.title The title of the plot.
#' 
#' @param ggplot_call The function used to create the plot. Replace
#'   \code{volcano_plot_ggplot_call} with your own function to customize
#'   the appearance of the plot.
#'
#' @return A \code{ggplot} object or a \code{plotly} object.
#'
#' @seealso \code{\link{de_analysis}}
#' 
#' @examples
#' # See help(de_analysis) for examples.
#' 
#' @export
#' 
volcano_plot <-
  function (de, k, labels, y = c("z","f0"),
            do.label = volcano_plot_do_label_default,
            ymin = 1e-6, ymax = Inf, max.overlaps = Inf,
            plot.title = paste("topic",k),
            ggplot_call = volcano_plot_ggplot_call) {
  if (!inherits(de,"topic_model_de_analysis"))
    stop("Input \"de\" should be an object of class ",
         "\"topic_model_de_analysis\"")
  y <- match.arg(y)
  n <- nrow(de$postmean)
  if (missing(labels)) {
    if (!is.null(rownames(de$postmean)))
      labels <- rownames(de$postmean)
    else
      labels <- as.character(seq(1,n))
  }
  if (!(is.character(labels) & length(labels) == n))
    stop("Input argument \"labels\", when specified, should be a character ",
         "vector with one entry per log-fold change estimate (column of ",
         "the counts matrix)")
  dat <- compile_volcano_plot_data(de,k,y,ymin,ymax,labels,do.label)
  return(ggplot_call(dat,y,plot.title,max.overlaps))
}

#' @rdname volcano_plot
#'
#' @param x An object of class \dQuote{topic_model_de_analysis},
#'   usually an output from \code{\link{de_analysis}}.
#'
#' @param \dots Additional arguments passed to \code{volcano_plot}.
#' 
#' @importFrom graphics plot
#' 
#' @method plot topic_model_de_analysis
#'
#' @export
#'
plot.topic_model_de_analysis <- function (x, ...)
  volcano_plot(x,...)

#' @rdname volcano_plot
#'
#' @param file Save the interactive volcano plot to this HTML
#'   file using \code{\link[htmlwidgets]{saveWidget}}.
#' 
#' @param width Width of the plot in pixels. Passed as argument
#'   \dQuote{width} to \code{\link[plotly]{plot_ly}}.
#'
#' @param height Height of the plot in pixels. Passed as argument
#'   \dQuote{height} to \code{\link[plotly]{plot_ly}}.
#'
#' @param plot_ly_call The function used to create the plot. Replace
#'   \code{volcano_plot_ly_call} with your own function to customize
#'   the appearance of the interactive plot.
#'
#' @importFrom htmlwidgets saveWidget
#' 
#' @export
#' 
volcano_plotly <- function (de, k, file, labels, y = c("z","f0"),
                            ymin = 1e-6, ymax = Inf, width = 500,
                            height = 500, plot.title = paste("topic",k),
                            plot_ly_call = volcano_plot_ly_call) {

  # Check and process input arguments.
  if (!inherits(de,"topic_model_de_analysis"))
    stop("Input \"de\" should be an object of class ",
         "\"topic_model_de_analysis\"")
  y <- match.arg(y)
  n <- nrow(de$postmean)
  if (missing(labels)) {
    if (!is.null(rownames(de$postmean)))
      labels <- rownames(de$postmean)
    else
      labels <- as.character(seq(1,n))
  }
  if (!(is.character(labels) & length(labels) == n))
    stop("Input argument \"labels\", when specified, should be a character ",
         "vector with one entry per log-fold change estimate (column of ",
         "the counts matrix)")
  dat <- compile_volcano_plot_data(de,k,y,ymin,ymax,labels)
  p <- volcano_plot_ly_call(dat,y,plot.title,width,height)
  if (!missing(file))
    saveWidget(p,file,selfcontained = TRUE,title = plot.title)
  return(p)
}

#' @rdname volcano_plot
#'
#' @param lfc A vector of log-fold change estimates.
#'
#' @param y A vector of the same length as \code{lfc}.
#'
#' @importFrom stats quantile
#' 
#' @export
#' 
volcano_plot_do_label_default <- function (lfc, y)
  y >= quantile(y,0.999,na.rm = TRUE) |
  lfc <= quantile(lfc,0.001,na.rm = TRUE) |
  lfc >= quantile(lfc,0.999,na.rm = TRUE)

#' @rdname volcano_plot
#'
#' @param dat A data frame passed as input to
#'   \code{\link[ggplot2]{ggplot}}, containing, at a minimum, columns
#'   \dQuote{f0}, \dQuote{postmean}, \dQuote{y}, \dQuote{lfsr} and
#'   \dQuote{label}.
#' 
#' @param font.size Font size used in plot.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 waiver
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
volcano_plot_ggplot_call <- function (dat, y, plot.title, max.overlaps = Inf,
                                      font.size = 9) {
    dat$lfsr <- cut(dat$lfsr,c(-1,0.001,0.01,0.05,Inf))
  if (y == "z")
    ybreaks <- c(0,1,2,5,10,20,50,100,200,500,1000,2000,5000,1e4,2e4,5e4)
  else
    ybreaks <- 10^seq(-8,0)
  return(ggplot(dat,aes_string(x = "postmean",y = "y",color = "lfsr",
                               label = "label")) +
         geom_point(shape = 20,na.rm = TRUE) +
         geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                         segment.color = "darkgray",segment.size = 0.25,
                         min.segment.length = 0,max.overlaps = max.overlaps,
                         na.rm = TRUE) +
         scale_y_continuous(trans = ifelse(y == "z","sqrt","log10"),
                            breaks = ybreaks) +
         scale_color_manual(values = c("deepskyblue","gold","orange","coral")) +
         labs(x = "log-fold change",
              y = ifelse(y == "z","|posterior z-score|","null estimate"),
                  title = plot.title) +
         theme_cowplot(font.size) +
         theme(plot.title = element_text(size = font.size,face = "plain")))
}

#' @rdname volcano_plot
#'
#' @importFrom plotly plot_ly
#' @importFrom plotly hide_colorbar
#' @importFrom plotly layout
#' 
#' @export
#' 
volcano_plot_ly_call <- function (dat, y, plot.title, width, height) {
  dat$fill <- cut(dat$lfsr,c(-1,0.001,0.01,0.05,Inf))
  if (y == "z")
    dat$y <- sqrt(dat$y)
  p <- plot_ly(data = dat,x = ~postmean,y = ~y,color = ~fill,
               colors = c("deepskyblue","gold","orange","coral"),
               text = ~sprintf(paste0("%s\nmean(null): %0.2e\n",
                                      "logFC(lower): %+0.3f\n",
                                      "logFC(mean): %+0.3f\n",
                                      "logFC(upper): %0.3f\n",
                                      "posterior z-score: %+0.3f\n",
                                      "lfsr: %0.2e"),
                               label,f0,lower,postmean,upper,z,lfsr),
               type = "scatter",mode = "markers",hoverinfo = "text",
               width = width,height = height)
  p <- hide_colorbar(p)
  p <- layout(p,
              xaxis = list(title = "log-fold change",zeroline = FALSE,
                           showgrid = FALSE),
              yaxis = list(title=ifelse(y == "z","sqrt|posterior z-score|",
                           "null estimate"),zeroline = FALSE,showgrid = FALSE,
                           type = ifelse(y == "z","identity","log")),
              hoverlabel = list(bgcolor = "white",bordercolor = "black",
                                font = list(color = "black",family = "arial",
                                            size = 12)),
              font = list(family = "arial",size = 12),
              showlegend = FALSE,title = plot.title)
  return(p)
}

# This is used by volcano_plot and volcano_plotly to compile the data
# frame passed to 'ggplot'.
compile_volcano_plot_data <- function (de, k, y, ymin, ymax, labels,
                                       do.label = NULL) {
  if (y == "z")
    y <- abs(de$z[,k])
  else
    y <- de$f0
  if (all(is.na(de$lfsr))) {
    message("lfsr is not available, probably because \"shrink.method\" was ",
            "not set to \"ash\"; lfsr in plot should be ignored")
    lfsr <- 0
  } else
    lfsr <- de$lfsr[,k]
  dat <- data.frame(label    = labels,
                    f0       = de$f0,
                    lower    = de$lower[,k],
                    postmean = de$postmean[,k],
                    upper    = de$upper[,k],
                    z        = de$z[,k],
                    y        = clamp(y,ymin,ymax),
                    lfsr     = lfsr,
                    stringsAsFactors = FALSE)
  dat <- dat[order(dat$lfsr,decreasing = TRUE,na.last = FALSE),]
  if (!is.null(do.label))
    dat$label[which(!do.label(dat$postmean,dat$y))] <- ""
  return(dat)
}
