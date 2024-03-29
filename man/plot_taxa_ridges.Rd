% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_taxa_ridges.R
\name{plot_taxa_ridges}
\alias{plot_taxa_ridges}
\title{Taxon change point ridge plots}
\usage{
plot_taxa_ridges(
  titan.out,
  z1 = TRUE,
  z2 = TRUE,
  z1_fill_low = "lightblue",
  z1_fill_high = "black",
  z2_fill_low = "pink",
  z2_fill_high = "red",
  xlabel = "Environmental Gradient",
  n_ytaxa = 90,
  printspp = FALSE,
  grid = TRUE,
  d_lines = FALSE,
  trans = "identity",
  xaxis = !grid,
  xlim,
  breaks = ggplot2::waiver(),
  bw,
  rel_heights,
  ...,
  axis.text.x,
  axis.text.y,
  axis.title.x,
  axis.title.y
)
}
\arguments{
\item{titan.out}{A TITAN output object.}

\item{z1}{A logical specifying whether decreasing taxa should be plotted.}

\item{z2}{A logical specifying whether increasing taxa should be plotted.}

\item{z1_fill_low, z1_fill_high, z2_fill_low, z2_fill_high}{Respective fill
colors passed to \code{\link[=scale_fill_gradient]{scale_fill_gradient()}}}

\item{xlabel}{A character string for the x axis label.}

\item{n_ytaxa}{The maximum number of taxa to be plotted.}

\item{printspp}{A logical specifying whether the sppmax table should be
printed.}

\item{grid}{The \code{grid} argument of \code{\link[=theme_ridges]{theme_ridges()}}. Setting this to
\code{FALSE} removes horizontal lines from the plot, see examples.}

\item{d_lines}{Argument to pass to \code{\link[=geom_density_ridges]{geom_density_ridges()}}.
Short for density lines, a logical determining whether the median change
point values from estimated from each taxon's density function, is plotted.
If \code{TRUE}, a dashed, black vertical line will be plotted.}

\item{trans}{a scale transformation to be applied to the x-axis through
ggplot2. e.g. \code{"log10"} or \code{"sqrt"}.}

\item{xaxis}{Logical; should the x-axis be plotted?}

\item{xlim}{x axis limits, e.g. \code{xlim = c(0,10)}.}

\item{breaks}{Argument to pass to \code{\link[=geom_density_ridges]{geom_density_ridges()}}.
Determines values and labels of breaks (ticks) on the x-axis.}

\item{bw}{The bandwidth of used in the kernel density estimate; see
\code{\link[=density]{density()}}.}

\item{rel_heights}{Argument to pass to \code{\link[cowplot:plot_grid]{cowplot::plot_grid()}} as an override
for the defaults}

\item{...}{Arguments to pass to \code{\link[=geom_density_ridges]{geom_density_ridges()}}}

\item{axis.text.x, axis.text.y, axis.title.x, axis.title.y}{Font sizes of
respective axis text. Passed as the size argument to the ggplot2 thematic
elements of the same name. Defaults to current default ggplot2 theme.}
}
\value{
A plot of decreasing and/or increasing taxon-specific change points
along the environmental gradient.
}
\description{
Taxon change point ridge plots
}
\note{
Should not be used with output objects from TITAN v1.0.
}
\examples{

data(glades.titan)


# basic usage
plot_taxa_ridges(glades.titan)


\dontrun{ reduce R CMD check time

# manipulating the x axis thematic components
plot_taxa_ridges(glades.titan, grid = FALSE)
plot_taxa_ridges(glades.titan, xaxis = TRUE)


# applying scale transformations
plot_taxa_ridges(glades.titan, trans = "log10", xlim = c(1,150))


# removing z2
plot_taxa_ridges(glades.titan, z2 = FALSE, trans = "sqrt", xlim = c(0, 150))


# more styling
plot_taxa_ridges(glades.titan,
  axis.text.x = 12,
  axis.text.y = 4,
  axis.title.x = 14,
  xlabel = expression(paste("Surface water total phosphorus ("*mu*"g/l)"))
)

}



}
\references{
Baker, ME and RS King.  2010. A new method for detecting and
interpreting biodiversity and ecological community thresholds.  Methods in
Ecology and Evolution 1(1): 25:37.

King, RS and ME Baker  2010. Considerations for identifying and
interpreting ecological community thresholds. Journal of the North American
Benthological Association 29(3):998-1008.
}
\seealso{
\code{\link[=plot_sumz]{plot_sumz()}}, \code{\link[=plot_cps]{plot_cps()}}
}
\author{
M. Baker, R. King, D. Kahle
}
