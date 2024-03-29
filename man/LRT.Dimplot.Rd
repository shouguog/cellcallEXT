% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sankey3D.R
\name{LRT.Dimplot}
\alias{LRT.Dimplot}
\title{get score of Triple relation corresponding specific cellA-cellB}
\usage{
LRT.Dimplot(
  object,
  fontSize = 12,
  nodeWidth = 30,
  nodePadding = 5,
  height = NULL,
  width = 1200,
  sinksRight = FALSE,
  DIY.color = FALSE,
  color.DIY = NULL
)
}
\arguments{
\item{object}{a Cellwave objects or a dataframe about score of Triple relation corresponding specific cellA-cellB}

\item{fontSize}{the font size of text.}

\item{nodeWidth}{the node Width of sankey graph.}

\item{nodePadding}{the padding of node.}

\item{height}{the height of graph, default is NULL.}

\item{width}{the width of graph, default is 1200.}

\item{sinksRight}{boolean. If TRUE, the last nodes are moved to the right border of the plot.}

\item{DIY.color}{boolean. If TRUE, set the parameter color.DIY with your color-setting, default is FALSE.}

\item{color.DIY}{a color dataframe, rownames is cell type, value is color, default is NULL.}
}
\description{
get score of Triple relation corresponding specific cellA-cellB
}
