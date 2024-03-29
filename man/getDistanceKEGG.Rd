% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codeAllCaseControl_v1.5.R
\name{getDistanceKEGG}
\alias{getDistanceKEGG}
\title{abstract path from kegg pathway file}
\usage{
getDistanceKEGG(data, method)
}
\arguments{
\item{data}{a dataframe store L-R-TF relation in KEGG and the distance between LR and TF}

\item{method}{the method to calculate multiple path corresponding specific L-R-TF, default is "mean"}
}
\value{
the distance between LR and TF
}
\description{
abstract path from kegg pathway file
}
