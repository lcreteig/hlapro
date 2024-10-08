% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsa.R
\name{read_lum_csv}
\alias{read_lum_csv}
\title{Parse a Luminex single-antigen csv file into a data frame}
\usage{
read_lum_csv(csv_filepath, lots_path)
}
\arguments{
\item{csv_filepath}{A character path to the csv file.}

\item{lots_path}{A character path to the folder that stores the lot-specific
file (\code{.eds}) for the kit that was used to run the single antigen bead assay.
If the folder contains multiple lot files, the right one is read in
automatically (based on the header of the csv file).}
}
\value{
A data frame with (a selection of) the contents from the \code{.csv} and
\code{.eds} file. Its contents match the table that would be produced by loading
these files into Immucor's MATCH IT!® Antibody Analysis Software.
}
\description{
\code{read_lum_csv()} reads in a csv with raw Luminex results from a single
antigen bead assay, together with its associated lot file from Immucor, and
returns all the information therein in a single data frame.
}
\examples{
\dontrun{
dat <- read_lum_csv("~/LSA1-001.csv", "~/lot_files/")
}
}
