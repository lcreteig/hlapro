% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downscaling.R
\name{get_serology}
\alias{get_serology}
\title{Retrieve serological equivalents of an HLA-allele}
\usage{
get_serology(allele)
}
\arguments{
\item{allele}{A string or character vector with (an) HLA allele(s).}
}
\value{
A string or character vector of the same length as \code{allele},
with the corresponding serology if it exists, or \code{NA} if none exists.
}
\description{
\code{get_serology()} takes in a string or character vector of HLA alleles. The
corresponding split-level (if it exists) or broad-level (if no split exists)
allele is looked up in \link{etrl_hla} and returned. If no such alleles exist,
\code{NA} is returned instead.
}
\details{
This function uses the EuroTransplant Reference Laboratory HLA (\link{etrl_hla})
tables to do the lookup. These tables define several alleles for each protein
(e.g. \code{HLA-A*01:01}, \code{HLA-A*01:02}), and their serological equivalents. All
others are grouped into an XX code (\code{HLA-A*01:XX}), which is also mapped to
a serological equivalent.
\subsection{Workings}{

All entered alleles will first be reduced to the two-field level (see
\code{\link[=reduce_to_nth_field]{reduce_to_nth_field()}}). If this reduced allele occurs in the lookup table,
the corresponding rows are returned. If not, it will be converted into an XX
code, and that row will be returned.

If the allele is already at the serological broad- or split-level, the lookup
will be performed using those respective columns in \link{etrl_hla}.
}

\subsection{Exceptions}{

If the allele has a suffix (e.g. \verb{HLA-C*01:37N}), it has no serological
equivalent, and hence will also return nothing.
}
}
\examples{
get_serology("A24") # is a serological split; returns itself ("A24")
get_serology("A9") # is a serological broad; returns itself ("A9")

get_serology("A*01:01:01:50") # has a broad-level equivalent only ("A1")
get_serology("A*23:01:01:11") # has a split equivalent ("A23")
# as well as a broad ("A9"); only the former is returned

# Vectors also work:
get_serology(c("A24", "A*01:XX", "B*15:15"))
}
\seealso{
\itemize{
\item \link{get_broad}: for looking up the broad-level equivalent of a split allele
\item \link{get_split}: for looking up the serological split-level equivalent of an
allele
\item \link{get_public}: for looking up the public epitope of an allele
\item \link{etrl_hla}: the lookup table that's used in this function
}
}
