% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperFunctions.R
\name{fromListWithNames}
\alias{fromListWithNames}
\title{fromListWithNames}
\usage{
fromListWithNames(
  input = list(Set1 = c("A", "B", "C", "D", "E"), Set2 = c("A", "B", "D", "E", "F"), Set3
    = c("B", "D", "E", "F", "G", "H"))
)
}
\arguments{
\item{input}{list of named character vectors containing the sets}
}
\description{
alternative to UpSetR::fromList that maintains element names in row.names of the output detection matrix.
Required to track members in intersect sets. Grabbed & modified from https://github.com/hms-dbmi/UpSetR/issues/85.
}
