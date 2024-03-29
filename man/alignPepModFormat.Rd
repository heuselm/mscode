% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignModPepPtmFormat.R
\name{alignPepModFormat}
\alias{alignPepModFormat}
\title{alignPepModFormat}
\usage{
alignPepModFormat(
  Modpeps_in = c("_M[Oxidation]ADEAGSEADHEGTHSTK_",
    "_KDGEAGAQ[Deamidated]GPP[Oxidation]GP_")
)
}
\arguments{
\item{Modpeps_in}{Vector of Modified peptide sequences in mascot [Modname] notation format}
}
\value{
Vector of Modified peptide sequences with (UniMod:x) notation format
}
\description{
Align modification format from Mascot [Oxidation] to UniMod (UniMod:35) format
}
