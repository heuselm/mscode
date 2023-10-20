#' alignPepModFormat
#' @description Align modification format from Mascot [Oxidation] to UniMod (UniMod:35) format
#' @param Modpeps_in Vector of Modified peptide sequences in mascot [Modname] notation format
#' @return Vector of Modified peptide sequences with (UniMod:x) notation format
#' @import data.table
#' @export

alignPepModFormat <- function(Modpeps_in = c("_M[Oxidation]ADEAGSEADHEGTHSTK_", "_KDGEAGAQ[Deamidated]GPP[Oxidation]GP_"))
  {
  # copy
  Modpeps_out = copy(Modpeps_in)# use name of sets_list object in filenames

  # replace most frequent mod syntax PD/Mascot --> UniMod
  Modpeps_out = gsub("\\[Acetylated\\]","\\(UniMod:1\\)", Modpeps_out)
  Modpeps_out = gsub("\\[Acetylation\\]","\\(UniMod:1\\)", Modpeps_out)
  Modpeps_out = gsub("\\[Acetyl\\]","\\(UniMod:1\\)", Modpeps_out)
  Modpeps_out = gsub("\\[Amidated \\(Any C-term\\)\\]","\\(UniMod:2\\)", Modpeps_out)
  Modpeps_out = gsub("\\[Amidated\\]","\\(UniMod:2\\)", Modpeps_out)
  Modpeps_out = gsub("\\[Carbamidomethyl\\]","\\(UniMod:4\\)", Modpeps_out)
  Modpeps_out = gsub("\\[Carbamylation\\]","\\(UniMod:5\\)", Modpeps_out)
  Modpeps_out = gsub("\\[Deamidated\\]","\\(UniMod:7\\)", Modpeps_out)
  Modpeps_out = gsub("\\[Phosphorylation\\]","\\(UniMod:21\\)", Modpeps_out)
  Modpeps_out = gsub("\\[Phospho\\]","\\(UniMod:21\\)", Modpeps_out)
  Modpeps_out = gsub("\\[Oxidation\\]","\\(UniMod:35\\)", Modpeps_out)

  return(Modpeps_out)
}
