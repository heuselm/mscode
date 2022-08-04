#' convertMzID2tsv
#' @description Convert mzIdentML files to human-readable .tsv via the R/mzID package.
#' @param mzidfile path to .mzid file to be converted
#' @return Writes mzidfile.tsv next to input file and returns a list of mzID R object and the corresponding flattened data.frame
#' @import mzID
#' @export

convertMzID2tsv = function(mzidfile){
  filemz = mzID(mzidfile)
  # assign(paste0(basename(mzidfile),"_mzID"), filemz)
  filemz_df = flatten(filemz)
  # assign(paste0(basename(mzidfile),"_df"), filemz_df)
  write.table(filemz_df, file = paste0(mzidfile, ".tsv"), sep = "\t", row.names = F)
  return(list("mzID_object" = filemz,
              "mzID_flat_dataframe" = filemz_df))
}
