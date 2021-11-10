#' @title splitFragPipeProteinTable
#' Splits FragPipe combined_protein.tsv table into one table per quantification type
#' @param combined_protein_file path to combined folder location. Character, either of length 1; light crosslinker only is assumed;
#' if two folders are provided (length 2) it is assumed that an isotopically labelled XL reagent was used; with two separate searches performed
#' for each channel. The first folder is assumed to be from the light search and the second folder is assumed to be from the heavy search.
#' @param write_split_tables Logical, whether split tables should be written, or only returned as elements of a list in R. Default: TRUE
#' @return a list of the separate tables
#'
#' @import data.table
#' @export

splitFragPipeProteinTable <- function(combined_protein_tsv= "combined_protein.tsv",
  write_split_tables = TRUE){

  # read combined folder
  combined_protein = fread(combined_protein_tsv)
  names(combined_protein) = gsub(" ", "_", names(combined_protein))

  # isolate protein description columns
  # Note the 'indistinguishable proteins' column is at the very end of the table
  combined_protein_descr = combined_protein[, c(1:17,ncol(combined_protein))]

  # isolate quants
  combined_protein_quant = combined_protein[, c(18:ncol(combined_protein)-1)]

  quanttypes = c("Total_Spectral_Count", "Razor_Spectral_Count", "Unique_Ion_Count", "Total_Intensity",
               "Razor_Intensity", "Unique_Spectral_Count", "Total_Ion_Count", "Razor_Ion_Count", "Unique_Intensity")

  # Generate result container list
  res = list()

  # iterate over quanttypes and assemble separate tables
  for (i in quanttypes){
    combined_protein_quant_i =
      combined_protein_quant[, grep(i, names(combined_protein_quant)), with = F]
    dt_i = cbind(combined_protein_descr, combined_protein_quant_i)
    if(write_split_tables){
      fwrite(dt_i, file = paste0("combined_protein_", i,".csv"))
    }
    res[[i]] = dt_i
  }

  # return result list
  return(res)

}
