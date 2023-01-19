#' IDOverlapHeatmap
#' @description Calculates and visualizes in heatmap the % overlap (intersect/union) of ID sets from a proteomics matrix
#' @param proteomicsMatrix Quantitative proteomics matrix with analyte ID in row.names and measurement id as column names.
#' Either NAs or 0s for missing values
#' @return A distance matrix with ID set overlap in percent
#' @export

IDOverlapHeatmap <- function(proteomicsMatrix, save_pdf = TRUE, pdf_suffix =  "objname", size = 7){
  IDOverlapMatrix = calculateIDOverlapMatrix(proteomicsMatrix)

  p = corrplot.mixed(round(m.log10.IDsim*100,0), is.corr = F, tl.pos = 'l',
                     lower = "number",
                     upper = "color",
                     tl.col = 'black',
                     lower.col = 'black',
                     number.digits = 0,
                     number.cex = 0.7)
  print(p)

  if(save_pdf){
    if (pdf_suffix == "objname"){
      pdf_suffix_used = deparse(substitute(proteomicsMatrix))
    } else{
      pdf_suffix_used = pdf_suffix
    }
  }
  dev.copy(pdf, file = paste0("IDOverlapHeatmap_",pdf_suffix_used, ".pdf"), height = size, width = size)
  dev.off

  return(list("IDOverlapMatrix" = IDOverlapMatrix,
              'IDOverlapHeatmap' = p))
}
