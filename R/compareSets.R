#' compareSets
#' @description Comparison of sets, e.g. peptides or proteins detected by a given analysis.
#' Draws both UpSetR and roughly area-proportional Venn Diagrams via the eulerr package.
#' For the comparison of more sets, see compareSets function.
#' @param sets_list List of Character vectors with the respective sets. List names must be unique!
#' @param plot_pdf selection whether pdf plots are created or whether the output occurs exclusively within Rconsole.
#' @param venn_quantities whether quantities (overlap size numbers) are plotted in Venn Diagrams. This is
#' recommended because the eulerr Venn Diagrams overlap areas are only with error proportional to actual overlap size.
#' @param plot_title character selection either "R" or "pdf"
#' @return returns a list of set comparison objects and plots from both UpSetR and Eulerr
#' @import eulerr
#' @import UpSetR
#' @import grid
#' @import gridExtra
#' @export

compareSets <- function(sets_list = list("Set1" = c("A","B","C","D","E"),
                                         "Set2" = c("A","B","D","E","F"),
                                         "Set3" = c("B","D","E","F","G","H")),
                        venn_quantities = TRUE,
                        plot_pdf = TRUE,
                        pdf_analysis_tag = "current_analysis",
                        plot_title = "set comparison")
  {
  # use name of sets_list object in filenames
  list_obj_name = deparse(substitute(sets_list))

  # Get set analysis intersect groups
  detection_sets = getSetAnalysisGroups(sets_list)
  barplot(sapply(detection_sets, length), main = "intersect set sizes")

  if(plot_pdf){
    dev.copy(pdf, file = paste0("compareSets_intersectsetsizes",list_obj_name,".pdf"),
             height = 5+length(sets_list),
             width = 5+length(sets_list))
    dev.off()
  }

  # convert to detection matrix
  detection_matrix = fromListWithNames(sets_list)

  # UpSetR plot
  upset_obj = UpSetR::upset(detection_matrix)
  print(upset_obj)
  grid::grid.text(paste(plot_title,"\nUpSet plot"), x = 0.13, y = 0.95,  gp=gpar(fontsize=20))

  # plot pdf if desired
  if(plot_pdf){
    dev.copy(pdf, file = paste0("compareSets_upset",list_obj_name,".pdf"),
                                height = 5+length(sets_list),
                                width = 5+length(sets_list))
    dev.off()
  }

  # Venn Diagram (Euler)
  euler_obj = eulerr::euler(detection_matrix)
  plot_venn = plot(euler_obj, quantities = venn_quantities, main = plot_title)
  plot(plot_venn)

  # plot pdf if desired
  if(plot_pdf){
    dev.copy(pdf, file = paste0("compareSets_Venn",list_obj_name,".pdf"),
             height = 5+length(sets_list),
             width = 5+length(sets_list))
    dev.off()
  }

  # Assemble results
  res = list("detection_matrix" = detection_matrix,
             "detection_sets" = detection_sets,
             "plot_UpSet" = upset_obj,
             "plot_Venn" = plot_venn,
             "eulerr_object" = euler_obj)
  return(res)
}
