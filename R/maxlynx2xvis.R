#' @title maxlynx2xvis
#' simplify maxlynx crosslinks table for visualization of inter- and intra-protein crosslinks via xvis
#' (https://xvis.genzentrum.lmu.de/CrossVisNoLogin.php)
#' @param combined_folder_location path to combined folder location. Character, either of length 1; light crosslinker only is assumed;
#' if two folders are provided (length 2) it is assumed that an isotopically labelled XL reagent was used; with two separate searches performed
#' for each channel. The first folder is assumed to be from the light search and the second folder is assumed to be from the heavy search.
#' @param expDesign In case there is no filled-in experimentalDesignTemplate.txt in the
#' combined folder, here, a data.table with column names 'Name' and 'Experiment'
#' should minimally be provided. Additional columns will be propagated to the output xvis_xls.csv
#' @param xvis_filter Filter string which Crosslink types should be included in the output xvis_xls.csv; defaults to "InterProXL|IntraProXL" (both included)
#' @param xvis_remove_decoys Logical, whether t-d, d-t and d-d XLs shall be removed before export to xvis_xls.csv. Defaults to TRUE.
#' @param plot_pdf Logical, whether overview pdf should be generated, highlighting the ID numbers of and overlap between channels
#' @return a list of 1) merged data.table of crosslinks underlying xvis_xls.csv; 2) Crosslink MSMS count table underlying summary plot;
#' 3) Summary ggplot
#' @export


maxlynx2xvis = function(combined_folder_location = "path/to/combined/",
                        expDesign = NULL,
                        xvis_filter = "InterProXL|IntraProXL",
                        xvis_remove_decoys = TRUE,
                        plot_pdf = TRUE){

  # Read XLMSMS.txt files
  xls_light = fread(paste0(combined_folder_location[1],"/txt/crosslinkMsms.txt"))
  prot_light = fread(paste0(combined_folder_location[1],"/txt/proteinGroups.txt"))
  if (length(combined_folder_location) == 2){
    xls_heavy = fread(paste0(combined_folder_location[2],"/txt/crosslinkMsms.txt"))
    prot_heavy = fread(paste0(combined_folder_location[2],"/txt/proteinGroups.txt"))
  } else {
    error("combined_folder_location parameter must be provided as character of length 1-2")
  }

  # Read experimental Design
  if (is.null(expDesign)){
    expDes = fread(paste0(combined_folder_location[1],"experimentalDesignTemplate.txt"))
  } else if (class(expDesign) == "character"){
    expDes = fread(expDesign)
  } else if (class(expDesign) == "data.frame"){
    expDes = as.data.table(expDesign)
  } else {
    expDes = expDesign
  }

  # Combine
  if (length(combined_folder_location) == 1){
    xls = xls_light[, label:="light"]
    prot = prot_light[, label:="light"]
  } else {
    xls = rbind(xls_light[, label:="light"],
    xls_heavy[, label:='heavy'])

    prot = rbind(prot_light[, label:="light"],
    prot_heavy[, label:='heavy'])

  }
  names(xls) = gsub(" ", "_", names(xls))
  names(prot) = gsub(" ", "_", names(prot))

  # Annotate
  xls = merge(xls, expDes, by.x = "Raw_file", by.y = "Name")

  # Count N and plot
  xls_n = unique(xls[, .N, .(Raw_file, label, Crosslink_Product_Type, Experiment)])

  p1 = ggplot(xls_n) + geom_bar(aes(Experiment, y = N, fill = label), stat = "identity", position = "dodge") +
  facet_wrap(~Crosslink_Product_Type, scales = "free_y") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("MaxLynx XL-MS result overview")
  plot(p1)

  if(plot_pdf){
    ggsave(device = "pdf", filename = "MaxLynx_Results.pdf", plot = p1, width = 7, height = 7)
  }

  # Add id column
  xls[, Id:=paste(Sequence1, Sequence2,
    paste0("a",gsub(";","", Pep_InterLink1)),
    paste0("a",gsub(";","", Pep_InterLink1)),
    sep = "-"), row.names(xls)]

  # Assemble & output xl tables for xvis:
  xls[, AbsPos1:=as.numeric(gsub(";","", Pro_InterLink1)), row.names(xls)]
  xls[, AbsPos2:=as.numeric(gsub(";","", Pro_InterLink2)), row.names(xls)]

  # Subset xls table for xvis
  xvis_xls = unique(xls[grep(xvis_filter, Crosslink_Product_Type)])

  # remove decoys
  if (xvis_remove_decoys == TRUE){
    xvis_xls = xvis_xls[grep("REV_", Proteins1, invert = T)][grep("REV_", Proteins2, invert = T)]
  }
  setnames(xvis_xls, "Proteins1", "Protein1")
  setnames(xvis_xls, "Proteins2", "Protein2")
  setnames(xvis_xls, "Mass_error_[ppm]", "Mass_error_ppm")
  fwrite(xvis_xls, "xvis_xls.csv")

  # assemble Protein lengths
  xvis_lengths = unique(data.table("Protein" = unique(c(xvis_xls$Protein1, xvis_xls$Protein2))))
  xvis_lengths = merge(xvis_lengths, unique(prot[, .(Protein_IDs, Sequence_length)]), by.x = "Protein", by.y = "Protein_IDs", all.x = T, all.y = F)
  setnames(xvis_lengths, "Sequence_length", "Length")
  fwrite(xvis_lengths, "xvis_lengths.csv")

  res = list("xls" = xls,
              "xls_n" = xls_n,
              "summary_ggplot" = plot)
  return(res)
}

