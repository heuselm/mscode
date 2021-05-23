#' compareDiannResultSets
#' @description Comparison of DIANN resultsets stored side-by-side in a result folder,
#' e.g. Results from pass1 analysis and pass2 re-analysis (Match Between Runs)
#' or results obtained using different input spectral libraries. Compares
#' the ID numbers from stats.tsv, but also the ID numbers obtained on protein group level
#' by counting from the quantitative tables. Compares sets on identified proteins in Venn Diagrams
#' and their quantitative patterns in side-by-side heatmap. Further, pairwise sample-sample correlation
#' for a subset of runs is implemented, however currently not experimental-design-aware.
#' @param result_folder Path to the result folder. All resultsets identified by their *.stats.tsv file will be analyzed.
#' @param Venn_diagram Whether to calculate and draw Venn diagrams from the sets of analytes detected
#' @param compare_precursor_level Whether to include precursor level comparisons in Venn and heatmap. This might
#' be slow for calculation of heatmap, defaults to FALSE.
#' @param heatmap Whether to calculate and draw heatmaps, defaults to true. These are the most informative plots..
#' @param sample_correlation Whether to calculate sample-sample correlation and output heatmap and xy scatters on subset
#' @return Returns a combined data.table of merged quantitative results with extra column resultset
#' @import eulerr
#' @import ggplot2
#' @import pheatmap
#' @import corrplot
#' @import psych
#' @export

compareDiannResultSets = function(result_folder = "./",
                                  Venn_diagram = TRUE,
                                  compare_precursor_level = FALSE,
                                  heatmap = TRUE,
                                  sample_correlation = TRUE){

  # Not a great coder, forgive me (for now) ;-)
  prev_wd = getwd()
  setwd(result_folder)

  ## Load results:
  list.files()[grep(".stats.tsv", list.files())]
  n_resultsets = length(list.files()[grep(".stats.tsv", list.files())])
  stats = data.table()
  quant = data.table()

  for (file in list.files()[grep(".stats.tsv", list.files())]){
    message("loading ", file)
    stats_i = fread(file)
    qres_i = gsub(".stats", "", file)
    message("loading ", gsub(".stats", "", file))
    quant_i = fread(gsub(".stats", "", file))
    stats = rbind(stats, stats_i[, resultset:=gsub(".tsv", "", qres_i)], fill = TRUE)
    quant = rbind(quant, quant_i[, resultset:=gsub(".tsv", "", qres_i)], fill = TRUE)
  }

  n_runs = length(unique(stats$File.Name))

  # # Annotation table
  # fread("../../Annotation.txt")
  # annotation = fread("../../Annotation.txt")

  #############################################################################################
  # Compare Stats (i.e. re-produce plots similar to DIANN report with res-sets side-by-side) ##
  #############################################################################################
  # Raw IDs from stats.tsv files
  stats[, File.Path:=File.Name]
  stats[, File.Name:=unlist(strsplit(File.Name, split = "\\\\"))
        [length(unlist(strsplit(File.Name, split = "\\\\")))], File.Name]

  pdf("01_Stats_comparison.pdf", width = (n_runs/2)+3)

  for (i in names(stats)[2:length(names(stats))]){
    p = ggplot(stats, aes_string("File.Name", y = i, fill = "resultset")) +
      geom_bar(stat = "identity", position = "dodge") +
      ggtitle(paste(i, "per DIANN result set (stats)")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    plot(p)
  }

  dev.off()

  ggplot(stats, aes(File.Name, y = Precursors.Identified, fill = resultset)) +
    geom_bar(stat = "identity", position = "dodge") +
    ggtitle("Precursors.Identified per DIANN result set (stats)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave("01_Stats_Precursors.Identified.pdf", width = 12, height = 6)

  ggplot(stats, aes(File.Name, y = Proteins.Identified, fill = resultset)) +
    geom_bar(stat = "identity", position = "dodge") +
    ggtitle("Proteins.Identified per DIANN result set (stats)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave("01_Stats_Proteins.Identified.pdf", width = 12, height = 6)

  # No, write a combined table with these values side-by side
  stats_melt = melt(stats, id.vars = c("File.Name", "resultset"))
  stats_recast = dcast(stats_melt, File.Name~variable+resultset)
  fwrite(stats_recast, "01_Stats_comparison_wide.csv")

  ##########################################################################################
  # Compare Quantitative Data (including protein group counts that aren't assessed above) ##
  ##########################################################################################
  # Raw IDs from stats.tsv files

  # flatten filenames
  # quant[grep("sp\\|", Protein.Group), Protein.Group:=unlist(strsplit(Protein.Group, split = "\\|"))[2], .(Protein.Group)]
  # quant[grep("sp\\|", Protein.Ids), Protein.Ids:=unlist(strsplit(Protein.Ids, split = "\\|"))[2], .(Protein.Ids)]

  quant[, File.Path:=File.Name]
  quant[, File.Name:=unlist(strsplit(File.Name, split = "\\\\"))[length(unlist(strsplit(File.Name, split = "\\\\")))], File.Name]

  # plot ids per run per resultset
  ggplot(quant[, length(unique(Precursor.Id)), .(resultset, File.Name)], aes(File.Name, y = V1, fill = resultset)) +
    geom_bar(stat = "identity", position = "dodge") + ylab("N.Precursors") +
    ggtitle("Precursors.Identified per DIANN result set (quantdata)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave("02_Identifications_Precursors_perRun_quantdata.pdf", width = 12, height = 6)

  ggplot(quant[, length(unique(Protein.Group)), .(resultset, File.Name)], aes(File.Name, y = V1, fill = resultset)) +
    geom_bar(stat = "identity", position = "dodge") + ylab("N.Protein.Groups") +
    ggtitle("Protein.Groups.Identified per DIANN result set (quantdata)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave("02_Identifications_ProteinGroups_perRun_quantdata.pdf", width = 12, height = 6)

  # Venn of IDs, per result set raw
  if(Venn_diagram){
    IDs_raw = list()
    for (i in unique(quant$resultset)){
      prots.i = quant[resultset == i, unique(Protein.Group)]
      IDs_raw[[i]] = prots.i
    }
    str(IDs_raw)
    p = plot(euler(IDs_raw), quantities = TRUE, main = "Protein Groups IDd per resultset", adjust_labels = T)
    plot (p)
    dev.copy2pdf(file = "03_Venn_ProteinGroup_Sets_per_resultset.pdf", height = 7, width = 7)
  }

  if(compare_precursor_level){
    IDs_raw_prec = list()
    for (i in unique(quant$resultset)){
      prec.i = quant[resultset == i, unique(Precursor.Id)]
      IDs_raw_prec[[i]] = prec.i
    }
    str(IDs_raw_prec)
    plot(euler(IDs_raw_prec), quantities = TRUE, main = "Precursor IDd per resultset", adjust_labels = T)
    dev.copy2pdf(file = "03_Venn_Precursor_Sets_per_resultset.pdf", height = 7, width = 7)
  }

  # intensity distributions
  ggplot(quant) + geom_density(aes(log10(Precursor.Quantity), y = ..count.., col = resultset), alpha = 0.3) +
    theme_bw() + ggtitle("N-scaled Precursor Intensity distributions")
  ggsave("04_IntensityDistGlobal.pdf", width = 5, height = 4)

  ggplot(quant) + geom_violin(aes(x = File.Name, y = log10(Precursor.Quantity), fill = resultset), alpha = 0.3) +
    theme_bw() + ggtitle("Precursor Intensity distributions") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave("04_IntensityDistPerRun.pdf", width = 20, height = 10)

  # count peptides per protein
  quant[, peptides_per_protein:=length(unique(Modified.Sequence)), .(Protein.Group, File.Name, resultset)]

  ggplot(quant[, length(unique(Modified.Sequence)), .(Protein.Group,resultset)], aes(V1, fill = resultset)) +
    geom_histogram(position = position_dodge(width = 0.7, preserve = "total"), breaks = c(c(1:30), 1000)) + ylab("N") +
    coord_cartesian(xlim = c(1,30)) + xlab("Peptides (Modified.Sequence) per Protein.Group") +
    ggtitle("Peptides per protein per DIANN result set (quantdata)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_continuous(breaks = c(1:30))
  ggsave("05_Peptides_per_protein.pdf", width = 7, height = 7)

  #############################
  # Comparative heatmaps ######
  #############################

  ###################
  ## protein level ##
  ###################
  quant.prot = dcast(unique(quant[,.(resultset,File.Name,PG.Quantity, Protein.Group)]),
                     Protein.Group~paste(resultset,File.Name), value.var = "PG.Quantity")

  quant.prot.m = as.matrix(quant.prot[,2:ncol(quant.prot)])
  row.names(quant.prot.m) = quant.prot$Protein.Group

  quant.prot.m.log10 = log10(quant.prot.m)
  quant.prot.m.log10.clean = copy(quant.prot.m.log10)
  quant.prot.m.log10.clean[is.infinite(quant.prot.m.log10.clean)] = 0
  quant.prot.m.log10.clean[is.na(quant.prot.m.log10.clean)] = 0
  hist(rowSums(quant.prot.m.log10.clean))

  ## Visualize quantitative data in heatmap, protein level
  # prep annotation
  colnames(quant.prot.m.log10.clean)
  ann_col = data.frame("rn" = colnames(quant.prot.m.log10.clean))
  ann_col$resultset = sapply(ann_col$rn, function(x){strsplit(x, split = " ")[[1]][1]})
  ann_col$File.Name = sapply(ann_col$rn, function(x){strsplit(x, split = " ")[[1]][2]})
  row.names(ann_col) = ann_col$rn
  ann_col$rn = NULL

  # plot heatmap
  n_runs = length(unique(quant$File.Name))
  gap_pos = n_runs * 1:(n_resultsets-1)

  p = pheatmap(quant.prot.m.log10.clean,
               show_rownames = F,
               annotation_col = ann_col,
               cluster_cols = FALSE,
               cluster_rows = TRUE,
               gaps_col = gap_pos)

  print(p)

  dev.copy(pdf, "06_Heatmap_proteinlevel_log10.pdf", height = 24, width = n_resultsets*8)
  dev.off()

  #####################
  ## precursor level ##
  #####################
  if(compare_precursor_level){
    quant.prec = dcast(unique(quant[,.(resultset,File.Name,Precursor.Quantity, Protein.Group, Precursor.Id)]),
                       Protein.Group~paste(resultset,File.Name), value.var = "Precursor.Quantity", fun.aggregate = sum)

    quant.prec.m = as.matrix(quant.prec[,2:ncol(quant.prec)])
    row.names(quant.prec.m) = quant.prec$Protein.Group

    quant.prec.m.log10 = log10(quant.prec.m)
    quant.prec.m.log10.clean = copy(quant.prec.m.log10)
    quant.prec.m.log10.clean[is.infinite(quant.prec.m.log10.clean)] = 0
    quant.prec.m.log10.clean[is.na(quant.prec.m.log10.clean)] = 0
    hist(rowSums(quant.prec.m.log10.clean))

    ## Visualize quantitative data in heatmap, protein level
    # prep annotation
    colnames(quant.prec.m.log10.clean)
    ann_col = data.frame("rn" = colnames(quant.prec.m.log10.clean))
    ann_col$resultset = sapply(ann_col$rn, function(x){strsplit(x, split = " ")[[1]][1]})
    ann_col$File.Name = sapply(ann_col$rn, function(x){strsplit(x, split = " ")[[1]][2]})
    row.names(ann_col) = ann_col$rn
    ann_col$rn = NULL

    # plot heatmap
    gap_pos = n_runs * 1:(n_resultsets-1)

    p = pheatmap(quant.prec.m.log10.clean,
                 show_rownames = F,
                 annotation_col = ann_col,
                 cluster_cols = FALSE,
                 cluster_rows = TRUE,
                 gaps_col = gap_pos)

    print(p)

    dev.copy(pdf, "06_Heatmap_precursorlevel_log10.pdf", height = 24, width = n_resultsets*8)
    dev.off()

  }

  ###############################
  ## Sample-sample correlation ##
  ###############################
  if (sample_correlation){
    # Corrplot across all runs/analyses
    quant.prot.m.log10.clean.cor = cor(quant.prot.m.log10.clean)
    corrplot(quant.prot.m.log10.clean.cor, method = "color", is.corr = F, type = 'lower', tl.pos = "l", title = "Pairwise Pearson's R")

    pdf("07_Sample_Correlation_CorrPlot_proteinlevel_global.pdf", height = 12, width = 12)
    corrplot(quant.prot.m.log10.clean.cor, method = "color", is.corr = F, type = 'lower', tl.pos = "l", title = "Pairwise Pearson's R")
    dev.off()

    # For each resultset, plot the correlation for the first 4 samples
    resultsets = unique(dt.raw.a$resultset)
    i = 1
    for (i in seq_along(resultsets)){
      dt.i = quant.prot.m.log10.clean[, grep(resultsets[i], colnames(quant.prot.m.log10.clean)) ]

      # Single-resultset corrplot
      dt.i.c = cor(dt.i)
      pdf(paste0("07_Sample_Correlation_CorrPlot_proteinlevel_",resultsets[i],".pdf"), height = 12, width = 12)
      corrplot.mixed(dt.i.c, lower = "color", upper = "number", mar = c(3,3,3,3), tl.pos = "l")
      dev.off()

      # Single-resultset corr pairplot, S01-04
      pdf(paste0("07_Sample_Correlation_CorrPairPlot_proteinlevel_s1to4_",resultsets[i],".pdf"), height = 10, width = 10)
      pairs.panels(dt.i[rowSums(dt.i[,1:4])>0,1:4], smoother = T, main = resultsets[i])
      dev.off()

      # Single-resultset corr pairplot, S01-02
      pdf(paste0("07_Sample_Correlation_CorrPairPlot_proteinlevel_1to2_AllVsCompleteObs",resultsets[i],".pdf"), height = 6, width = 6)
      pairs.panels(dt.i[dt.i[,1]>0 & dt.i[,2]>0,1:2], smoother = T, main = paste(resultsets[i], " -- Shared Observations"))
      pairs.panels(dt.i[rowSums(dt.i[,1:2])>0 ,1:2], smoother = T, main = paste(resultsets[i], " -- All Observations"))
      dev.off()
    }

    if(compare_precursor_level){
      # Corrplot across all runs/analyses
      quant.prec.m.log10.clean.cor = cor(quant.prec.m.log10.clean)
      corrplot(quant.prec.m.log10.clean.cor, method = "color", is.corr = F, type = 'lower', tl.pos = "l", title = "Pairwise Pearson's R")

      pdf("07_Sample_Correlation_CorrPlot_precursorlevel_global.pdf", height = 12, width = 12)
      corrplot(quant.prec.m.log10.clean.cor, method = "color", is.corr = F, type = 'lower', tl.pos = "l", title = "Pairwise Pearson's R")
      dev.off()

      # For each resultset, plot the correlation for the first 4 samples
      resultsets = unique(dt.raw.a$resultset)
      i = 1
      for (i in seq_along(resultsets)){
        dt.i = quant.prec.m.log10.clean[, grep(resultsets[i], colnames(quant.prec.m.log10.clean)) ]

        # Single-resultset corrplot
        dt.i.c = cor(dt.i)
        pdf(paste0("07_Sample_Correlation_CorrPlot_precursorlevel_",resultsets[i],".pdf"), height = 12, width = 12)
        corrplot.mixed(dt.i.c, lower = "color", upper = "number", mar = c(3,3,3,3), tl.pos = "l")
        dev.off()

        # Single-resultset corr pairplot, S01-04
        pdf(paste0("07_Sample_Correlation_CorrPairPlot_precursorlevel_s1to4_",resultsets[i],".pdf"), height = 10, width = 10)
        pairs.panels(dt.i[rowSums(dt.i[,1:4])>0,1:4], smoother = T, main = resultsets[i])
        dev.off()

        # Single-resultset corr pairplot, S01-02
        pdf(paste0("07_Sample_Correlation_CorrPairPlot_precursorlevel_1to2_AllVsCompleteObs",resultsets[i],".pdf"), height = 6, width = 6)
        pairs.panels(dt.i[dt.i[,1]>0 & dt.i[,2]>0,1:2], smoother = T, main = paste(resultsets[i], " -- Shared Observations"))
        pairs.panels(dt.i[rowSums(dt.i[,1:2])>0 ,1:2], smoother = T, main = paste(resultsets[i], " -- All Observations"))
        dev.off()
      }
    }
  }
  return(quant)
  setwd(prev_wd) #Go back to original wd, if it was different
}

