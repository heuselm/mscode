#' plotTicProfiles
#' @description A function to extract and plot TIC profiles from Thermo.raw mass spectrometry data files

#' @param path_to_raw Path folder containing the .raw files
#' @param study_design Optional, if plots by condition are desired.
#' Study design table as data.table or path to tab-separated .txt with mandatory columns:
#' \itemize{
#' \item filename: raw file name including .raw extension
#' \item condition: String, biological condition (e.g. "treated" and "untreated")
#' \item replicate: Replicate number (integer).
#' }
#' @param write_csv Character, whether .csv output shall be written
#'\itemize{
#' \item next to raw files: "with_rawfiles"
#' \item in the current working directory: "workdir"
#' \item not at all: "none"
#' }
#'
#' @return list with elements tic_table (data.table) and plot (ggplot)

#' @import rawrr data.table ggplot2 plotly htmlwidgets

#' @export

plotTicProfiles <- function(path_to_raw = "../data/dda/",
  study_design = NULL, write_csv = "none") {

  files = list.files(path_to_raw)[grep(".raw$", list.files(path_to_raw))]
  nruns = length(files)

  # read and collect tic profiles using rawrr
  tics = data.table()
  for (i in files){
    message("extracting TIC chromatogram from ",i)
    tic_i = rawrr::readChromatogram(rawfile = paste0(path_to_raw, i), type = "tic")
    tic_i_dt = data.table("Retention_time_min" = as.numeric(tic_i$times),
      "Intensity" = tic_i$intensities, file = i)
    tics = rbind(tics,tic_i_dt)
  }

  if (write_csv == "with_rawfiles"){
    # write table
    fwrite(tics, paste0(path_to_raw, "TICs.csv"))
  } else if (write_csv == "workdir") {
    fwrite(tics,"TICs.csv")
  }

  # Generic TIC plots
  tics[, Retention_time_min:=as.numeric(Retention_time_min)]
  tics[, file:=unlist(strsplit(file, split = "-"))[length(unlist(strsplit(file, split = "-")))], file]

  tics_all = ggplot(tics) + geom_line(aes(x = Retention_time_min, y = Intensity, group = file)) +
    facet_wrap(~file) + theme_bw() + theme(legend.position = "none") + ggtitle("TIC_all")
  plot(tics_all)
  ggsave("TICs_all.pdf", width = nruns/2, height = nruns/2)

  pdf("TICs_separate.pdf", width = 5, height = 3)
  for (i in unique(tics$file)){
    p = ggplot(tics[file == i]) + geom_line(aes(x = Retention_time_min, y = Intensity, group = file)) + theme_bw() +
      theme(legend.position = "none") + ggtitle(paste("TIC", i))
    plot(p)
  }
  dev.off()

  # Interactive plot
  p_int = ggplot(tics) + geom_line(aes(x = Retention_time_min, y = Intensity, group = file, color = file)) +
  theme_bw() + theme(legend.position = "bottom") + ggtitle(paste0("Total Ion Chromatograms\n", getwd(),"\n"))
  htmlwidgets::saveWidget(ggplotly(p_int), file = "TICs_all_interactive.html")

  res = list("tic_table" = tics,
    "plot" = tics_all)

  # Condition-informed plot
  if (!is.null(study_design)){
    tics = merge(tics, study_design, by.x = "file", by.y = "filename", all.y = F)
    tics_all_cond = ggplot(tics) + geom_line(aes(x = Retention_time_min, y = Intensity, group = replicate, col = replicate)) +
    facet_wrap(~condition) + theme_bw() + theme(legend.position = "bottom") + ggtitle("TICs_all")
    plot(tics_all_cond)
    ggsave("TICs_all_by_condition.pdf", width = nruns, height = nruns/2)

    res = list("tic_table" = tics,
    "plot_all" = tics_all,
      "plot_all_cond" = tics_all_cond)
  }

  return(res)
}


