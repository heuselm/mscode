ggally_rmse <- function (data, mapping, ..., stars = TRUE, method = "pearson",
                         use = "complete.obs", display_grid = FALSE, digits = 3, title_args = list(...),
                         group_args = list(...), justify_labels = "right", align_percent = 0.5,
                         title = "RMSE", alignPercent = warning("deprecated. Use `align_percent`"),
                         displayGrid = warning("deprecated. Use `display_grid`"))
{
  if (!missing(alignPercent)) {
    warning("`alignPercent` is deprecated. Please use `align_percent` if alignment still needs to be adjusted")
    align_percent <- alignPercent
  }
  if (!missing(displayGrid)) {
    warning("`displayGrid` is deprecated. Please use `display_grid`")
    display_grid <- displayGrid
  }
  na.rm <- if (missing(use)) {
    NA
  }
  else {
    (use %in% c("complete.obs", "pairwise.complete.obs",
                "na.or.complete"))
  }
  ggally_statistic(data = data, mapping = mapping, na.rm = na.rm,
                   align_percent = align_percent, display_grid = display_grid,
                   title_args = title_args, group_args = group_args, justify_labels = justify_labels,
                   justify_text = "left", sep = if ("colour" %in% names(mapping))
                     ": "
                   else ":\n", title = title, text_fn = function(x, y) {
                     # if (is_date(x)) {
                     #   x <- as.numeric(x)
                     # }
                     # if (is_date(y)) {
                     #   y <- as.numeric(y)
                     # }
                     rmse <- sqrt(mean((x-y)^2))
                     rmse_txt <- formatC(rmse, digits = digits, format = "f")
                     rmse_txt
                   })
}
