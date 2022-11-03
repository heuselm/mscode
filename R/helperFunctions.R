# Auxiliary functions

#' fromListWithNames
#' @description alternative to UpSetR::fromList that maintains element names in row.names of the output detection matrix.
#' Required to track members in intersect sets. Grabbed & modified from https://github.com/hms-dbmi/UpSetR/issues/85.
#' @param input list of named character vectors containing the sets
#' @export
fromListWithNames <- function (input = list("Set1" = c("A","B","C","D","E"),
                                            "Set2" = c("A","B","D","E","F"),
                                            "Set3" = c("B","D","E","F","G","H"))) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

#' getSetAnalysisGroups
#' @description Tracks members in intersect sets. Grabbed & modified from https://github.com/hms-dbmi/UpSetR/issues/85.
#' @param input list of named character vectors containing the sets
getSetAnalysisGroups <- function (listInput = list("Set1" = c("A","B","C","D","E"),
                                            "Set2" = c("A","B","D","E","F"),
                                            "Set3" = c("B","D","E","F","G","H")), sort = TRUE) {
  # listInput could look like this:
  # $one
  # [1] "a" "b" "c" "e" "g" "h" "k" "l" "m"
  # $two
  # [1] "a" "b" "d" "e" "j"
  # $three
  # [1] "a" "e" "f" "g" "h" "i" "j" "l" "m"
  listInputmat    <- fromListWithNames(listInput) == 1
  #     one   two three
  # a  TRUE  TRUE  TRUE
  # b  TRUE  TRUE FALSE
  #...
  # condensing matrix to unique combinations elements
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  # going through all unique combinations and collect elements for each in a list
  for (i in 1:nrow(listInputunique)) {
    currentRow <- listInputunique[i,]
    myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
    myelements
    # attr(,"groups")
    #   one   two three
    # FALSE FALSE  TRUE
    #  f  i
    # 12 13
  }
  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  return(grouplist)
  # save element list to facilitate access using an index in case rownames are not named
}
