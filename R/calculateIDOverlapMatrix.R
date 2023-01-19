#' calculateIDOverlapMatrix
#' @description Calculates % overlap (intersect/union) of ID sets from a proteomics matrix with identifiers in the row.names
#' @param proteomicsMatrix Quantitative proteomics matrix with analyte ID in row.names and measurement id as column names.
#' Either NAs or 0s for missing values
#' @return A distance matrix with ID set overlap in percent
#' @export

calculateIDOverlapMatrix <- function(proteomicsMatrix) {
  # Calculate the number of rows and columns in the matrix
  numRows <- nrow(proteomicsMatrix)
  numCols <- ncol(proteomicsMatrix)

  # Replace NA values with 0s
  proteomicsMatrix[is.na(proteomicsMatrix)] <- 0

  # Initialize the distance matrix with all 0s
  distanceMatrix <- matrix(0, numCols, numCols)
  colnames(distanceMatrix) = colnames(proteomicsMatrix)
  rownames(distanceMatrix) = colnames(proteomicsMatrix)

  # Loop through each pair of columns
  for (i in 1:numCols) {
    for (j in 1:numCols) {
      # Calculate the overlap between the two columns as the number of non-zero values in both columns
      overlap <- sum(proteomicsMatrix[,i] != 0 & proteomicsMatrix[,j] != 0)

      # Calculate the overlap in percent as the overlap divided by the total number of non-zero values in both columns
      overlapPercent <- overlap / (sum(proteomicsMatrix[,i] != 0) + sum(proteomicsMatrix[,j] != 0) - overlap)

      # Store the overlap percent in the distance matrix
      distanceMatrix[i,j] <- overlapPercent
    }
  }

  # Return the distance matrix
  return(distanceMatrix)
}
