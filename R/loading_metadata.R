#Read in annotations .csv file & match file names with annotations
#' Read in annotations & match sample names in counts dataframe with metadata
#'
#' This function works to read in the a metadata .csv file as a dataframe &
#' -then matches a certain sample name with its associated label & replaces
#' -a 'messy' label with the clean label (e.g. sample1_trimmed_sorted.txt to sample1_treated).
#' Finally, the data are sorted to ensure all data is in the same order in both count data frame & metadata.
#' NOTE: This function uses stringr::str_detect so match_col need to be unique for each sample to act as a unique
#' identifier for the function, otherwise it will fail & match incorrectly.
#'
#' @param count_data Data frame. The dataframe read in from featureCounts- or HTSeq-generated reads.
#' @param metadata_path Character. The path to the metadata .csv file.
#' @param match_col Character. The matching column name by which the counts data & metadata will be matched (essentially a barcode to match samples). MUST BE PRESENT IN METADATA FILE.
#' @param clean_col Character. The column by which samples in counts data will be renames by (the clean labelling of samples). MUST BE PRESENT IN METADATA FILE.
#'
#' @returns A a named list containing both the cleaned up & matched counts data & metadata.
#' @export
#'
match_annotations <- function(count_data, metadata_path, match_col, clean_col) {

  #Read in annotation .csv file as a dataframe
  metadata <- as.data.frame(utils::read.csv(metadata_path, header = TRUE, sep = ","))

  #Store sample names needed for cleaning
  current_cols <- colnames(count_data)
  new_cols <- current_cols

  #Loop through every row in the metadata & store replacement sample names
  for (i in seq_len(nrow(metadata))) {
    pattern <- metadata[[match_col]][i]
    replacement <- metadata[[clean_col]][i]

    matches <- stringr::str_detect(current_cols, pattern)
    new_cols[matches] <- replacement
  }

  colnames(count_data) <- new_cols

  #Generate warning if any of the rows are left unchanged after the loop
  unchanged <- current_cols == new_cols
  if (any(unchanged)) {
    warning(paste(
      "WARNING: The following columns could not be mapped to the metadata & were left unchanged:",
      paste(current_cols[unchanged], collapse = ",")
    ))
  }

  #Reorder metadata rows to ensure each row corresponds to the appropriate sample
  idx <- match(colnames(count_data), metadata[[clean_col]])

  #Stop running if any columns in the counts dataframe cannot be found in the metadata
  if (any(is.na(idx))) {
    stop("ERROR: Some columns within the counts data could not be found in the metadata. Please check if columns are matched.")
  }

  metadata_sorted <- metadata[idx, ]

  if (!all(colnames(count_data) == metadata_sorted[[clean_col]])) {
    stop("ERROR: Final alignments do not match. Columns & rows of count data & metadata do not match.")
  }

  #Ensure row names of metadata match the clean column
  rownames(metadata_sorted) <- metadata_sorted[[clean_col]]

  if (!all(rownames(metadata_sorted) == metadata_sorted[[clean_col]])) {
    stop("ERROR: Something went wrong. Metadata row names do not match the specified column.")
  }

  message("Yay! Count data & metadata are now cleaned & aligned!")

  return(list(
    counts_final = count_data,
    metadata_final = metadata_sorted
  ))
}

#Subset samples for condition-wise examination
#' Subsetting samples of interest for condition-wise analysis
#'
#' This function allows for the user to subset their counts data by condition/sample name-
#' -by providing the function a character vector of the samples names of interest.
#' It is able to match the data with the metadata & order it so that all data are in order.
#'
#' @param sample_list Character vector. A list of specific samples by which analysis will be undertaken (e.g. all control & treated in one cell line).
#' @param count_data Data frame. The cleaned count data frame from the match_annotations function.
#' @param metadata Data frame. The aligned metadata from the match_annotations function.
#' @param matching_sample_names Character. The column within the metadata data frame by which samples will be subsetted (i.e. cleaned sample names/clean_col).
#'
#' @returns A named list containing the subsetted & matched counts data frame & metadata.
#' @export
#'
subset_of_interest <- function(sample_list, count_data, metadata, matching_sample_names) {

  #Subset data to contain only samples of interest
  samples_of_interest <- count_data[, sample_list, drop = FALSE]

  #Subset & order metadata annotations to exactly match count data
  samples_of_interest_labels <- metadata[metadata[[matching_sample_names]] %in% sample_list, ]

  #Order metadata annotations to match count data
  idx <- match(colnames(samples_of_interest), samples_of_interest_labels[[matching_sample_names]])

  if (any(is.na(idx))) {
    stop("ERROR: Some samples in samples of interest were not found in the metadata.")
  }

  samples_of_interest_labels <- samples_of_interest_labels[idx, ]

  #Undergo sanity check to ensure if rownames of metadata perfectly matches the column/sample names of the count data
  if (!all(samples_of_interest_labels[[matching_sample_names]] == colnames(samples_of_interest))) {
    stop("ERROR: Something went wrong, metadata rows do not match count data columns.")
  }

  message("Yay! Subsetted samples having matched & ordered count data & metadata annotations")

  return(list(
    matched_data_of_interest = samples_of_interest,
    matched_labels_of_interest = samples_of_interest_labels
  ))
}
