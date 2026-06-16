#Load in HTSeq-generated counts
#' Read in & combine all output files from HTSeq
#'
#' This function works to combine all HTSeq output files from a folder of a specified directory.
#' It combines all files into a single dataframe, setting the row names as the Geneids & removing summary rows at the end.
#' Allows for data to be ready for downstream analysis
#'
#' @param dir_path Character. Path to folder containing all .txt files from HTSeq output.
#' @param file_pattern Character. The regex pattern to match files (default is "\\.txt$")
#'
#' @returns A dataframe of raw counts for each gene ID (as rows) & sample (as columns).
#' @export
#'
read_htseq_data <- function(dir_path, file_pattern = "\\.txt$") {

  htseq_files <- list.files(
    path = dir_path,
    pattern = file_pattern,
    full.names = TRUE
  )

  if (length(htseq_files) == 0) {
    stop("ERROR: No files detected in defined file path.")
  }

  #Extract file names, cleaning sample/file names will be dependent case-by-case
  sample_names <- basename(htseq_files)
  sample_names <- tools::file_path_sans_ext(sample_names)

  #Set names of HTSeq files as sample names so they can be
  #-matched with the metadata
  names(htseq_files) <- sample_names

  #Generate starting dataframe to enter raw values
  htseq_df <- tibble::tibble(
    Sample = sample_names,
    file_path = htseq_files
  )

  #Combine & format all of the HTSeq files by Gene ID & Count
  htseq_data <- htseq_df|>
    dplyr::mutate(data = purrr::map(file_path, ~ readr::read_tsv(.x,
                                                                 col_names = c("Geneid", "Count"),
                                                                 col_types = "cd",
                                                                 show_col_types = FALSE))) |>
    tidyr::unnest(cols = c(data)) |>
    dplyr::select(-file_path) |>
    tidyr::pivot_wider(names_from = Sample, values_from = Count, values_fill = list(Count = 0)) |>
    dplyr::filter(!stringr::str_detect(Geneid, "^__")) |>
    dplyr::filter(!is.na(Geneid)) |>
    as.data.frame()

  #Set row names as Gene ID
  rownames(htseq_data) <- htseq_data$Geneid
  htseq_data$Geneid <- NULL

  return(htseq_data)
}

#Load in featureCounts-generated counts
#' Read in & format featureCounts generated counts
#'
#' This function imports the featureCounts output .txt file & formats it by-
#' -removing any N/A counts, irrelevant metadata columns, & formatting so that-
#' -rows correspond to gene IDs & columns correspond to samples.
#'
#' @param dir_path Character. Path to .txt file containing featureCounts-generated reads.
#'
#' @returns A dataframe containing all raw counts, where rows correspond to gene IDs, & columns correspond to samples.
#' @export
#'
read_featurecounts_data <- function(dir_path) {

  featurecounts_data <- utils::read.delim(dir_path, skip = 1, header = TRUE, stringsAsFactors = FALSE)

  #Filter out any genes that are N/A & convert the data into a dataframe
  featurecounts_data <- featurecounts_data |>
    dplyr::filter(!is.na(Geneid)) |>
    as.data.frame()

  #Set rownames as gene IDs
  rownames(featurecounts_data) <- featurecounts_data$Geneid

  #Remove metadata columns
  featurecounts_data <- featurecounts_data |>
    dplyr::select(-Geneid, -Chr, -Start, -End, -Strand, -Length)

  #Remove any N/A values & replace with 0
  featurecounts_data <- featurecounts_data |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, 0))) |>
    as.data.frame()
  return(featurecounts_data)
}
