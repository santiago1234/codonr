#' Get file path of fastq data
#'
#' @param sname character SampleName to get data from
#' @param dir_to_files path to dir with fq
#' @param sample_report tibble: sample report data
#' @importFrom magrittr %>%
#' @return character file path to fastq files
get_files_for_sample <- function(sname, dir_to_files, sample_report) {
  sample_report %>%
    dplyr::filter(SampleName == sname) %>%
    dplyr::pull(Output) %>%
    purrr::map_chr(function(x) file.path(dir_to_files, x))
}


#' Read Files and Saves it with SampleName
#' This functions reads the provided fastq file(s) and saves it to
#' `fout`, by merging the files. It is similar to doing `cat f1.fq f2.fq >files.fq`
#'in the command line
#'
#' @param fqfiles character vector with path to fastq file(s)
#' @param fout name of output file
#' @param ...
#'
#' @return file size
get_files <- function(fqfiles, fout, ...) {
  fqfiles %>%
    purrr::map(
      function(x) {
        fq <- ShortRead::readFastq(x)
        ShortRead::writeFastq(fq, fout, mode="a", ...)
      }
    )

}


#' Get NGS-data (fastq) from LIMS
#' I use this function to move and rename data from LIMS. LIMS data is split into
#' different sub-fastq files for each sample. This function takes the SampleReport.csv data,
#' collects all the subfiles associated with a sample name, merges into a single fastq, and
#' saves the result into a file called as the sample name.
#' @param sreport_path character, the path to sample SampleReport.csv,
#' this file is usually in the same directory where the data is.
#' @param dir_to_files The directory to the flow cell data
#' @param output_dir Where to save the fastq files?
#' @param snames optional character vector to select specific samples
#' @param ...
#'
#' @return NULL
#' @export
#'
#' @examples
#' \donotrun{
#' # this codes will only work in the hickory cluster
#' # provided the data exist
#' sr <- "/n/analysis/Bazzini/arb/MOLNG-2424/HGCC7BGX7/Sample_Report.csv"
#' dir <- "/n/analysis/Bazzini/arb/MOLNG-2424/HGCC7BGX7/"
#' # you should use compress=FALSE, to get uncompress file
#' get_fq_files_from_LIMS(sr, dir, output_dir = "~/", snames = "K562_PCR_A", compress=FALSE)
#' }
get_fq_files_from_LIMS <- function(sreport_path, dir_to_files,
                                   output_dir="./", snames, ...) {
  ## load sample report
  sample_report <- readr::read_csv(sreport_path)
  sample_names <- unique(sample_report$SampleName)

  if (!missing(snames)) {
    sample_names <- sample_names[sample_names %in% snames]
  }
  if(length(sample_names) < 1) stop("snames not found")

  sample_names %>%
    purrr::map(
      function(mysample_name) {
        outname <- file.path(output_dir, stringr::str_c(mysample_name, ".fq"))
        get_files_for_sample(sname = mysample_name, dir_to_files, sample_report) %>%
          get_files(fout = outname, ...)
      }
    )
}
