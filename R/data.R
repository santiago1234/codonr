#' Codon Optimality Code
#'
#' A data set containing the codon optimality code as defined in
#' \href{http://emboj.embopress.org/content/35/19/2087.long}{Bazzini2016}
#'
#' @format A data frame with 64 rows and 2 variables:
#' \describe{
#'   \item{codon}{3 nucleotide letter codon}
#'   \item{optimality}{optimality definition (stop, neutral, optimal, non-optimal)}
#' }
#'
"optimality_code_embo2016"

#' @title Train data
#' @description Data used to train mRNA stability predictor, this data set
#' was splitted such that the same gene is not in the training and testing
#' sets at the same time
#' @format A data frame with 67775 rows and 7 variables:
#' \describe{
#'   \item{\code{gene_id}}{character ensembl gene id}
#'   \item{\code{specie}}{character the specie}
#'   \item{\code{cell_type}}{character cell type where mRNA stability comes}
#'   \item{\code{datatype}}{character technique used to generate mRNA stability profile}
#'   \item{\code{decay_rate}}{double mRNA decay rate, scaled for each stability profile}
#'   \item{\code{utrlenlog}}{double log length of 3' UTR}
#'   \item{\code{coding}}{character ORF star to stop codon}
#'}
#' @source \url{https://github.com/santiago1234/MZT-rna-stability/blob/master/results/19-04-30-PredictiveModelDecayAllSpecies/19-04-30-EDA/EDAanalysos.md}
"train_set"

#' @title Test data
#' @description Data used to evaluate mRNA stability predictor, this data set
#' was splitted such that the same gene is not in the training and testing
#' sets at the same time
#' @format A data frame with 67775 rows and 7 variables:
#' \describe{
#'   \item{\code{gene_id}}{character ensembl gene id}
#'   \item{\code{specie}}{character the specie}
#'   \item{\code{cell_type}}{character cell type where mRNA stability comes}
#'   \item{\code{datatype}}{character technique used to generate mRNA stability profile}
#'   \item{\code{decay_rate}}{double mRNA decay rate, scaled for each stability profile}
#'   \item{\code{utrlenlog}}{double log length of 3' UTR}
#'   \item{\code{coding}}{character ORF star to stop codon}
#'}
#' @source \url{https://github.com/santiago1234/MZT-rna-stability/blob/master/results/19-04-30-PredictiveModelDecayAllSpecies/19-04-30-EDA/EDAanalysos.md}
"test_set"
