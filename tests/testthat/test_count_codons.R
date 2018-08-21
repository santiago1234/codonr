context("count codons")

dna <- "ATGGCAAGGCCCAAAGTGTTTTTCGATCTGACCGCCGGCGGCAGTCCTGTTGGAAGGGTGGTAATGGAG"
output <- count_codons(dna)

test_that('number of counts equals number of codons', {
  expect_equal(sum(output$n), nchar(dna) / 3)
})

test_that('correct number of counts', {
  dna1 <- "ACGacgACGtttACG"
  outpu1 <- count_codons(dna1)
  expect_equal(dplyr::filter(outpu1, codon == "ACG")$n, 4)
  expect_false("CGA" %in% outpu1$codon)
  expect_error(count_codons("A"))
})


context("codon composition")
# test codon composition function -----------------------------------------

## create today data set
set.seed(10)
n <- 100
## generate fake data set for testing
Reduce(function(x, y) paste0(x, y), sample(c("A", "C", "G", "T"), n*3, replace = T))
cds <- replicate(
  n,
  expr = Reduce(
    function(x, y) paste0(x, y),
    sample(c("A", "C", "G", "T"), rpois(n = 1, lambda = 100)*3, replace = T)
    )
  )
ids <- paste0("id_", rnorm(100))
data <- tibble::tibble(ids = ids, cds = cds)
output <- codon_composition(data, id_col = "ids", orf_col = "cds")

test_that('correct output', {
  expect_equal(nrow(data), nrow(output))
  expect_false(any(duplicated(output$id_col)))
  expect_error(codon_composition(dplyr::bind_rows(data, data), , id_col = "ids", orf_col = "cds"))
  # total sum equal length
  result <-
    tidyr::gather(output, key="codon", value="n", -id_col) %>%
    dplyr::group_by(id_col) %>%
    dplyr::summarise(total_codons = sum(n)) %>%
    dplyr::inner_join(data, by=c("id_col" = "ids")) %>%
    dplyr::mutate(expected_codons = nchar(cds) / 3)
  expect_true(all(result$total_codons == result$expected_codons))

})


