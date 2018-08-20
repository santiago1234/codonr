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


