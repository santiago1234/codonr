sequence <- "ACGCACGCACGCACGC"
k <- 4
counts <- count_kmers(sequence, k = k)

test_that("correct counts", {
  expect_equal(sum(counts$n), nchar(sequence) - k + 1)
  acgc <- dplyr::filter(counts, kmer == "ACGC")$n
  expect_equal(acgc, 4)
})

test_that("correct length", {
  expect_true(all(nchar(counts$kmer) == k))
})
