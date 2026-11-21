# -----------------------------------------------------------------------------
# File:    tests/testthat/test-real_data.R
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Unit tests of R/real_data.R: the reader on a made-up domain file
#          and on the six domains, the frequencies of a small vector, and
#          the hotspot counts of short chains on made-up domains.
# Inputs:  data/domains/<domain>.txt
# Outputs: None (made-up files are written to temporary folders).
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

test_that("read_domain reads N, the counts and the counts with position 0", {
  folder <- new_folder()
  writeLines(c("3", "0:0:", "1:2:GENE1|p.A1B,GENE2|p.C2D", "2:0:", "3:5:GENE3|p.E3F"),
             file.path(folder, "made.txt"))
  domain <- read_domain("made", folder)
  expect_equal(domain$N, 3)
  expect_equal(domain$counts, c(2, 0, 5))
  expect_equal(domain$with_position_0, c(0, 2, 0, 5))
})

test_that("read_domain reads the six domains of the paper", {
  expect_equal(DOMAINS, c("cd00054", "cd00083", "cd00590", "smart00179", "smart00209", "cd00063"))
  for (name in DOMAINS) {
    domain <- read_domain(name, PATHS$domains)
    expect_length(domain$counts, domain$N)
    expect_equal(domain$with_position_0, c(0, domain$counts), info = name)
    expect_true(all(domain$counts >= 0 & domain$counts == round(domain$counts)), info = name)
  }
  expect_equal(read_domain("cd00083", PATHS$domains)$N, 93)
})

test_that("count_frequencies counts the positions at 0, 1, ..., max", {
  expect_equal(count_frequencies(c(0, 0, 2, 3, 3)), c("0" = 2, "1" = 0, "2" = 1, "3" = 2))
  expect_equal(sum(count_frequencies(DOMAIN_COUNTS$made1)), length(DOMAIN_COUNTS$made1))
})

test_that("domain_bayes_rejections sets its own seed and finds the large counts", {
  x <- DOMAIN_COUNTS$made1
  set.seed(1)
  first <- domain_bayes_rejections(x, "NP", SMALL_SETTINGS)
  set.seed(2)
  runif(5)
  second <- domain_bayes_rejections(x, "NP", SMALL_SETTINGS)
  expect_equal(first, second)   # the state of the random numbers before the call does not matter
  expect_equal(first, sum(x > 10))
  # the same count from the sampler, seeded here
  set.seed(SMALL_SETTINGS$real_seed)
  fit <- bayes_local_fdr(x, "NP", 1.5, sampler_settings(lambda_init = 1, theta_init = 0.5),
                         draws = 40, burn = 10, version = "2018-03-25")
  expect_equal(first, sum(fit$frequency[fit$fdr <= 0.05]))
  expect_true(domain_bayes_rejections(x, "GP", SMALL_SETTINGS) %in% 0:length(x))
})

test_that("domain_rejections has one row per domain and one column per method", {
  folder <- new_folder()
  for (name in names(DOMAIN_COUNTS)) write_domain_file(folder, name, DOMAIN_COUNTS[[name]])
  paths <- list(domains = folder)
  table <- domain_rejections(names(DOMAIN_COUNTS), c("NP", "GP"), SMALL_SETTINGS, paths, cores = 1)
  expect_named(table, c("domain", "N", "EB", "NP", "GP"))
  expect_equal(table$domain, c("made1", "made2"))
  expect_equal(table$N, as.numeric(lengths(DOMAIN_COUNTS)))
  for (i in 1:2) {
    x <- DOMAIN_COUNTS[[i]]
    expect_equal(table$EB[i], eb_rejections(x, 0.05))
    expect_equal(table$NP[i], domain_bayes_rejections(x, "NP", SMALL_SETTINGS))
    expect_equal(table$GP[i], domain_bayes_rejections(x, "GP", SMALL_SETTINGS))
  }
  # the same table from two parallel processes
  expect_equal(domain_rejections(names(DOMAIN_COUNTS), c("NP", "GP"), SMALL_SETTINGS, paths, cores = 2),
               table)
})

test_that("domain_rejections stops when a chain fails", {
  folder <- new_folder()
  write_domain_file(folder, "made1", DOMAIN_COUNTS$made1)
  expect_error(suppressWarnings(domain_rejections("made1", "XX", SMALL_SETTINGS, list(domains = folder))))
})
