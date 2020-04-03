library(VariantAnnotation)
library(testthat)

get_working_dir <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
  } else {
    # 'source'd via R console
    return("./analyses/genomicData")
  }
}
working_dir = get_working_dir()
source(file.path(working_dir, "uniq_mutation_functions.R"))

beddf = matrix(nrow = 2, ncol = 5)
colnames(beddf) = c("Chromosome", "Start_Position", "End_Position", "CENTER", 'threshold_key')
beddf[,'Chromosome'] = c("1", "1")
beddf[,'Start_Position'] = c(1, 3)
beddf[,'End_Position'] = c(3, 5)
beddf[,'CENTER'] = c("SAGE", "TEST")
beddf[,'threshold_key'] = c("SAGE", "TEST")
beddf = as.data.frame(beddf, stringsAsFactors = F)
beddf$Start_Position <- as.numeric(beddf$Start_Position)
beddf$End_Position <- as.numeric(beddf$End_Position)

mafdf = matrix(nrow = 2, ncol = 3)
colnames(mafdf) = c("Chromosome", "Start_Position", "End_Position")
mafdf[,'Chromosome'] = c("1", "1")
mafdf[,'Start_Position'] = c(2, 4)
mafdf[,'End_Position'] = c(2, 4)
mafdf = as.data.frame(mafdf, stringsAsFactors = F)
mafdf$Start_Position <- as.numeric(mafdf$Start_Position)
mafdf$End_Position <- as.numeric(mafdf$End_Position)


test_that("Check overlapping mutations, no threshold", {
  overlap_mafdf = check_mutation_overlap(beddf, mafdf)
  expect_equal(overlap_mafdf, mafdf)
})

test_that("Check threshold_key works", {
  overlap_mafdf = check_mutation_overlap(beddf, mafdf, threshold_key = "threshold_key")
  expect_equal(overlap_mafdf, mafdf)
})


test_that("Threshold check, at least 2 centers cover a bed region so no overlap returned", {
  overlap_mafdf = check_mutation_overlap(beddf, mafdf, threshold = 2)
  expect_true(nrow(overlap_mafdf) == 0)
})


test_that("Check some overlap mutations", {
  beddf$Start_Position <- c(1, 3)
  beddf$End_Position <- c(1, 5)
  overlap_mafdf = check_mutation_overlap(beddf, mafdf)
  expect_equal(overlap_mafdf, mafdf[2,])
})


test_that("Threshold check, at least 2 centers cover a bed region", {
  beddf$Start_Position <- c(3, 3)
  beddf$End_Position <- c(5, 5)
  overlap_mafdf = check_mutation_overlap(beddf, mafdf, threshold = 2)
  expect_equal(overlap_mafdf, mafdf[2,])
})

clinicaldf = matrix(nrow = 3, ncol = 2)
colnames(clinicaldf) = c("ONCOTREE_CODE", "CENTER")
clinicaldf[,'ONCOTREE_CODE'] = c("TEST", "TEST", "FOO")
clinicaldf[,'CENTER'] = c("SAGE", "TEST", "SAGE")
clinicaldf = as.data.frame(clinicaldf, stringsAsFactors = F)

test_that("Check code is returned when covered by n centers", {
  codes = find_codes_coveredby_n_centers(clinicaldf, threshold = 2)
  expect_equal(codes, c("TEST"))
})

test_that("Check code is not returned when not covered by n centers", {
  codes = find_codes_coveredby_n_centers(clinicaldf, threshold = 3)
  expect_true(length(codes) == 0)
})

test_that("Check code is not returned when not covered by n centers", {
  codes = find_codes_coveredby_n_centers(clinicaldf, threshold = 3)
  expect_true(length(codes) == 0)
})

test_that("Check all codes returned", {
  codes = find_codes_coveredby_n_centers(clinicaldf, threshold = 1)
  expect_true(all( c("FOO", "TEST") %in% codes))
})


