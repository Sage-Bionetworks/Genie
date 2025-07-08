# tests for mergecheck_functions.R

source("../../mergecheck_functions.R")

library(testthat)


test_that("get_database_to_synapse_mapping_synid_gets_correct_synid_for_testing", {
  result <- get_database_to_synapse_mapping_synid(
    testing = TRUE, staging = FALSE
  )
  expect_equal(result, "syn11600968")
})

test_that("get_database_to_synapse_mapping_synid_gets_correct_synid_for_staging", {
  result <- get_database_to_synapse_mapping_synid(
    testing = FALSE, staging = TRUE
  )
  expect_equal(result, "syn12094210")
})

test_that("get_database_to_synapse_mapping_synid_gets_correct_synid_for_prod", {
  result <- get_database_to_synapse_mapping_synid(
    testing = FALSE, staging = FALSE
  )
  expect_equal(result, "syn10967259")
})

test_that("get_database_to_synapse_mapping_synid_throws_error", {
  expect_error(get_database_to_synapse_mapping_synid(testing = TRUE, staging = TRUE), 
    "Mutation in cis only available in staging or testing mode not both")
})



