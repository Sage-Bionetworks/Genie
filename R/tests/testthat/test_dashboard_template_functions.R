# tests for dashboard_template_functions.R

source("../../dashboard_template_functions.R")

library(synapser)
library(testthat)

sample_counts_table <- function() {
  data <- data.frame(
    Center = factor(c("GOLD","SAGE", "TEST"), 
                      levels = c("GOLD", "SAGE", "TEST")),
    Counts = c(1, 2, 1)
  )
  return(data)
}

empty_counts_table <- function() {
  data <- data.frame(
    Center = logical(),
    Counts = logical()
  )
  return(data)
}

sample_maf_table <- function() {
  data <- data.frame(
    Center = c("TEST", "TEST", "SAGE", "SAGE", "GOLD", "BRONZE"),
    Tumor_Sample_Barcode = c("SAGE1", "SAGE2", "SAGE3", "SAGE4", "SAGE5", "SAGE6"),
    Annotation_Status = c("SUCCESS", "FAILED", "FAILED", "FAILED", "FAILED", "SUCCESS")
  )
  return(data)
}

sample_maf_table_no_failed_annotations <- function() {
  data <- data.frame(
    Center = c("TEST", "SAGE", "GOLD"),
    Tumor_Sample_Barcode = c("SAGE1", "SAGE2", "SAGE3"),
    Annotation_Status = c("SUCCESS", "SUCCESS", "SUCCESS")
  )
  return(data)
}


test_that("get_syn_id_from_mapped_database_gets_correct_value", {
  synLogin()
  result <- get_syn_id_from_mapped_database(
    database_name = "main",
    database_synid_mappingid = "syn11600968"
  )
  expect_equal(result, "syn7208886")
})


test_that("get_failed_annotation_table_counts_returns_expected_output", {
  result <- get_failed_annotation_table_counts(
    maf_data=sample_maf_table(), 
    group_by_cols="Center", 
    counts_col_name="Counts")
  expect_equal(result, sample_counts_table())
})

test_that("get_failed_annotation_table_counts_returns_empty_table_with_no_failed_annotations", {
  result <- get_failed_annotation_table_counts(
    maf_data=sample_maf_table_no_failed_annotations(), 
    group_by_cols="Center", 
    counts_col_name="Counts")
  expect_equal(result, empty_counts_table())
})

