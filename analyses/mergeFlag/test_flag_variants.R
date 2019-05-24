get_working_dir <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
  } else {
    # 'source'd via R console
    return("./")
  }
}
working_dir = get_working_dir()
source(file.path(working_dir, "mergecheck_functions.R"))
library(testthat)
library(VariantAnnotation)
genieMutData = matrix(nrow = 2, ncol = 13)
colnames(genieMutData) = c("Chromosome", "Hugo_Symbol", "Start_Position", "End_Position", "Reference_Allele",
                           "Tumor_Seq_Allele2", "t_depth", 't_alt_count', "Tumor_Sample_Barcode", 
                           "Protein_position", "HGVSp_Short", "Variant_Classification", "Center")


genieMutData[,'Chromosome'] = c("1", "1")
genieMutData[,'Start_Position'] = c("1", "3")
genieMutData[,'End_Position'] = c("1", "3")
genieMutData[,'Reference_Allele'] = c("A", "C")
genieMutData[,'Tumor_Seq_Allele2'] = c("T", "G")
genieMutData[,'t_alt_count'] = c(1, 2)
genieMutData[,'t_depth'] = c(100, 100)
genieMutData[,'Tumor_Sample_Barcode'] = c("SAGE1", "SAGE1")
genieMutData[,'Protein_position'] =  c("3/4", "4/4")
genieMutData[,'Hugo_Symbol'] = c("ROS1", "ROS1")
genieMutData[,'Center'] = c("SAGE", "SAGE")

genieMutData = as.data.frame(genieMutData, stringsAsFactors = F)
genieClinData = data.frame(SAMPLE_ID = c("SAGE1"),
                           stringsAsFactors = F)

test_that("Mutations are flagged, same starts and ends", {
  tbl = flag_variants_to_merge(genieMutData, genieClinData, c("SAGE1"), upload=F)
  expected = genieMutData[, c("Center", "Tumor_Sample_Barcode", "Hugo_Symbol",
                              "Variant_Classification", "Chromosome", "Start_Position",
                              "Reference_Allele", "Tumor_Seq_Allele2","t_depth",
                              "HGVSp_Short")]
  expected$t_alt_count_num = c(1, 2)
  # Change everything to factors or else its difficult to validate is absurd
  # This is ok because we don't do anything with factors after the function is ran
  expected <- data.frame( lapply( expected , factor ))
  tbl <- data.frame( lapply( tbl , factor ))
  expect_equal(tbl, expected[,colnames(tbl)])
})

genieMutData$Start_Position[2] = "15"
genieMutData$End_Position[2] = "15"
test_that("Mutations are not flagged", {
  tbl = flag_variants_to_merge(genieMutData, genieClinData, c("SAGE1"), upload=F)
  expect_equal(nrow(tbl), 0)
  expect_equal(colnames(tbl), c("Center", "Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSp_Short",
                                "Variant_Classification", "Chromosome", "Start_Position",
                                "Reference_Allele", "Tumor_Seq_Allele2", "t_alt_count_num", "t_depth"))
  
})

genieMutData$Start_Position = c("1", "10")
genieMutData$End_Position = c("5", "12")
test_that("Mutations are flagged, different starts and ends", {
  tbl = flag_variants_to_merge(genieMutData, genieClinData, c("SAGE1"), upload=F)
  expected = genieMutData[, c("Center", "Tumor_Sample_Barcode", "Hugo_Symbol",
                              "Variant_Classification", "Chromosome", "Start_Position",
                              "Reference_Allele", "Tumor_Seq_Allele2","t_depth",
                              "HGVSp_Short")]
  expected$t_alt_count_num = c(1, 2)
  # Change everything to factors or else its difficult to validate is absurd
  # This is ok because we don't do anything with factors after the function is ran
  expected <- data.frame( lapply( expected , factor ))
  tbl <- data.frame( lapply( tbl , factor ))
  expect_equal(tbl, expected[,colnames(tbl)])
})


genieMutData$Start_Position = c("1", "10")
genieMutData$End_Position = c("5", "12")
test_that("Mutations not flagged, different starts and ends", {
  tbl = flag_variants_to_merge(genieMutData, genieClinData, c("SAGE1"), upload=F)
  expected = genieMutData[, c("Center", "Tumor_Sample_Barcode", "Hugo_Symbol",
                              "Variant_Classification", "Chromosome", "Start_Position",
                              "Reference_Allele", "Tumor_Seq_Allele2","t_depth",
                              "HGVSp_Short")]
  expected$t_alt_count_num = c(1, 2)
  # Change everything to factors or else its difficult to validate is absurd
  # This is ok because we don't do anything with factors after the function is ran
  expected <- data.frame( lapply( expected , factor ))
  tbl <- data.frame( lapply( tbl , factor ))
  expect_equal(tbl, expected[,colnames(tbl)])
})

