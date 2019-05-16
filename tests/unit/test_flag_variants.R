source("flag_variants.R")
library(testthat)
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
genieMutData = as.data.frame(genieMutData, stringsAsFactors = F)
genieClinData = data.frame(SAMPLE_ID = c("SAGE1"),
                           stringsAsFactors = F)

test_that("Mutations are flagged", {
  tbl = flag_variants_to_merge(genieMutData, genieClinData, c("SAGE1"))
  expected = genieMutData[, c("Center", "Tumor_Sample_Barcode", "Hugo_Symbol",
                              "Variant_Classification", "Chromosome", "Start_Position",
                              "Reference_Allele", "Tumor_Seq_Allele2","t_depth",
                              "HGVSp_Short")]
  expected$t_alt_count_num = c(1, 2)
  expect_equal(tbl, expected)
})
tbl$Center
expected$Center

