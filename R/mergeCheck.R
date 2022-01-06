# Do not use scientific notation
options(scipen=999)
library(argparse)
parser <- ArgumentParser()
parser$add_argument("--testing",
                    action = "store_true",
                    help = "Use testing files")
parser$add_argument("--syn_user",
                    help = "Synapse username")
parser$add_argument("--syn_pass",
                    help = "Synapse password")
args <- parser$parse_args()
genie_user <- args$syn_user
genie_pass <- args$syn_pass
testing <- args$testing

library(synapser)
library(VariantAnnotation)

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
# Run quick unit test to make sure that the function hasn't be changed
source(file.path(working_dir, "test_flag_variants.R"))

# login to synapse
tryCatch({
  synLogin()
}, error = function(err) {
  #genieUser = Sys.getenv("GENIE_USER")
  #geniePass = Sys.getenv("GENIE_PASS")
  synLogin(genie_user, genie_pass)
})

# limits
variant_limit = 100000
tbl_size_limit = 500

#testing = as.logical(args[1])
if (testing) {
  databaseSynIdMappingId = 'syn11600968'
} else {
  databaseSynIdMappingId = 'syn10967259'
}
databaseSynIdMapping = synTableQuery(sprintf('select * from %s', databaseSynIdMappingId),
                                     includeRowIdAndRowVersion = F)
databaseSynIdMappingDf = synapser::as.data.frame(databaseSynIdMapping)
sampleSynId = databaseSynIdMappingDf$Id[databaseSynIdMappingDf$Database == "sample"]
mutationsInCisSynId =
  databaseSynIdMappingDf$Id[databaseSynIdMappingDf$Database == "mutationsInCis"]
mafSynId = databaseSynIdMappingDf$Id[databaseSynIdMappingDf$Database == "vcf2maf"]

centersTable = synTableQuery(sprintf('select distinct CENTER from %s', sampleSynId),
                             includeRowIdAndRowVersion = F)
centers = synapser::as.data.frame(centersTable)
centers = centers$CENTER

centerMappingSynId = databaseSynIdMappingDf$Id[databaseSynIdMappingDf$Database == "centerMapping"]
centerMapping = synTableQuery(sprintf('select * from %s where release is true',
                                      centerMappingSynId),
                              includeRowIdAndRowVersion = F)
centerMappingDf = synapser::as.data.frame(centerMapping)
centerMappingDf$mutationInCisFilter = as.character(centerMappingDf$mutationInCisFilter)
# Have this here for now until the annotations are changed in the center mapping folder
centerMappingDf$mutationInCisFilter[centerMappingDf$mutationInCisFilter == "TRUE"] = "ON"
centerMappingDf$mutationInCisFilter[centerMappingDf$mutationInCisFilter == "FALSE"] = "OFF"

# ON - Filter is on: TOSS SAMPLES
# OFF - Filter is off: KEEP SAMPLES
# FLAG - Filter is on, but don't toss samples: FLAG BUT KEEP SAMPLES

for (center in centers) {
  print(center)
  # read aggregated clinical data
  clin_query_string = sprintf("select SAMPLE_ID from %s where CENTER = '%s'",
                              sampleSynId, center)
  genieClinTable = synTableQuery(clin_query_string, includeRowIdAndRowVersion = F)
  genieClinData = synapser::as.data.frame(genieClinTable)
  maf_query_string = sprintf("select Tumor_Sample_Barcode, count(Tumor_Sample_Barcode) from %s where Center = '%s' group by Tumor_Sample_Barcode",
                             mafSynId, center)
  mafSampleCountTable = synTableQuery(maf_query_string, includeRowIdAndRowVersion = F)
  mafSampleCount = synapser::as.data.frame(mafSampleCountTable)
  mafSampleCount =
    mafSampleCount[mafSampleCount$Tumor_Sample_Barcode %in% genieClinData$SAMPLE_ID,]
  total = 0
  splitBySamples = c()
  samplesToQuery = c()
  for (sample in mafSampleCount$Tumor_Sample_Barcode) {
    total = total + mafSampleCount[mafSampleCount$Tumor_Sample_Barcode == sample, "COUNT(Tumor_Sample_Barcode)"]
    samplesToQuery = c(sample,samplesToQuery)
    if (total > variant_limit || sample == mafSampleCount$Tumor_Sample_Barcode[nrow(mafSampleCount)]) {
      splitBySamples = c(splitBySamples, paste(samplesToQuery,collapse = "','"))
      total = 0
      samplesToQuery = c()
    }
  }
  for (querySamples in splitBySamples) {
    samplesToRun = unlist(strsplit(querySamples, "','"))
    # read aggregated MAF file
    genieMutTable = synTableQuery(sprintf("SELECT Center,Tumor_Sample_Barcode,Hugo_Symbol,HGVSp_Short,Variant_Classification,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2,t_depth,t_alt_count,End_Position,Protein_position FROM %s where Tumor_Sample_Barcode in ('%s') and inBED is true and Annotation_Status = 'SUCCESS'",
                                          mafSynId, querySamples),
                                  includeRowIdAndRowVersion = F)
    #genieMutTable = synTableQuery(sprintf("SELECT Center,Tumor_Sample_Barcode,Hugo_Symbol,HGVSp_Short,Variant_Classification,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2,t_depth,t_alt_count,End_Position,Protein_position FROM %s where Tumor_Sample_Barcode in ('%s')", mafSynId, querySamples),includeRowIdAndRowVersion=F)
    
    genieMutData = synapser::as.data.frame(genieMutTable)
    flag_variants_to_merge(genieMutData, genieClinData, samplesToRun, upload = TRUE)
    #write.csv(rbind(annotated_df[is.na(annotated_df$Flag),],new_rows), "Missing_variant_annotation.csv", row.names=F)
  }#FOR LOOP END FOR QUERY SAMPLES
}
