library(argparse)
parser <- ArgumentParser()
parser$add_argument("release",
                    help = "Release version (ie. 5.3-consortium)")
parser$add_argument("--database_synid_mappingid",
                    default = "syn10967259",
                    help = "Database configuration mapping Synapse id")
parser$add_argument("--syn_user",
                    help = "Synapse username")
parser$add_argument("--syn_pass",
                    help = "Synapse password")
args <- parser$parse_args()

library(synapser)
library(VariantAnnotation)
library(knitr)
library(glue)
library(data.table)
library(dplyr)

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
source(file.path(working_dir, "test_uniq_mutations.R"))


# Start of module
release <- args$release
genie_user <- args$syn_user
genie_pass <- args$syn_pass
database_synid_mappingid <- args$database_synid_mappingid

# Login to Synapse
tryCatch({
  synLogin()
}, error = function(err) {
  synLogin(genie_user, genie_pass)
})

# Get release folder synapse id
release_folder_synid <- get_release_folder_synid(database_synid_mappingid, release)
# Get release files mapping
release_files_mapping = get_file_mapping(release_folder_synid)

# Does this number have to be 5 panels from different sites?
NUMBER_CENTERS_COVER_REGION = 6
NUMBER_CENTERS_WITH_CODE = 10
UNIQUE_MUTATIONS_FOLDER_SYNID = "syn18455620"

release_uniq_mutations_ent = synStore(Folder(release,
                                             parentId = UNIQUE_MUTATIONS_FOLDER_SYNID))
center_with_code_folder_ent = synStore(Folder(sprintf("%s_centers_with_code",
                                                      NUMBER_CENTERS_WITH_CODE),
                                              parentId = release_uniq_mutations_ent$properties$id))
uniq_mutation_folder_ent = synStore(Folder(sprintf("regions_covered_by_%s_centers",
                                                   NUMBER_CENTERS_COVER_REGION),
                                           parentId = center_with_code_folder_ent$properties$id))


bed_synid = release_files_mapping[['genomic_information.txt']]
maf_synid = release_files_mapping[['data_mutations_extended.txt']]
sample_synid = release_files_mapping[['data_clinical_sample.txt']]

bed_ent = synGet(bed_synid, followLink=T)
maf_ent = synGet(maf_synid, followLink=T)
sample_ent = synGet(sample_synid, followLink=T)

clinicaldf = fread(sample_ent$path, skip=4)
clinicaldf$CENTER = sapply(strsplit(clinicaldf$SEQ_ASSAY_ID, "-"), function(x) x[1])
clinicaldf = clinicaldf[,c("SAMPLE_ID", "CENTER", "SEQ_ASSAY_ID", "ONCOTREE_CODE")]

mafdf = fread(maf_ent$path)
mafdf = mafdf[,c("Hugo_Symbol", "Chromosome", "Start_Position",
                 "End_Position", "Center", "Tumor_Sample_Barcode",
                 "HGVSp_Short")]

beddf = fread(bed_ent$path)
beddf$CENTER = sapply(strsplit(beddf$SEQ_ASSAY_ID, "-"),
                      function(x) x[1])

# Get number of oncotree codes per center
# It doesn't make sense to check if a specific code is seen across a n number of panels
# because a center could have 5 panels... So the thresholding should be done on
# a specific code seen across n number of centers
codes_per_center = table(clinicaldf$ONCOTREE_CODE, clinicaldf$CENTER)
# Gets the number of centers with a certain codes
codes_count_across_centers = apply(codes_per_center, 1, function(x) {
  sum(x > 0)
})

# Number of panels that has a specific oncotree code
# number_centers_with_code =  floor(0.75*length(unique(clinicaldf$CENTER)))

codes_above_threshold_center = codes_count_across_centers[
  codes_count_across_centers >= NUMBER_CENTERS_WITH_CODE]

unique_mutation_folder = tempdir()

unique_mutation_files = find_unique_mutations(clinicaldf,
                                              codes_above_threshold_center,
                                              mafdf, beddf,
                                              threshold=NUMBER_CENTERS_COVER_REGION,
                                              dir = unique_mutation_folder)

all_unique_mutationsdf = data.frame()
for (mut_file in unique_mutation_files) {
  mutdf = fread(mut_file)
  all_unique_mutationsdf = rbind(all_unique_mutationsdf, mutdf)
}

# Write out each centers unique mutations into center staging folder
processing = all_unique_mutationsdf %>%
  group_by(CENTER) %>%
  group_walk(~ write_and_store_mutations(.x, database_synid_mappingid,
                                         .y$CENTER))
