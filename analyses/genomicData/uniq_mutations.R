library(argparse)
parser <- ArgumentParser()
parser$add_argument("release",
                    help = "Release version (ie. 5.3-consortium)")
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
library(knitr)
library(glue)
library(readr)

tryCatch({
  synLogin()
}, error = function(err) {
  #genieUser = Sys.getenv("GENIE_USER")
  #geniePass = Sys.getenv("GENIE_PASS")
  synLogin(genie_user, genie_pass)
})

get_file_mapping = function(release_folder_synid) {
  release_ent = synGet(release_folder_synid)
  print(release_ent$properties$name)
  release_files = synGetChildren(release_folder_synid)
  release_files_list = as.list(release_files)
  file_mapping = sapply(release_files_list, function(release_file) {
    mapping = c()
    if (release_file$name == "genie_combined.bed") {
      release_file$name = "genomic_information.txt"
    }
    mapping[release_file$name] = release_file$id
    mapping
  })
}

check_mutation_overlap <- function(beddf, mutdf, threshold=NA, threshold_key="CENTER") {
  # Must check on this overlap.....
  # Must actually check if the unique mutation is actually unique and seen by other sites
  # parametrize,  2 of the regions cover... or 5 of the regions cover....
  maf_vr = GRanges(seqnames = Rle(paste0("chr", mutdf$Chromosome)),
                   ranges = IRanges(start = mutdf$Start_Position,
                                    end = mutdf$End_Position))
  seqlevels(maf_vr) = sort(seqlevels(maf_vr))
  
  bed_gr = GRanges(seqnames = Rle(paste0("chr",beddf$Chromosome)),
                   ranges = IRanges(start = beddf$Start_Position,
                                    end = beddf$End_Position))
  seqlevels(bed_gr) = sort(seqlevels(bed_gr))
  # Number of other centers that must at least have this overlap
  if (!is.na(threshold)) {
    overlap_with_threshold = sapply(c(1:nrow(mutdf)), function(numrow) {
      overlap = bed_gr %over% maf_vr[numrow]
      
      if (length(unique(beddf[[threshold_key]][overlap])) >= threshold) {
        return(T)
      } else{
        return(F)
      } 
    })
    mutdf[overlap_with_threshold, ]
  } else {
    mutdf[maf_vr %over% bed_gr,]
  }
}

find_unique_mutations_panel <- function(clinicaldf, codes_above_threshold,
                                        mafdf, beddf, dir, threshold=NA) {
  unique_mutation_files = c()
  for (code in names(codes_above_threshold)) {
    print(code)
    uniq_muts_df = data.frame()
    
    # Get oncotree code
    samples = clinicaldf$SAMPLE_ID[clinicaldf$ONCOTREE_CODE == code]
    mutationdf = mafdf[mafdf$Tumor_Sample_Barcode %in% samples, ]
    mutationdf$mutation = paste(mutationdf$Hugo_Symbol, mutationdf$HGVSp_Short)
    merged_mafdf = merge.data.frame(mutationdf, clinicaldf,
                                    by.x = "Tumor_Sample_Barcode",
                                    by.y = "SAMPLE_ID")
    mutation_per_panel = table(merged_mafdf$mutation, merged_mafdf$SEQ_ASSAY_ID)
    
    unique_mutation_check = apply(mutation_per_panel, 1, function(x) {
      sum(x > 0)
    })
    unique_mutation = unique_mutation_check[unique_mutation_check == 1]
    merged_mafdf = merged_mafdf[merged_mafdf$mutation %in% names(unique_mutation), ]
    # Grab bed regions for the panels involved with this mutation
    newbeddf = beddf[beddf$SEQ_ASSAY_ID %in% unique(merged_mafdf$SEQ_ASSAY_ID), ]
    # Ensures the maf in bed
    uniq_muts_df = check_mutation_overlap(newbeddf,
                                          merged_mafdf,
                                          threshold = threshold,
                                          threshold_key = "SEQ_ASSAY_ID")
    write.csv(uniq_muts_df,
              paste0(dir, code, "_unique_mutations.csv"),
              row.names = F)
    unique_mutation_files = c(unique_mutation_files,
                              paste0(dir, code, "_unique_mutations.csv"))
  }
  unique_mutation_files
}

if (testing) {
  database_synid_mappingid = 'syn11600968'
} else{
  database_synid_mappingid = 'syn10967259'
}

database_synid_mapping = synTableQuery(sprintf('select * from %s',
                                               database_synid_mappingid))
database_synid_mappingdf = synapser::as.data.frame(database_synid_mapping)
release_folder_fileview_synid = database_synid_mappingdf$Id[database_synid_mappingdf$Database == "releaseFolder"]

choose_from_release = synTableQuery(paste(sprintf("select distinct(name) as releases from %s",
                                                  release_folder_fileview_synid),
                                          "where name not like 'Release%' and name <> 'case_lists'"))
releases = synapser::as.data.frame(choose_from_release)
if (!any(releases$releases %in% release)) {
  stop(sprintf("Must choose correct release: %s",
               paste0(releases$releases, collapse=", ")))
}

release_folder = synTableQuery(sprintf("select id from %s where name = '%s'",
                                       release_folder_fileview_synid, release))
release_folder_synid = release_folder$asDataFrame()$id

# Get release files mapping
RELEASE_FILES_MAPPING = get_file_mapping(release_folder_synid)
NUMBER_PANELS_COVER_REGION = 5

BED_SYNID = RELEASE_FILES_MAPPING[['genomic_information.txt']]
MAF_SYNID = RELEASE_FILES_MAPPING[['data_mutations_extended.txt']]
CLINICAL_SAMPLE_SYNID = RELEASE_FILES_MAPPING[['data_clinical_sample.txt']]
bed_ent = synGet(BED_SYNID, followLink=T)
maf_ent = synGet(MAF_SYNID, followLink=T)
clinical_ent = synGet(CLINICAL_SAMPLE_SYNID, followLink=T)

CLINICALDF = readr::read_tsv(clinical_ent$path, comment="#")
CLINICALDF$CENTER = sapply(strsplit(CLINICALDF$SEQ_ASSAY_ID, "-"), function(x) x[1])
CLINICALDF = CLINICALDF[,c("SAMPLE_ID", "CENTER", "SEQ_ASSAY_ID", "ONCOTREE_CODE")]

MAFDF = readr::read_tsv(maf_ent$path, comment="#")
MAFDF = MAFDF[,c("Hugo_Symbol", "Chromosome", "Start_Position",
                 "End_Position", "Center", "Tumor_Sample_Barcode",
                 "HGVSp_Short")]

BEDDF = readr::read_tsv(bed_ent$path, comment="#")
BEDDF$CENTER = sapply(strsplit(BEDDF$SEQ_ASSAY_ID, "-"), function(x) x[1])

# Get number of oncotree codes per seq assay id
codes_per_panel = table(CLINICALDF$ONCOTREE_CODE, CLINICALDF$SEQ_ASSAY_ID)
# Gets the number of centers with a certain codes
codes_count_across_panels = apply(codes_per_panel, 1, function(x) {
  sum(x > 0)
})
hist(codes_count_across_panels,
     main="Distribution of codes seen across panels",
     xlab = "Codes seen across x panels")
# Based on the histogram, choose 15 as the threshold
# number_center_with_panels = 15
number_panels_with_code =  floor(0.75*length(unique( CLINICALDF$SEQ_ASSAY_ID)))

CODES_ABOVE_THRESHOLD_PANEL = codes_count_across_panels[
  codes_count_across_panels >= number_panels_with_code]


unique_mutation_folder = "~/Genie/analyses/genomicData/unique_muts_panel300/"
dir.create(unique_mutation_folder)
unique_mutation_files = list.files(unique_mutation_folder,
                                   full.names = T)

if (length(unique_mutation_files) == 0) {
  unique_mutation_files = find_unique_mutations_panel(CLINICALDF,
                                                      CODES_ABOVE_THRESHOLD_PANEL,
                                                      MAFDF,
                                                      BEDDF,
                                                      threshold=NUMBER_PANELS_COVER_REGION,
                                                      dir = unique_mutation_folder)
} 


panel_unique_mutation_counts = 
  matrix(nrow = length(unique(CLINICALDF$ONCOTREE_CODE)),
         ncol = length(unique(CLINICALDF$SEQ_ASSAY_ID)),
         dimnames = list(unique(CLINICALDF$ONCOTREE_CODE),
                         unique(CLINICALDF$SEQ_ASSAY_ID)))
# Need to write out the panel id
for (mut_file in unique_mutation_files) {
  mutdf = read.csv(mut_file)
  code = sub(".*/(.*)_unique_mutations.csv","\\1", mut_file)
  code_count = table(mutdf$SEQ_ASSAY_ID)
  panel_unique_mutation_counts[code, names(code_count)] = code_count
}
panel_unique_mutation_counts[is.na(panel_unique_mutation_counts)] = 0


row_sum = apply(panel_unique_mutation_counts, 1, sum)
# barplot(sort(row_sum), las = 2)

# row_above_0 = row_sum[row_sum > 0]
# barplot(sort(row_above_0), las = 2)
# barplot(sort(log(row_above_0)), las = 2)

#Plot codes above average of counts
row_sum_avg = row_sum[row_sum > mean(row_sum)]
barplot(sort(row_sum_avg, decreasing = T),
        main = "Unique number of mutations per code",
        las = 2,
        cex.names = 0.8)

barplot(sort(log(row_sum_avg), decreasing = T),
        main = "Unique number of mutations per code (log)",
        las = 2,
        cex.names = 0.8)
