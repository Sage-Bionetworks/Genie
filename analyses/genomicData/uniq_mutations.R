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
release <- args$release
genie_user <- args$syn_user
genie_pass <- args$syn_pass
testing <- args$testing

library(synapser)
library(VariantAnnotation)
library(knitr)
library(glue)
library(data.table)

tryCatch({
  synLogin()
}, error = function(err) {
  synLogin(genie_user, genie_pass)
})

#' Gets the filename to synapse id mapping of a releease
#' 
#' @param release_folder_synid A synapse id
#' @return A named vector \code{c("filename"= "synapseid")}
#' @examples
#' get_file_mapping("syn12345")
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

check_mutation_overlap <- function(beddf_temp, mutdf_temp,
                                   threshold=NA,
                                   threshold_key="CENTER") {
  # Must check on this overlap.....
  # Must actually check if the unique mutation is actually unique and seen by other sites
  # parametrize,  2 of the regions cover... or 5 of the regions cover....
  maf_vr = GRanges(seqnames = Rle(paste0("chr", mutdf_temp$Chromosome)),
                   ranges = IRanges(start = mutdf_temp$Start_Position,
                                    end = mutdf_temp$End_Position))
  seqlevels(maf_vr) = sort(seqlevels(maf_vr))
  
  bed_gr = GRanges(seqnames = Rle(paste0("chr",beddf_temp$Chromosome)),
                   ranges = IRanges(start = beddf_temp$Start_Position,
                                    end = beddf_temp$End_Position))
  seqlevels(bed_gr) = sort(seqlevels(bed_gr))
  # Number of other centers that must at least have this overlap
  if (!is.na(threshold)) {
    overlap_with_threshold = sapply(c(1:nrow(mutdf_temp)), function(numrow) {
      # This is the suppress the 'Found more than one class "DataFrame" in cache; using the first'
      # Error
      overlap = suppressMessages(bed_gr %over% maf_vr[numrow])
      
      if (length(unique(beddf_temp[[threshold_key]][overlap])) >= threshold) {
        return(T)
      } else{
        return(F)
      } 
    })
    mutdf_temp[overlap_with_threshold, ]
  } else {
    mutdf_temp[maf_vr %over% bed_gr,]
  }
}

find_unique_mutations_panel <- function(clinicaldf_panel, codes_above_threshold,
                                        mafdf_panel, beddf_panel, directory, threshold=NA) {
  unique_mutation_files = c()
  for (code in names(codes_above_threshold)) {
    print(code)
    uniq_muts_df = data.frame()
    
    # Get oncotree code
    samples = clinicaldf_panel$SAMPLE_ID[clinicaldf_panel$ONCOTREE_CODE == code]
    mutationdf = mafdf_panel[mafdf_panel$Tumor_Sample_Barcode %in% samples, ]
    mutationdf$mutation = paste(mutationdf$Hugo_Symbol, mutationdf$HGVSp_Short)
    merged_mafdf = merge.data.frame(mutationdf, clinicaldf_panel,
                                    by.x = "Tumor_Sample_Barcode",
                                    by.y = "SAMPLE_ID")
    # Need to check unique mutation across panels
    # The issue with this is what if a center has 2 panels and each of those panels
    # has one mutation
    
    mutation_per_panel = table(merged_mafdf$mutation, merged_mafdf$SEQ_ASSAY_ID)
    
    unique_mutation_check = apply(mutation_per_panel, 1, function(x) {
      sum(x > 0)
    })
    unique_mutation = unique_mutation_check[unique_mutation_check == 1]
    merged_mafdf = merged_mafdf[merged_mafdf$mutation %in% names(unique_mutation), ]
    # Grab bed regions for the panels involved with this mutation
    newbeddf = beddf_panel[beddf_panel$SEQ_ASSAY_ID %in% unique(merged_mafdf$SEQ_ASSAY_ID), ]
    # Ensures the maf in bed
    if (nrow(merged_mafdf) > 0) {
      uniq_muts_df = check_mutation_overlap(newbeddf,
                                            merged_mafdf,
                                            threshold = threshold)
      if (nrow(uniq_muts_df) > 0) {
        new_code = sub("/", '_', code)
        write.csv(uniq_muts_df,
                  paste0(directory, new_code, "_unique_mutations.csv"),
                  row.names = F)
        unique_mutation_files = c(unique_mutation_files,
                                  paste0(directory, new_code, "_unique_mutations.csv"))
      }
    }
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
                                       release_folder_fileview_synid, release),
                               includeRowIdAndRowVersion=F)
release_folder_synid = release_folder$asDataFrame()$id

# Get release files mapping
RELEASE_FILES_MAPPING = get_file_mapping(release_folder_synid)

# Does this number have to be 5 panels from different sites?
NUMBER_CENTERS_COVER_REGION = 5
NUMBER_CENTERS_WITH_CODE = 5
UNIQUE_MUTATIONS_FOLDER_SYNID = "syn18455620"


release_uniq_mutations_ent = synStore(Folder(release,
                                             parentId = UNIQUE_MUTATIONS_FOLDER_SYNID))
center_with_code_folder_ent = synStore(Folder(sprintf("%s_centers_with_code",
                                                      NUMBER_CENTERS_WITH_CODE),
                                              parentId = release_uniq_mutations_ent$properties$id))
uniq_mutation_folder_ent = synStore(Folder(sprintf("regions_covered_by_%s_centers",
                                                   NUMBER_CENTERS_COVER_REGION),
                                           parentId = center_with_code_folder_ent$properties$id))


BED_SYNID = RELEASE_FILES_MAPPING[['genomic_information.txt']]
MAF_SYNID = RELEASE_FILES_MAPPING[['data_mutations_extended.txt']]
CLINICAL_SAMPLE_SYNID = RELEASE_FILES_MAPPING[['data_clinical_sample.txt']]
bed_ent = synGet(BED_SYNID, followLink=T)
maf_ent = synGet(MAF_SYNID, followLink=T)
clinical_ent = synGet(CLINICAL_SAMPLE_SYNID, followLink=T)

CLINICALDF = fread(clinical_ent$path, skip=4)
CLINICALDF$CENTER = sapply(strsplit(CLINICALDF$SEQ_ASSAY_ID, "-"), function(x) x[1])
CLINICALDF = CLINICALDF[,c("SAMPLE_ID", "CENTER", "SEQ_ASSAY_ID", "ONCOTREE_CODE")]

MAFDF = fread(maf_ent$path)
MAFDF = MAFDF[,c("Hugo_Symbol", "Chromosome", "Start_Position",
                 "End_Position", "Center", "Tumor_Sample_Barcode",
                 "HGVSp_Short")]

BEDDF = fread(bed_ent$path)
BEDDF$CENTER = sapply(strsplit(BEDDF$SEQ_ASSAY_ID, "-"), function(x) x[1])

# Get number of oncotree codes per center
# It doesn't make sense to check if a specific code is seen across a n number of panels
# because a center could have 5 panels... So the thresholding should be done on
# a specific code seen across n number of centers
codes_per_center = table(CLINICALDF$ONCOTREE_CODE, CLINICALDF$CENTER)
# Gets the number of centers with a certain codes
codes_count_across_centers = apply(codes_per_center, 1, function(x) {
  sum(x > 0)
})

length(unique(CLINICALDF$CENTER))
# Number of panels that has a specific oncotree code
# number_centers_with_code =  floor(0.75*length(unique(CLINICALDF$SEQ_ASSAY_ID)))


CODES_ABOVE_THRESHOLD_CENTER = codes_count_across_centers[
  codes_count_across_centers >= NUMBER_CENTERS_WITH_CODE]

unique_mutation_folder = sprintf("%s/%s_centers_with_code/%s_centers_cover_region/",
                                 release,
                                 NUMBER_CENTERS_WITH_CODE,
                                 NUMBER_CENTERS_COVER_REGION)
dir.create(unique_mutation_folder, recursive = T)
unique_mutation_files = list.files(unique_mutation_folder,
                                   full.names = T)

if (length(unique_mutation_files) == 0) {
  unique_mutation_files = find_unique_mutations_panel(CLINICALDF,
                                                      CODES_ABOVE_THRESHOLD_CENTER,
                                                      MAFDF,
                                                      BEDDF,
                                                      threshold=NUMBER_CENTERS_COVER_REGION,
                                                      dir = unique_mutation_folder)
}

panel_unique_mutation_counts = 
  matrix(nrow = length(unique(CLINICALDF$ONCOTREE_CODE)),
         ncol = length(unique(CLINICALDF$SEQ_ASSAY_ID)),
         dimnames = list(unique(CLINICALDF$ONCOTREE_CODE),
                         unique(CLINICALDF$SEQ_ASSAY_ID)))
# Need to write out the panel id
for (mut_file in unique_mutation_files) {
  mutdf = fread(mut_file)
  code = sub(".*/(.*)_unique_mutations.csv","\\1", mut_file)
  code_count = table(mutdf$SEQ_ASSAY_ID)
  panel_unique_mutation_counts[code, names(code_count)] = code_count
  synStore(File(mut_file, parentId = uniq_mutation_folder_ent$properties$id))
}
panel_unique_mutation_counts[is.na(panel_unique_mutation_counts)] = 0


row_sum = apply(panel_unique_mutation_counts, 1, sum)

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
