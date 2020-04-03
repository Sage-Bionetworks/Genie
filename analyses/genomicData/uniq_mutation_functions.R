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
  file_mapping
}


#' Checks mutation overlap with GRanges
#'
#' @param beddf_temp A bed dataframe
#' @param mutdf_temp A mutation dataframe
#' @param threshold Threshold
#' @param threshold_key Column to check thresholds, CENTER or SEQ_ASSAY_ID
#' @return A mutation dataframe that overlaps with specified threshold and part of the bed file
#' @examples
#' check_mutation_overlap(beddf, mutdf, threshold=3, threshold_key="CENTER)
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
      # suppressMessages 'Found more than one class "DataFrame" in cache; using the first' Error
      overlap = suppressMessages(bed_gr %over% maf_vr[numrow])
      # Only if a region is covered by `threshold` number of CENTERs or SEQ_ASSAY_IDs
      # return T
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


#' Finds unique mutations on a per code basis
#'
#' @param clinicaldf A clinical dataframe
#' @param mafdf A mutation dataframe
#' @param code An oncotree code
#' @return A mutation dataframe containing unique mutations for that code
#' @examples
#' get_code_unique_mutations(clinicaldf, mafdf, code)
get_code_unique_mutations <- function(clinicaldf, mafdf, code) {
  samples = clinicaldf$SAMPLE_ID[clinicaldf$ONCOTREE_CODE == code]
  mutationdf = mafdf[mafdf$Tumor_Sample_Barcode %in% samples, ]
  mutationdf$mutation = paste(mutationdf$Hugo_Symbol, mutationdf$HGVSp_Short)
  merged_mafdf = merge.data.frame(mutationdf, clinicaldf,
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
  merged_mafdf
}


#' Finds all unique mutations for list of oncotree codes
#'
#' @param clinicaldf_panel A clinical dataframe
#' @param oncotree_codes List of oncotree codes
#' @param mafdf_panel A mutation dataframe
#' @param beddf_panel A bed dataframe
#' @param threshold X number of centers that must cover a bed region
#' @return A list of unique mutation files per oncotree code
#' @examples
#' find_unique_mutations(clindf, codes, mafdf, beddf, threshold=3)
find_unique_mutations <- function(clinicaldf_panel, oncotree_codes,
                                  mafdf_panel, beddf_panel, threshold=NA) {
  # unique_mutation_files = c()
  all_unique_mutationsdf = data.frame()
  for (code in oncotree_codes) {
    print(code)
    merged_mafdf = get_code_unique_mutations(clinicaldf_panel, mafdf_panel, code)
    # Grab bed regions for the panels involved with this mutation
    newbeddf = beddf_panel[beddf_panel$SEQ_ASSAY_ID %in% unique(merged_mafdf$SEQ_ASSAY_ID), ]
    # Ensures the maf in bed
    if (nrow(merged_mafdf) > 0) {
      uniq_muts_df = check_mutation_overlap(newbeddf,
                                            merged_mafdf,
                                            threshold = threshold)
      all_unique_mutationsdf = rbind(all_unique_mutationsdf, uniq_muts_df)
    }
  }
  all_unique_mutationsdf
}

#' Finds oncotree codes that are covered by n number of centers
#'
#' @param clinicaldf Clinical dataframe
#' @param threshold The n number of centers
#' @return List of oncotree codes that meet the criteria
#' @examples
#' find_codes_coveredby_n_centers(clinicaldf, 4)
find_codes_coveredby_n_centers <- function(clinicaldf, threshold) {
  codes_per_center = table(clinicaldf$ONCOTREE_CODE, clinicaldf$CENTER)
  # Gets the number of centers with a certain codes
  code_count = apply(codes_per_center, 1, function(x) {
    sum(x > 0)
  })

  oncotree_codes = code_count[code_count >= threshold]
  names(oncotree_codes)
}


#' Gets release folder synapse id given a release name
#'
#' @param database_synid_mappingid Synapse id of the database synid mapping table
#' @param release GENIE release name (ie. 7.2-public)
#' @return Synapse id of release folder
#' @examples
#' get_release_folder_synid("syn10967259", "7.0-public")
get_release_folder_synid <- function(database_synid_mappingid, release) {
  database_synid_mapping = synTableQuery(glue('select * from {synid}',
                                              synid = database_synid_mappingid))
  database_synid_mappingdf = synapser::as.data.frame(database_synid_mapping)
  release_folder_ind = database_synid_mappingdf$Database == "releaseFolder"
  release_folder_fileview_synid = database_synid_mappingdf$Id[release_folder_ind]

  choose_from_release = synTableQuery(glue("select distinct(name) as releases from {synid} where ",
                                           "name not like 'Release%' and name <> 'case_lists'",
                                           synid = release_folder_fileview_synid))
  releases = synapser::as.data.frame(choose_from_release)
  if (!any(releases$releases %in% release)) {
    stop(glue("Must choose correct release: {releases}",
              releases = paste0(releases$releases, collapse=", ")))
  }

  release_folder = synTableQuery(glue("select id from {synid} where name = '{release}'",
                                      synid = release_folder_fileview_synid,
                                      release = release),
                                 includeRowIdAndRowVersion=F)
  release_folder$asDataFrame()$id
}


#' Writes a centers unique mutations and store them into their staging directory
#'
#' @param df Unique mutation dataframe
#' @param database_synid_mapping Synapse id of the database synid mapping table
#' @param center GENIE center (ie. DUKE)
#' @param release GENIE release (ie. 7.0-public)

#' @return Synapse id of stored file
#' @examples
#' write_and_store_mutations(df, "syn10967259", "DUKE", "7.0-public")
write_and_store_mutations <- function(df, database_synid_mappingid, center) {
  filename = paste0(tempdir(), "/", center, "_unique_mutations.tsv")
  write.table(df, filename, quote = F, sep = "\t", row.names = F)

  database_synid_mapping = synTableQuery(glue('select * from {synid}',
                                              synid = database_synid_mappingid))
  database_synid_mappingdf = synapser::as.data.frame(database_synid_mapping)

  center_mapping_ind = database_synid_mappingdf$Database == "centerMapping"
  center_mapping_synid = database_synid_mappingdf$Id[center_mapping_ind]

  staging_folders = synTableQuery(glue("select stagingSynId from {synid} where center = '{center}'",
                                       synid = center_mapping_synid,
                                       center = center),
                                  includeRowIdAndRowVersion = F)
  ent = synStore(File(filename, parentId = staging_folders$asDataFrame()$stagingSynId))
  ent$id
}
