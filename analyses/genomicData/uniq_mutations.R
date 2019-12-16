library(synapser)
library(VariantAnnotation)
library(knitr)
library(glue)


synLogin()


MAF_SYNID = "syn21348156"
CLINICAL_SAMPLE_SYNID = "syn21348153"
BED_SYNID = "syn21348159"

clinical_ent = synGet(CLINICAL_SAMPLE_SYNID)

clinical_query = glue::glue("select SAMPLE_ID, CENTER, SEQ_ASSAY_ID, ONCOTREE_CODE ",
                            "from {CLINICAL_SAMPLE_TABLE_SYNID}")
codes = synTableQuery(clinical_query, includeRowIdAndRowVersion = F)
CODESDF = codes$asDataFrame()


check_mutation_overlap <- function(beddf, mutdf, threshold=NA) {
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
      
      if (length(unique(beddf$CENTER[overlap])) >= threshold) {
        return(T)
      } else{
        return(F)
      } 
    })
    mutdf[overlap_with_threshold, ]
  } else {
    mutdf[maf_vr %over% bed_gr,]
  }

  #code
  #overlapping = beddf[bed_gr %over% maf_vr,]
  #paste(unique(overlapping$CENTER),collapse = ",")
}


find_unique_mutations_panel <- function(codesdf, codes_above_threshold, dir) {
  unique_mutation_files = c()
  for (code in names(codes_above_threshold)) {
    print(code)
    uniq_muts_df = data.frame()
    # Get oncotree code
    samples = codesdf$SAMPLE_ID[codesdf$ONCOTREE_CODE == code]
    samples_str = paste(samples, collapse = "','")
    query_str = glue::glue("select Hugo_Symbol, Chromosome, Start_Position, ",
                           "End_Position, Center, Tumor_Sample_Barcode from ",
                           "{MAF_TABLE_SYNID} where Tumor_Sample_Barcode in ('{samples_str}')")
    mutation = synTableQuery(query_str, includeRowIdAndRowVersion = F)
    mutationdf = mutation$asDataFrame()
    merged_mafdf = merge.data.frame(mutationdf, codesdf,
                                    by.x = "Tumor_Sample_Barcode",
                                    by.y = "SAMPLE_ID")
    symbol_per_panel = table(merged_mafdf$Hugo_Symbol, merged_mafdf$SEQ_ASSAY_ID)
    unique_mutation_check = apply(symbol_per_panel, 1, function(x) {
      sum(x > 0)
    })
    unique_mutation = unique_mutation_check[unique_mutation_check == 1]
    # Grab bed regions for the panels involved with this mutation
    
    query_string = glue::glue("select * from {BED_TABLE_SYNID} where SEQ_ASSAY_ID in ('",
                              paste(unique(merged_mafdf$SEQ_ASSAY_ID),collapse = "','"),
                              "')")
    bed = synTableQuery(query_string, includeRowIdAndRowVersion = F)
    beddf = bed$asDataFrame()
    mutdf = merged_mafdf[merged_mafdf$Hugo_Symbol %in% names(unique_mutation),]
    uniq_muts_df = check_mutation_overlap(beddf, mutdf)
      
    write.csv(uniq_muts_df,
              paste0(dir, code, "_unique_mutations.csv"),
              row.names = F)
    unique_mutation_files = c(unique_mutation_files,
                              paste0(dir, code, "_unique_mutations.csv"))
  }
  unique_mutation_files
}


codes_per_panel = table(CODESDF$ONCOTREE_CODE, CODESDF$SEQ_ASSAY_ID)
# Gets the number of centers with a certain codes
codes_count_across_panels = apply(codes_per_panel, 1, function(x) {
  sum(x > 0)
})
hist(codes_count_across_panels,
     main="Distribution of codes seen across panels",
     xlab = "Codes seen across x panels")
# Based on the histogram, choose 15 as the threshold
number_center_with_panels = 15

codes_above_threshold = codes_count_across_panels[
  codes_count_across_panels >= number_center_with_panels]


# MUST RERUN THIS
unique_mutation_files = list.files("unique_muts_panel",full.names = T)

if (length(unique_mutation_files) == 0) {
    unique_mutation_files = find_unique_mutations_panel(
      CODESDF,
      codes_above_threshold,
      dir = "unique_muts_panel/")
} 



panel_unique_mutation_counts = 
  matrix(nrow = length(unique(CODESDF$ONCOTREE_CODE)),
         ncol = length(unique(CODESDF$SEQ_ASSAY_ID)),
         dimnames = list(unique(CODESDF$ONCOTREE_CODE),
                         unique(CODESDF$SEQ_ASSAY_ID)))
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