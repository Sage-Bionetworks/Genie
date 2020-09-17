# Update mutation in cis table
uploadToTable <- function(tbl, databaseSynId, subSetSamples, centerMappingDf) {
  # Old samples
  # filterCenters = centerMappingDf$center[centerMappingDf$mutationInCisFilter == TRUE]
  keepCenters = centerMappingDf$center[centerMappingDf$mutationInCisFilter == "OFF"]
  flagCenters = centerMappingDf$center[centerMappingDf$mutationInCisFilter == "FLAG"]
  query_string = sprintf("SELECT * FROM %s where Tumor_Sample_Barcode in ('%s')",
                         databaseSynId,
                         subSetSamples)
  annotated <- synapser::synTableQuery(query_string)
  annotated_df <- synapser::as.data.frame(annotated)
  annotated_df$HGVSp_Short[is.na(annotated_df$HGVSp_Short)] <- ""
  samples = paste(annotated_df$Tumor_Sample_Barcode,
                  "Sp:", annotated_df$HGVSp_Short,
                  "str:", annotated_df$Start_Position,
                  "ref:", annotated_df$Reference_Allele,
                  "tumor:", annotated_df$Tumor_Seq_Allele2)
  # New samples
  #tbl$HGVSp_Short[is.na(tbl$HGVSp_Short)] <- ""
  #This has to be done because HGVSp_Short is a factor
  nodup_tbl <- tbl[!duplicated(tbl),]
  hgvsp = as.character(nodup_tbl$HGVSp_Short)
  hgvsp[is.na(hgvsp)] <- ""
  new_samples = paste(nodup_tbl$Tumor_Sample_Barcode,
                      "Sp:", hgvsp,
                      "str:", nodup_tbl$Start_Position,
                      "ref:", nodup_tbl$Reference_Allele,
                      "tumor:", nodup_tbl$Tumor_Seq_Allele2)
  #Any data that is not the current output but in the database will be changed to FIXED
  if (nrow(annotated_df) > 0 ) {
    # FLAG variants first
    if (any(annotated_df$Center %in% flagCenters)) {
      annotated_df$Flag[annotated_df$Center %in% flagCenters] = "FLAG"
    }
    if (any(annotated_df$Flag[!samples %in% new_samples] == "TOSS")) {
      annotated_df$Flag[!samples %in% new_samples][
        annotated_df$Flag[!samples %in% new_samples] == "TOSS"] = "FIXED"
    } 
    if (any(annotated_df$Center %in% keepCenters)) {
      annotated_df$Flag[annotated_df$Center %in% keepCenters] = "KEEP"
    }
    synapser::synStore(Table(databaseSynId, annotated_df))
  }
  
  #Append any new data
  new_rows = nodup_tbl[!new_samples %in% samples,]
  #Even when nodup_tbl is empty, it can be subsetted, causing one NA row to be uploaded.
  if (nrow(new_rows) > 0 && nrow(nodup_tbl) > 0) {
    #By default everything is labelled TOSS unless the default filter is FALSE
    #NOTE:  If a center is set from default filter TRUE -> FALSE then switches back,
    #       This center's mutationInCis KEEP data will have to be erased.
    #       And processing will have to be run again after mutationsInCis is run on this center
    #       Due to the centers mutationInCis input file.
    new_rows$Flag = "TOSS"
    new_rows$Flag[new_rows$Center %in% keepCenters] = "KEEP"
    new_rows$Flag[new_rows$Center %in% flagCenters] = "FLAG"
    tableToAppend <- synapser::Table(databaseSynId, new_rows)
    synapser::synStore(tableToAppend)
  }
}

# Flag variants to merge
flag_variants_to_merge <- function(genieMutData, genieClinData, samplesToRun, upload = TRUE) {
  
  genieMutData <- data.frame( lapply( genieMutData , factor ))
  # Create factors for clinical data
  
  SAMPLE_ID = genieClinData$SAMPLE_ID[genieClinData$SAMPLE_ID %in% samplesToRun]
  genieClinData = data.frame(SAMPLE_ID)
  genieClinData <- data.frame( lapply( genieClinData , factor ))
  # all inclusive list of samples should be from the clinical data file
  # therefore factor levels for Tumor_Sample_Barcode in the MAF should be set to that of SAMPLE_ID of clinical data table 
  # check that no samples are listed in the MAF that are not listed in the clinical data file
  # reversing the order of the inputs would tell you which samples are submitted that have no entries (no mutations) in the MAF
  if (length(setdiff(levels(genieMutData$Tumor_Sample_Barcode),
                     levels(genieClinData$SAMPLE_ID))) == 0) {
    genieMutData$Tumor_Sample_Barcode = factor(genieMutData$Tumor_Sample_Barcode,
                                               levels = levels(genieClinData$SAMPLE_ID))
  }
  
  # records with count data that preclude a VAF estimate - set VAF to 100% (1/1, alt/depth)
  # genieMutData$t_depth <-  as.numeric(genieMutData$t_depth) DONT DO THIS LINE...
  genieMutData$t_depth <-  as.numeric(levels(genieMutData$t_depth))[genieMutData$t_depth]
  noVAF.idx = which((genieMutData$t_depth == 0) | is.na(genieMutData$t_depth))
  genieMutData$t_alt_count_num = 
    as.numeric(levels(genieMutData$t_alt_count))[genieMutData$t_alt_count]
  #if (length(noVAF.idx) > 0) {
  genieMutData$t_alt_count_num[noVAF.idx] = 1
  genieMutData$t_depth[noVAF.idx] = 1
  #}
  genieMutData$Start_Position <- as.numeric(as.character(genieMutData$Start_Position))
  genieMutData$End_Position <- as.numeric(as.character(genieMutData$End_Position))
  genieMutData$Tumor_Seq_Allele2 <- as.character(genieMutData$Tumor_Seq_Allele2)
  #invalid class “VRanges” object: if 'alt' is 'NA', then 'altDepth' should be 'NA'
  genieMutData$t_alt_count_num[is.na(genieMutData$Tumor_Seq_Allele2)] <- NA
  
  # get VRanges for all variants called in the MAF
  mafVR = VRanges(seqnames = Rle(paste0("chr",genieMutData$Chromosome)),
                  ranges = IRanges(start = genieMutData$Start_Position,
                                   end = genieMutData$End_Position),
                  ref = genieMutData$Reference_Allele,
                  # alt = genieMutData$Tumor_Seq_Allele2,
                  altDepth = genieMutData$t_alt_count_num,
                  totalDepth = genieMutData$t_depth,
                  sampleNames = genieMutData$Tumor_Sample_Barcode)
  seqlevels(mafVR) = sort(seqlevels(mafVR))
  
  # precompute
  vaf = altDepth(mafVR)/totalDepth(mafVR)
  ord = order(mafVR)
  
  # start with empty table
  tbl = genieMutData[1, c("Center","Tumor_Sample_Barcode","Hugo_Symbol",
                          "HGVSp_Short","Variant_Classification","Chromosome",
                          "Start_Position","Reference_Allele","Tumor_Seq_Allele2",
                          "t_alt_count_num","t_depth")]
  tbl = tbl[-1,]
  
  # check for potential variants that may need to be evaluated for merge (cis/trans)
  genieMutData$Tumor_Sample_Barcode <- as.character(genieMutData$Tumor_Sample_Barcode)
  #genieClinData$SAMPLE_ID <- as.character(genieClinData$SAMPLE_ID)
  t = Sys.time()
  samplesRun = c()
  for (i in 1:length(samplesToRun)) {
    
    # get sample indices (in order from pre sort above)
    idx = ord[which(genieMutData$Tumor_Sample_Barcode[ord] == samplesToRun[i])]
    samplesRun = c(samplesRun, samplesToRun[i])
    # get length of idx
    l = length(idx)
    
    # if sample has more than one variant
    if (l > 1) {
      # get differences in BPs of variant sites
      dBP = distance(mafVR[idx[1:(l - 1)]], mafVR[idx[2:(l)]])
      # get difference in VAFs of variants
      dVAF = abs(diff(vaf[idx]))

      # potential matches - criteria of difference in BPs between of > 0 & < 6 bps difference, < 5% VAF difference
      pm = which((dBP > 0) & (dBP < 6) & (dVAF < .05))
      for (m in pm) {
        # calc difference in codon number
        codonDiff = abs(diff(as.numeric(sapply(strsplit(as.character(genieMutData$Protein_position[c(idx[m],idx[m + 1])]),split = "/"),"[",1))))
        if (is.na(codonDiff) | (codonDiff == 1)) {
          tbl = rbind(tbl,genieMutData[c(idx[m],idx[m + 1]),
                                       c("Center","Tumor_Sample_Barcode","Hugo_Symbol",
                                         "HGVSp_Short","Variant_Classification","Chromosome",
                                         "Start_Position","Reference_Allele",
                                         "Tumor_Seq_Allele2","t_alt_count_num","t_depth")])
        }
      }
    }
    if (upload) {
      # time ticker per 100 cases
      if ((i %% 100) == 0) {
        print(i)
        t[2] = Sys.time()
        print(t[2] - t[1])
        t[1] = t[2]
      }
      #HAve to keep track of how far it is in the chunk
      if (nrow(tbl) > tbl_size_limit || i == length(samplesToRun)) {
        uploadToTable(tbl, mutationsInCisSynId, paste(samplesRun, collapse = "','"), centerMappingDf)
        tbl = genieMutData[1,c("Center","Tumor_Sample_Barcode","Hugo_Symbol",
                               "HGVSp_Short","Variant_Classification","Chromosome",
                               "Start_Position","Reference_Allele","Tumor_Seq_Allele2",
                               "t_alt_count_num","t_depth")]
        tbl = tbl[-1,]
        samplesRun = c()
      }
    } else {
      return(tbl)
    }
  }
}

