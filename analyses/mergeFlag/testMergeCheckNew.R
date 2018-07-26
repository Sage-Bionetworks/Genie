# libraries
library(synapseClient)
library(VariantAnnotation)
# functions
uploadToTable <- function(tbl, databaseSynId, subSetSamples) {
  # Old samples
  annotated <- synTableQuery(sprintf("SELECT * FROM %s where Tumor_Sample_Barcode in ('%s')", databaseSynId, subSetSamples))
  annotated_df <- annotated@values
  annotated_df$HGVSp_Short[is.na(annotated_df$HGVSp_Short)] <- ""
  samples = paste(annotated_df$Tumor_Sample_Barcode,"Sp:",annotated_df$HGVSp_Short,"str:",annotated_df$Start_Position, "ref:", annotated_df$Reference_Allele,"tumor:",annotated_df$Tumor_Seq_Allele2)
  
  # New samples
  #tbl$HGVSp_Short[is.na(tbl$HGVSp_Short)] <- ""
  #This has to be done because HGVSp_Short is a factor
  nodup_tbl <- tbl[!duplicated(tbl),]
  hgvsp = as.character(nodup_tbl$HGVSp_Short)
  hgvsp[is.na(hgvsp)] <- ""
  new_samples = paste(nodup_tbl$Tumor_Sample_Barcode,"Sp:",hgvsp,"str:",nodup_tbl$Start_Position, "ref:", nodup_tbl$Reference_Allele,"tumor:",nodup_tbl$Tumor_Seq_Allele2)
  
  #Any data that is not the current output but in the database will be changed to CORRECTED
  if (nrow(annotated_df) >0 ){
    if (any(annotated@values$Flag[!samples %in% new_samples] == "TOSS")) {
      annotated@values$Flag[!samples %in% new_samples][annotated@values$Flag[!samples %in% new_samples] == "TOSS"] = "FIXED"
      synStore(annotated)
    }
  }
  
  #Append any new data
  new_rows = nodup_tbl[!new_samples %in% samples,]
  if (nrow(new_rows) > 0) {
    schema <- synGet(databaseSynId)
    new_rows$Flag = "TOSS"
    tableToAppend <- Table(schema, new_rows)
    table <- synStore(tableToAppend)
  }
}
# login to synapse
synapseLogin()

# limits
variant_limit = 100000
tbl_size_limit = 500


centersTable = synTableQuery('select distinct CENTER from syn7517674')
centers = centersTable@values
centers = c("VICC")
for (center in centers) {
  # read aggregated clinical data
  genieClinTable = synTableQuery(sprintf('select SAMPLE_ID from syn7517674 where CENTER = "%s"', center))
  genieClinData = genieClinTable@values

  
  databaseSynIdMapping = synTableQuery('SELECT * FROM syn10967259')
  databaseSynIdMappingDf = databaseSynIdMapping@values
  mafDatabaseSynId = databaseSynIdMappingDf$Id[databaseSynIdMappingDf$Database == "vcf2maf"]
  mafSampleCountTable = synTableQuery(sprintf('select Tumor_Sample_Barcode, count(Tumor_Sample_Barcode) from %s where Center = "%s" group by Tumor_Sample_Barcode',mafDatabaseSynId, center))
  mafSampleCount = mafSampleCountTable@values
  mafSampleCount = mafSampleCount[mafSampleCount$Tumor_Sample_Barcode %in%genieClinData$SAMPLE_ID,]
  total = 0
  splitBySamples = c()
  samplesToQuery = c()
  for (sample in mafSampleCount$Tumor_Sample_Barcode) {
    total = total + mafSampleCount[mafSampleCount$Tumor_Sample_Barcode == sample,"COUNT(Tumor_Sample_Barcode)"]
    samplesToQuery = c(sample,samplesToQuery)
    if (total > variant_limit || sample == mafSampleCount$Tumor_Sample_Barcode[nrow(mafSampleCount)]) {
      splitBySamples = c(splitBySamples, paste(samplesToQuery,collapse="','"))
      total = 0
      samplesToQuery = c()
    }
  }
  for (querySamples in splitBySamples) {
    samplesToRun = unlist(strsplit(querySamples,"','"))
    #SAMPLE_ID = tail(genieClinData$SAMPLE_ID,200)
    #genieClinData = data.frame(SAMPLE_ID)
    # read aggregated MAF file
    genieMutTable = synTableQuery(sprintf("SELECT Center,Tumor_Sample_Barcode,Hugo_Symbol,HGVSp_Short,Variant_Classification,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2,t_depth,t_alt_count,End_Position,Protein_position FROM %s where Tumor_Sample_Barcode in ('%s')", mafDatabaseSynId, querySamples))
    genieMutData = genieMutTable@values
    genieMutData <- data.frame( lapply( genieMutData , factor ))
    # Create factors for clinical data
    
    SAMPLE_ID = genieClinData$SAMPLE_ID[genieClinData$SAMPLE_ID %in% samplesToRun]
    genieClinData = data.frame(SAMPLE_ID)
    genieClinData <- data.frame( lapply( genieClinData , factor ))
    # all inclusive list of samples should be from the clinical data file
    # therefore factor levels for Tumor_Sample_Barcode in the MAF should be set to that of SAMPLE_ID of clinical data table 
    # check that no samples are listed in the MAF that are not listed in the clinical data file
    # reversing the order of the inputs would tell you which samples are submitted that have no entries (no mutations) in the MAF
    if (length(setdiff(levels(genieMutData$Tumor_Sample_Barcode),levels(genieClinData$SAMPLE_ID)))==0) {
      genieMutData$Tumor_Sample_Barcode = factor(genieMutData$Tumor_Sample_Barcode, levels=levels(genieClinData$SAMPLE_ID))
    }
    
    # records with count data that preclude a VAF estimate - set VAF to 100% (1/1, alt/depth)
    noVAF.idx = which((genieMutData$t_depth==0)|is.na(genieMutData$t_depth))
    genieMutData$t_alt_count_num = as.numeric(levels(genieMutData$t_alt_count))[genieMutData$t_alt_count]
    genieMutData$t_alt_count_num[noVAF.idx] = 1
    genieMutData$t_depth <-  as.numeric(levels(genieMutData$t_depth))[genieMutData$t_depth]
    genieMutData$t_depth[noVAF.idx] = 1
    genieMutData$Start_Position <- as.numeric(as.character(genieMutData$Start_Position))
    genieMutData$End_Position <- as.numeric(as.character(genieMutData$End_Position))
    genieMutData$Tumor_Seq_Allele2 <- as.character(genieMutData$Tumor_Seq_Allele2)
    #invalid class “VRanges” object: if 'alt' is 'NA', then 'altDepth' should be 'NA'
    genieMutData$t_alt_count_num[is.na(genieMutData$Tumor_Seq_Allele2)] <- NA
    
    # get VRanges for all variants called in the MAF
    mafVR = VRanges(seqnames=Rle(paste0("chr",genieMutData$Chromosome)),ranges=IRanges(start=genieMutData$Start_Position,end=genieMutData$End_Position),ref=genieMutData$Reference_Allele,alt=genieMutData$Tumor_Seq_Allele2,altDepth=genieMutData$t_alt_count_num,totalDepth=genieMutData$t_depth,sampleNames=genieMutData$Tumor_Sample_Barcode)
    seqlevels(mafVR) = sort(seqlevels(mafVR))
    
    # precompute
    vaf = altDepth(mafVR)/totalDepth(mafVR)
    ord = order(mafVR)
    
    # start with empty table
    tbl = genieMutData[1,c("Center","Tumor_Sample_Barcode","Hugo_Symbol","HGVSp_Short","Variant_Classification","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","t_alt_count_num","t_depth")]
    tbl = tbl[-1,]
    
    # check for potential variants that may need to be evaluated for merge (cis/trans)
    genieMutData$Tumor_Sample_Barcode <- as.character(genieMutData$Tumor_Sample_Barcode)
    #genieClinData$SAMPLE_ID <- as.character(genieClinData$SAMPLE_ID)
    t = Sys.time()
    samplesRun = c()
    for (i in 1:length(samplesToRun)) {

      # get sample indices (in order from pre sort above)
      idx = ord[which(genieMutData$Tumor_Sample_Barcode[ord]==samplesToRun[i])]
      samplesRun = c(samplesRun, samplesToRun[i])
      # get length of idx
      l = length(idx)
      
      # if sample has more than one variant
      if (l>1) {
        # get differences in BPs of variant sites
        dBP = distance(mafVR[idx[1:(l-1)]],mafVR[idx[2:(l)]])
        # get difference in VAFs of variants
        dVAF = abs(diff(vaf[idx]))
        
        # potential matches - criteria of difference in BPs between of > 0 & < 6 bps difference, < 5% VAF difference
        pm = which((dBP>0) & (dBP<6) & (dVAF<.05))
        for (m in pm) {
          # calc difference in codon number
          codonDiff = abs(diff(as.numeric(sapply(strsplit(as.character(genieMutData$Protein_position[c(idx[m],idx[m+1])]),split="/"),"[",1))))
          if (is.na(codonDiff)|(codonDiff==1)) {
            tbl = rbind(tbl,genieMutData[c(idx[m],idx[m+1]),c("Center","Tumor_Sample_Barcode","Hugo_Symbol","HGVSp_Short","Variant_Classification","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","t_alt_count_num","t_depth")])
          }
        }
      }
      
      # time ticker per 100 cases
      if ((i %% 100)==0) {
        print(i)
        t[2] = Sys.time()
        print(t[2]-t[1])
        t[1] = t[2]
      }
      #HAve to keep track of how far it is in the chunk
      if (nrow(tbl) > tbl_size_limit || i == length(samplesToRun)) {
        uploadToTable(tbl, databaseSynIdMappingDf$Id[databaseSynIdMappingDf$Database == "mutationsInCisTestNew"], paste(samplesRun, collapse="','"))
        tbl = genieMutData[1,c("Center","Tumor_Sample_Barcode","Hugo_Symbol","HGVSp_Short","Variant_Classification","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","t_alt_count_num","t_depth")]
        tbl = tbl[-1,]
        samplesRun = c()
      }
    }
   #write.csv(rbind(annotated_df[is.na(annotated_df$Flag),],new_rows), "Missing_variant_annotation.csv", row.names=F)
  }#FOR LOOP END FOR QUERY SAMPLES
}
