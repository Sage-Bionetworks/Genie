# libraries
library(synapser)
#pyExec("syn.table_query_timeout=50000")
library(VariantAnnotation)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("Must supply a boolean value")
}
# SAGE login
tryCatch({
  synLogin()# set user and password
}, error = function(err) {
  genieUser = Sys.getenv("GENIE_USER")
  geniePass = Sys.getenv("GENIE_PASS")
  synLogin(genieUser, geniePass)
})
testing = as.logical(args[1])
if (testing) {
  databaseSynIdMappingId = 'syn11600968'
} else {
  databaseSynIdMappingId = 'syn10967259'
}
databaseSynIdMapping = synTableQuery(sprintf('select * from %s', databaseSynIdMappingId),includeRowIdAndRowVersion=F)
databaseSynIdMappingDf = synapser::as.data.frame(databaseSynIdMapping)
patientSynId = databaseSynIdMappingDf$Id[databaseSynIdMappingDf$Database == "patient"]
sampleSynId = databaseSynIdMappingDf$Id[databaseSynIdMappingDf$Database == "sample"]
bedSynId = databaseSynIdMappingDf$Id[databaseSynIdMappingDf$Database == "bed"]
mafSynId = databaseSynIdMappingDf$Id[databaseSynIdMappingDf$Database == "vcf2maf"]
# read aggregated clinical data from tables
patient = synTableQuery(sprintf('SELECT * FROM %s', patientSynId),includeRowIdAndRowVersion=F)
sample = synTableQuery(sprintf('SELECT * FROM %s', sampleSynId),includeRowIdAndRowVersion=F)

patientData = synapser::as.data.frame(patient)
sampleData = synapser::as.data.frame(sample)
sampleData$AGE_AT_SEQ_REPORT_NUMERICAL <- NULL
patientData$BIRTH_YEAR_NUMERICAL <- NULL
patientData$CENTER <- NULL
genieClinData <- merge.data.frame(patientData, sampleData, by="PATIENT_ID")
# read aggregated BED file data
genieBed = synTableQuery(sprintf('SELECT * FROM %s', bedSynId),includeRowIdAndRowVersion=F)
genieBedData = synapser::as.data.frame(genieBed)

# read aggregated MAF file
# FILTERED OUT COMMON VARIANTS HERE
#genieMut = synTableQuery(sprintf('SELECT * FROM %s where FILTER <> "common_variant"', mafSynId))
genieMut = synTableQuery(sprintf('SELECT * FROM %s', mafSynId))
genieMutData = synapser::as.data.frame(genieMut)
#Only use samples that exist in the clinical sample data pool
genieMutData <- genieMutData[genieMutData$Tumor_Sample_Barcode %in% genieClinData$SAMPLE_ID,]
print(nrow(genieMutData))


originalCols = colnames(genieMutData)
# records with count data that preclude a VAF estimate - set VAF to 100% (1/1, alt/depth)    genieMutData$t_depth <-  as.numeric(genieMutData$t_depth)
genieMutData$t_depth <-  as.numeric(genieMutData$t_depth)
noVAF.idx = which((genieMutData$t_depth==0)|is.na(genieMutData$t_depth))
#keeps the order if factors exist
#genieMutData$t_alt_count_num = as.numeric(levels(genieMutData$t_alt_count))[genieMutData$t_alt_count]
genieMutData$t_alt_count_num = as.numeric(genieMutData$t_alt_count)
genieMutData$t_alt_count_num[noVAF.idx] = 1
genieMutData$t_depth_new = genieMutData$t_depth
genieMutData$t_depth_new[noVAF.idx] = 1

#genieMutData$t_depth_new[is.na(genieMutData$t_alt_count_num)] = 1
#genieMutData$t_alt_count_num[is.na(genieMutData$t_alt_count_num)] = 1
# get VRanges for all variants called in the MAF
genieMutData$t_alt_count_num[is.na(genieMutData$Tumor_Seq_Allele2)] = NA

mafVR = VRanges(seqnames=Rle(paste0("chr",genieMutData$Chromosome)),ranges=IRanges(start=genieMutData$Start_Position,end=genieMutData$End_Position),ref=genieMutData$Reference_Allele,alt=genieMutData$Tumor_Seq_Allele2,altDepth=genieMutData$t_alt_count_num,totalDepth=genieMutData$t_depth_new,sampleNames=genieMutData$Tumor_Sample_Barcode)
seqlevels(mafVR) = sort(seqlevels(mafVR))

# get GRanges for each bed
bedGR = split(genieBedData, factor(genieBedData$SEQ_ASSAY_ID))
bedGR = lapply(bedGR, function(x) {
  GR = GRanges(seqnames=Rle(paste0("chr",x$Chromosome)),ranges=IRanges(start=x$Start_Position,end=x$End_Position))
  seqlevels(GR) = sort(seqlevels(GR))
  return(GR)
})

# all inclusive list of samples should be from the clinical data file
# therefore factor levels for Tumor_Sample_Barcode in the MAF should be set to that of SAMPLE_ID of clinical data table 
# check that no samples are listed in the MAF that are not listed in the clinical data file
# reversing the order of the inputs would tell you which samples are submitted that have no entries (no mutations) in the MAF
#if (length(setdiff(levels(genieMutData$Tumor_Sample_Barcode),levels(genieClinData$SAMPLE_ID)))==0) {
#  genieMutData$Tumor_Sample_Barcode = factor(genieMutData$Tumor_Sample_Barcode, levels=levels(genieClinData$SAMPLE_ID))
#}

# add factor to MAF and set to FALSE, forcing the variant to bed match for TRUE below to clear filter
oldInBed = as.logical(genieMutData$inBED)
genieMutData$inBED = FALSE
# collect some stats in panStats and set inBED to TRUE if variant is in corresponding BED
seq_assays = unique(genieClinData$SEQ_ASSAY_ID[!is.na(genieClinData$SEQ_ASSAY_ID)])
panStats = data.frame()
for (panelName in seq_assays) {
  print(panelName)
  samples.idx = which(genieClinData$SEQ_ASSAY_ID==panelName)
  panStats[panelName,"samples"] = length(samples.idx)
  samples.idx = which(genieMutData$Tumor_Sample_Barcode %in% genieClinData$SAMPLE_ID[samples.idx])
  panStats[panelName,"samples with calls"] = length(unique(genieMutData$Tumor_Sample_Barcode[samples.idx]))
  panStats[panelName,"variants calls"] = length(samples.idx)
  panStats[panelName,"unique variants calls"] = length(unique(mafVR[samples.idx]))
  if (length(samples.idx) >0 && !is.null(bedGR[[panelName]])) {
    genieMutData$inBED[samples.idx] = mafVR[samples.idx] %over% bedGR[[panelName]]
  }
  panStats[panelName,"variants called out of BED"] = length(samples.idx[!genieMutData$inBED[samples.idx]])
  panStats[panelName,"unique variants called out of BED"] = length(unique(mafVR[samples.idx[!genieMutData$inBED[samples.idx]]]))
}
genieMutData$t_depth_new <- NULL
genieMutData$t_alt_count_num <- NULL
#Compare old inBED with new inBED column
#If there are differences, update only the diffs
updateMutData = genieMutData[genieMutData$inBED != oldInBed,]

if (nrow(updateMutData) > 0) { 
  synStore(Table(mafSynId, updateMutData[c("ROW_ID","ROW_VERSION","inBED")]))
}