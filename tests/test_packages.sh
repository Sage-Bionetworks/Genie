script_dir=$(dirname $0)
if [ ! -d "$script_dir/testing_package" ]; then
  mkdir $script_dir/testing_package
fi

cd $script_dir/testing_package
synapse get -r syn11601335

genie validate clinical data_clinical_supp_SAGE.txt SAGE --test
genie validate cna data_CNA_SAGE.txt SAGE --test
genie validate maf data_mutations_extended_SAGE.txt SAGE --test
genie validate vcf GENIE-SAGE-1-1.vcf SAGE --test
genie validate mutationsInCis mutationsInCis_filtered_samples.csv SAGE --test
genie validate clinicalSP nonGENIE_data_clinical.txt SAGE --test
genie validate mafSP nonGENIE_data_mutations_extended_SAGE.txt SAGE --test
genie validate bedSP nonGENIE_SAGE-AKT1.bed SAGE --test
genie validate patientRetraction patientRetraction.csv SAGE --test
genie validate bed SAGE-PANEL-1.bed SAGE --test
genie validate md SAGE_workflow.md SAGE --test
genie validate sampleRetraction sampleRetraction.csv SAGE --test
genie validate seg genie_data_cna_hg19_SAGE.seg SAGE --test
genie validate fusions data_fusions_SAGE.txt SAGE --test
genie validate vitalStatus vital_status.txt SAGE --test
genie validate patientCounts patient_counts.txt SAGE --test
