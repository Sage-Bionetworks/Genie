script_dir=$(dirname $0)
if [ ! -d "$script_dir/testing_package" ]; then
  mkdir $script_dir/testing_package
fi

cd $script_dir/testing_package

synapse get -r syn11601335

genie validate data_clinical_supp_SAGE.txt SAGE --test
genie validate data_CNA_SAGE.txt SAGE --test
genie validate data_mutations_extended_SAGE.txt SAGE --test
genie validate GENIE-SAGE-1-1.vcf SAGE --test
genie validate mutationsInCis_filtered_samples.csv SAGE --test
genie validate nonGENIE_data_clinical.txt SAGE --test
genie validate nonGENIE_data_mutations_extended_SAGE.txt SAGE --test
genie validate nonGENIE_SAGE-AKT1.bed SAGE --test
genie validate patientRetraction.csv SAGE --test
genie validate SAGE-PANEL-1.bed SAGE --test
genie validate SAGE_workflow.md SAGE --test
genie validate sampleRetraction.csv SAGE --test
genie validate genie_data_cna_hg19_SAGE.seg SAGE --test
genie validate data_fusions_SAGE.txt SAGE --test
genie validate vital_status.txt SAGE --test
