python GENIE.py main
python GENIE.py maf --vcf2mafPath ~/vcf2maf-1.6.14 --vepPath ~/vep --vepData ~/.vep --createNewMafDatabase
python GENIE.py vcf --vcf2mafPath ~/vcf2maf-1.6.14 --vepPath ~/vep --vepData ~/.vep

#python GENIE.py mafSP 
#python runAnnovar.py

# python MAFinBED.py > $logDir/MAFinBED.txt 2>&1
# python database_to_staging.py Jul-2017 ~/cbioportal/ > $logDir/database_to_staging.txt 2>&1
# python consortium_to_public.py Jul-2017 ~/cbioportal/ > $logDir/consortium_to_public.txt 2>&1

# #Run dashboard stuff
# python ../sage_dashboard/dashboardTableUpdater.py --release $RELEASE $STAGING > $logDir/dashboard.txt 2>&1
# Rscript ../sage_dashboard/public_dashboard.R "Database" $STAGING > $logDir/dashboardR.txt 2>&1
# Rscript ../sage_dashboard/clinicalImages.R $STAGING > $logDir/clinical.txt 2>&1