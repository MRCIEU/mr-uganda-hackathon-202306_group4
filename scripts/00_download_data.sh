###################################
## MRC-Uganda MR Hackathon
## -- Downloading the needed data
## by: David A Hughes
## date: Jun 29th 2023
###################################


## Download the SBP data from PanUKBB
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-SBP-both_sexes-auto_medadj_irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files_tabix/continuous-SBP-both_sexes-auto_medadj_irnt.tsv.bgz.tbi

## Download the BMI data from PanUKBB
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-21001-both_sexes-irnt.tsv.bgz
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files_tabix/continuous-21001-both_sexes-irnt.tsv.bgz.tbi

## Download the SBP data from UGR
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009053

## Download the BMI data from UGR
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009057





