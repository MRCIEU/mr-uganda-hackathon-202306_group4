###################################
## MRC-Uganda MR Hackathon
## -- Extracting the needed SNPs from
##    the GWAS files
## by: David A Hughes
## date: Jun 29th 2023
###################################

## load R into environment
module add lang/r/4.3.0-gcc


############################
## now working in R
############################
source("parameters.R")

### FILES USED from parameter file

#### YENGO INSTRUMENTS
## The Yengo BMI SNP Instrument file is "yengo_BMI_656.txt"
## The Yengo HEIGHT SNP Instrument file is "yengo_height_AFR_associations.txt"

#### UGR DATA
## The UGR BMI GWAS Summary Stats File is bmiannotated.txt.gz"
## The UGR SNP GWAS Summary Stats File is SBPannotated.txt.gz"

## The PAN-UKB BMI GWAS Summary Stats File is continuous-21001-both_sexes-irnt.tsv"
## The PAN-UKB SBP GWAS Summary Stats File is continuous-SBP-both_sexes-auto_medadj_irnt.tsv"



## Read in the yengo BMI snp data
yengo_data = read.table(yengo, header = TRUE, sep = "\t")
head(yengo_data)

## make a SNPID that looks like that in UGR
yengo_data$snpid = sapply(1:nrow(yengo_data), function(i){
		paste0(yengo_data$CHR[i], ":", yengo_data$POS[i], ":")	
	})

## Reorder the file
o = order(yengo_data$CHR, yengo_data$POS)
yengo_data = yengo_data[o,]

#### Load data.table
library(data.table)

############################################################
##
## Process the UGR BMI data set
##
############################################################
ugr_bmi_data = fread(ugr_bmi, header = TRUE, sep = " ")

start_time <- Sys.time()
chr_bp_a1_a2 = t( sapply(ugr_bmi_data$snpid, function(x){
	strsplit(x, split = ":")[[1]]
	}) )
end_time <- Sys.time()
end_time - start_time

colnames(chr_bp_a1_a2) = c("chr","bp","a1","a2")
chr_bp_a1_a2 = as.data.frame(chr_bp_a1_a2)
chr_bp_a1_a2[,1] = as.numeric(chr_bp_a1_a2[,1])
chr_bp_a1_a2[,2] = as.numeric(chr_bp_a1_a2[,2])


## Find the YENGO SNPs in this data
UGR_BMI_YENGO_data = c()
chrs = sort( unique(yengo_data$CHR) )
for(chr in chrs){
	cat(paste0("now processing chr ", chr, "\n"))
	## YENGO data for a chromsome
	w = which(yengo_data$CHR == chr)
	temp = yengo_data[w, ]
	
	## UGR data for a chromosome
	A = which(chr_bp_a1_a2$chr == chr)
	temp_chr_bp_a1_a2 = chr_bp_a1_a2[A, ]
	temp_ugr_bmi_data = ugr_bmi_data[A, ]

	## find the needed positions
	# B = which(temp_chr_bp_a1_a2$bp %in% temp$POS)
	B = match(temp$POS, temp_chr_bp_a1_a2$bp)
	
	temp_chr_bp_a1_a2 = temp_chr_bp_a1_a2[B, ]
	temp_ugr_bmi_data = temp_ugr_bmi_data[B,]
	temp_ugr_bmi_data = cbind(temp_ugr_bmi_data, temp_chr_bp_a1_a2)
	
	## 
	UGR_BMI_YENGO_data = rbind(UGR_BMI_YENGO_data, temp_ugr_bmi_data)
}

UGR_BMI_YENGO_data$rsid = yengo_data$SNP
## Make a ugr bmi data file with just the yengo SNPs
write.table(UGR_BMI_YENGO_data, file = "UGR_BMI_YENGO_data.txt", 
	col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

### remove data we are done with
rm(UGR_BMI_YENGO_data)
rm(chr_bp_a1_a2)
rm(ugr_bmi_data)

############################################################
##
## Process the UGR SBP data set
##
############################################################
ugr_sbp_data = fread(ugr_sbp, header = TRUE, sep = " ")


start_time <- Sys.time()
chr_bp_a1_a2 = t( sapply(ugr_sbp_data$snpid, function(x){
	strsplit(x, split = ":")[[1]]
	}) )
end_time <- Sys.time()
end_time - start_time

colnames(chr_bp_a1_a2) = c("chr","bp","a1","a2")
chr_bp_a1_a2 = as.data.frame(chr_bp_a1_a2)
chr_bp_a1_a2[,1] = as.numeric(chr_bp_a1_a2[,1])
chr_bp_a1_a2[,2] = as.numeric(chr_bp_a1_a2[,2])



## Find the YENGO SNPs in this data
UGR_SBP_YENGO_data = c()
chrs = sort( unique(yengo_data$CHR) )
for(chr in chrs){
	cat(paste0("now processing chr ", chr, "\n"))
	## YENGO data for a chromsome
	w = which(yengo_data$CHR == chr)
	temp = yengo_data[w, ]
	
	## UGR data for a chromosome
	A = which(chr_bp_a1_a2$chr == chr)
	temp_chr_bp_a1_a2 = chr_bp_a1_a2[A, ]
	temp_ugr_data = ugr_sbp_data[A, ]

	## find the needed positions
	# B = which(temp_chr_bp_a1_a2$bp %in% temp$POS)
	B = match(temp$POS, temp_chr_bp_a1_a2$bp)
	
	temp_chr_bp_a1_a2 = temp_chr_bp_a1_a2[B, ]
	temp_ugr_data = temp_ugr_data[B,]
	temp_ugr_data = cbind(temp_ugr_data, temp_chr_bp_a1_a2)
	
	## 
	UGR_SBP_YENGO_data = rbind(UGR_SBP_YENGO_data, temp_ugr_data)
}

UGR_SBP_YENGO_data$rsid = yengo_data$SNP


## Make a ugr bmi data file with just the yengo SNPs
write.table(UGR_SBP_YENGO_data, file = "UGR_SBP_YENGO_data.txt", 
	col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

rm(UGR_SBP_YENGO_data)
rm(chr_bp_a1_a2)
rm(ugr_sbp_data)


############################################################
##
## Process the PANUKB SBP data set
##
############################################################
pnukb_sbp_data = fread(pnukb_sbp, header = TRUE)


## Find the YENGO SNPs in this data
UKB_SBP_YENGO_data = c()
chrs = sort( unique(yengo_data$CHR) )
for(chr in chrs){
	cat(paste0("now processing chr ", chr, "\n"))
	## YENGO data for a chromsome
	w = which(yengo_data$CHR == chr)
	temp = yengo_data[w, ]
	
	## UGR data for a chromosome
	A = which(pnukb_sbp_data$chr == chr)
	temp_pnukb_sbp_data = pnukb_sbp_data[A, ]
	
	## find the needed positions
	# B = which(temp_chr_bp_a1_a2$bp %in% temp$POS)
	B = match(temp$POS, temp_pnukb_sbp_data$pos)
	temp_pnukb_sbp_data = temp_pnukb_sbp_data[B, ]
	
	## 
	UKB_SBP_YENGO_data = rbind(UKB_SBP_YENGO_data, temp_pnukb_sbp_data)
}


UKB_SBP_YENGO_data$rsid = yengo_data$SNP


## Make a ugr bmi data file with just the yengo SNPs
write.table(UKB_SBP_YENGO_data, file = "PANUKB_SBP_YENGO_data.txt", 
	col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

rm(UKB_SBP_YENGO_data)
rm(pnukb_sbp_data)



############################################################
##
## Process the PANUKB BMI data set
##
############################################################
pnukb_bmi_data = fread(pnukb_bmi, header = TRUE)

## Find the YENGO SNPs in this data
UKB_BMI_YENGO_data = c()
chrs = sort( unique(yengo_data$CHR) )
for(chr in chrs){
	cat(paste0("now processing chr ", chr, "\n"))
	## YENGO data for a chromsome
	w = which(yengo_data$CHR == chr)
	temp = yengo_data[w, ]
	
	## UGR data for a chromosome
	A = which(pnukb_bmi_data$chr == chr)
	temp_pnukb_bmi_data = pnukb_bmi_data[A, ]
	
	## find the needed positions
	# B = which(temp_chr_bp_a1_a2$bp %in% temp$POS)
	B = match(temp$POS, temp_pnukb_bmi_data$pos)
	temp_pnukb_bmi_data = temp_pnukb_bmi_data[B, ]
	
	## 
	UKB_BMI_YENGO_data = rbind(UKB_BMI_YENGO_data, temp_pnukb_bmi_data)
}


UKB_BMI_YENGO_data$rsid = yengo_data$SNP


## Make a ugr bmi data file with just the yengo SNPs
write.table(UKB_BMI_YENGO_data, file = "PANUKB_BMI_YENGO_data.txt", 
	col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

rm(UKB_BMI_YENGO_data)
rm(pnukb_bmi_data)


## Write the ordered Yengo data to file.
write.table(yengo_data[, 1:13], file = "yengo_data.txt", 
	col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)





###########################################
##
## Read in the yengo HEIGHT snp data
##
###########################################
yengo_height = read.table(yengo_height, header = TRUE, sep = "\t")
head(yengo_height)



############################################################
##
## Process the UGR SBP data set
##
############################################################
ugr_sbp_data = fread(ugr_sbp, header = TRUE, sep = " ")


start_time <- Sys.time()
chr_bp_a1_a2 = t( sapply(ugr_sbp_data$snpid, function(x){
	strsplit(x, split = ":")[[1]]
	}) )
end_time <- Sys.time()
end_time - start_time

colnames(chr_bp_a1_a2) = c("chr","bp","a1","a2")
chr_bp_a1_a2 = as.data.frame(chr_bp_a1_a2)
chr_bp_a1_a2[,1] = as.numeric(chr_bp_a1_a2[,1])
chr_bp_a1_a2[,2] = as.numeric(chr_bp_a1_a2[,2])



## Find the YENGO SNPs in this data
UGR_SBP_YENGO_data = c()
chrs = sort( unique(yengo_height$Chr) )
for(chr in chrs){
	cat(paste0("now processing chr ", chr, "\n"))
	## YENGO data for a chromsome
	w = which(yengo_height$Chr == chr)
	temp = yengo_height[w, ]
	
	## UGR data for a chromosome
	A = which(chr_bp_a1_a2$chr == chr)
	temp_chr_bp_a1_a2 = chr_bp_a1_a2[A, ]
	temp_ugr_data = ugr_sbp_data[A, ]

	## find the needed positions
	B = match(temp$BP_HG19, temp_chr_bp_a1_a2$bp)
	
	temp_chr_bp_a1_a2 = temp_chr_bp_a1_a2[B, ]
	temp_ugr_data = temp_ugr_data[B,]
	temp_ugr_data = cbind(temp_ugr_data, temp_chr_bp_a1_a2)
	
	## 
	UGR_SBP_YENGO_data = rbind(UGR_SBP_YENGO_data, temp_ugr_data)
}

UGR_SBP_YENGO_data$rsid = yengo_height$SNP


## Make a ugr bmi data file with just the yengo SNPs
write.table(UGR_SBP_YENGO_data, file = "UGR_SBP_YENGO_HEIGHT_SNPs_data.txt", 
	col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

rm(UGR_SBP_YENGO_data)
rm(chr_bp_a1_a2)
rm(ugr_sbp_data)



############################################################
##
## Process the PANUKB SBP data set
##
############################################################
pnukb_sbp_data = fread(pnukb_sbp, header = TRUE)


## Find the YENGO SNPs in this data
UKB_SBP_YENGO_data = c()
chrs = sort( unique(yengo_height$Chr) )
for(chr in chrs){
	cat(paste0("now processing chr ", chr, "\n"))
	## YENGO data for a chromsome
	w = which(yengo_height$Chr == chr)
	temp = yengo_height[w, ]
	
	## UGR data for a chromosome
	A = which(pnukb_sbp_data$chr == chr)
	temp_pnukb_sbp_data = pnukb_sbp_data[A, ]
	
	## find the needed positions
	B = match(temp$BP_HG19, temp_pnukb_sbp_data$pos)
	temp_pnukb_sbp_data = temp_pnukb_sbp_data[B, ]
	
	## 
	UKB_SBP_YENGO_data = rbind(UKB_SBP_YENGO_data, temp_pnukb_sbp_data)
}


UKB_SBP_YENGO_data$rsid = yengo_height$SNP


## Make a ugr bmi data file with just the yengo SNPs
write.table(UKB_SBP_YENGO_data, file = "PANUKB_SBP_YENGO_HEIGHT_SNPs_data.txt", 
	col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

rm(UKB_SBP_YENGO_data)
rm(pnukb_sbp_data)










