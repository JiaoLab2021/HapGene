#################### FPKM CPM TPM  ############################
rm(list=ls()) 
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
prefix <- args[2]

countdata <- read.table(input_file,skip = 1,sep="\t",header = T,row.names = 1, check.names = FALSE)
colnames(countdata) <- sapply(strsplit(basename(colnames(countdata)), "\\.bam"), function(x) x[[1]])

metadata <- countdata[,1:5]  
countdata <- countdata[,6:ncol(countdata)]  


# ------ FPKM Calculation------
kb <- metadata$Length / 1000  
fpkm <- round(t(t(countdata/kb)/colSums(countdata)*1000000), 2)
write.csv(fpkm,paste0(prefix,"_fpkm.csv"))

# ------ TPM Calculation------
kb <- metadata$Length / 1000 
rpk <- countdata / kb  
tpm <- round(t(t(rpk)/colSums(rpk) * 1000000), 2) 
write.csv(tpm,paste0(prefix,"_tpm.csv"))
