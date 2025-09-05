####################  featureCounts后计算FPKM CPM TPM  ############################
#setwd("/Users/icej/Desktop/")  #设置工作目录
rm(list=ls()) # 删除当前工作环境中的所有对象，避免之前对象造成冲突
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
prefix <- args[2]

countdata <- read.table(input_file,skip = 1,sep="\t",header = T,row.names = 1, check.names = FALSE) # skip = 1，跳过文件第一行注释信息，sep = "\t"表示使用制表符作为数据文件中的字段分隔符，header = T表示数据文件的第一行包含列名，row.names = 1表示数据文件的第一列作为行名。
colnames(countdata) <- sapply(strsplit(basename(colnames(countdata)), "\\.bam"), function(x) x[[1]])
# colnames(countdata) <- gsub("-", "_", colnames(counts))
# colnames(countdata) <- sub("^.*_(.*)_bam$", "\\1", colnames(countdata))  # 简化样本名
metadata <- countdata[,1:5]  # 提取基因信息count数据的1-5列
countdata <- countdata[,6:ncol(countdata)]  # 提取counts数，counts数据主体部分（就是将元素分开，只保留counts值部分）
#prefix <- "featurecounts"  # 设置输出文件前缀名

# ------ FPKM Calculation------
kb <- metadata$Length / 1000  # 提取Length列的值除以1000，单位是kb，不是bp
fpkm <- round(t(t(countdata/kb)/colSums(countdata)*1000000), 2)
write.csv(fpkm,paste0(prefix,"_fpkm.csv"))

# ------ TPM Calculation------
kb <- metadata$Length / 1000  # 提取Length列的值除以1000，单位是kb，不是bp
rpk <- countdata / kb  # 每千碱基reads，长度标准化
tpm <- round(t(t(rpk)/colSums(rpk) * 1000000), 2) # 每百万缩放因子，深度标准化
write.csv(tpm,paste0(prefix,"_tpm.csv"))
