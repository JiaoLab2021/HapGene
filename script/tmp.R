countdata <- read.table("featureCounts.txt",skip = 1,sep="\t",header = T,row.names = 1, check.names = FALSE) # skip = 1，跳过文件第一行注释信息，sep = "\t"表示使用制表符作为数据文件中的字段分隔符，header = T表示数据文件的第一行包含列名，row.names = 1表示数据文件的第一列作为行名。
colnames(countdata) <- sapply(strsplit(basename(colnames(countdata)), "\\.bam"), function(x) x[[1]])
# colnames(countdata) <- gsub("-", "_", colnames(counts))
# colnames(countdata) <- sub("^.*_(.*)_bam$", "\\1", colnames(countdata))  # 简化样本名
metadata <- countdata[,1:5]  # 提取基因信息count数据的1-5列
countdata <- countdata[,6:ncol(countdata)]  # 提取counts数，counts数据主体部分（就是将元素分开，只保留counts值部分）
prefix <- "featurecounts"  # 设置输出文件前缀名

# ------ FPKM Calculation------
kb <- metadata$Length / 1000  # 提取Length列的值除以1000，单位是kb，不是bp
fpkm <- t(t(countdata/kb)/colSums(countdata)*1000000) # 
write.csv(fpkm,paste0(prefix,"_fpkm.csv"))

# ------ TPM Calculation------
kb <- metadata$Length / 1000  # 提取Length列的值除以1000，单位是kb，不是bp
rpk <- countdata / kb  # 每千碱基reads，长度标准化
tpm <- t(t(rpk)/colSums(rpk) * 1000000) # 每百万缩放因子，深度标准化
write.csv(tpm,paste0(prefix,"_tpm.csv"))
# avg_tpm <- data.frame(avg_tpm=rowMeans(tpm))
# write.csv(avg_tpm,paste0(prefix,"_avg_tpm.csv"))

# ---------- CPM ----------
cpm <- t(t(countdata)/colSums(countdata) * 1000000) #参考cpm定义
write.csv(cpm,paste0(prefix,"_cpm.csv"))
# avg_cpm <- data.frame(avg_cpm=rowMeans(cpm))
# write.csv(avg_cpm,paste0(prefix,"_avg_cpm.csv"))`
