###load bed files from intersection of DNVs and 'active' enhancers, then determine burden in CHD vs control

#load in selected variants, ensure that each proband only counted once
case_fetal_enhancer <- read.delim("case750_new_fhp.bed", stringsAsFactors = F, header = F)
fetal_element_header <- c("dnv_chr", "dnv_pos", "dnv_pos1", "uniqueVarID", "nearest_gene", "upstream_gene", "downstream_gene", "fetal_element_chr", "fetal_element_start", "fetal_element_end", "fetal_element_ID", "x", "ATAC_chr", "ATAC_start", "ATAC_end", "ATAC_ID", "ATAC_height", "ATAC_strand", "y")
colnames(case_fetal_enhancer) <- fetal_element_header
case_fetal_enhancer <- case_fetal_enhancer[!duplicated(case_fetal_enhancer$uniqueVarID),]
case_fetal_enhancer$Blinded.ID <- sapply(strsplit(case_fetal_enhancer$uniqueVarID, ":"), "[", 1)
case_fetal_enhancer <- case_fetal_enhancer %>% filter(!Blinded.ID=="1-10377") #does not contribute to multihits
case_fetal_enhancer$fetal_element_gene <- sapply(strsplit(case_fetal_enhancer$fetal_element_ID, ";"), "[", 4)
con_fetal_enhancer <- read.delim("con1611_new_fhp.bed", stringsAsFactors = F, header = F) #6294
colnames(con_fetal_enhancer) <- fetal_element_header
con_fetal_enhancer <- con_fetal_enhancer[!duplicated(case_fetal_enhancer$uniqueVarID),]
con_fetal_enhancer$Blinded.ID <- sapply(strsplit(con_fetal_enhancer$uniqueVarID, ":"), "[", 1)
con_fetal_enhancer$fetal_element_gene <- sapply(strsplit(con_fetal_enhancer$fetal_element_ID, ";"), "[", 4)

casedata <- case_fetal_enhancer
casedata$Blinded.ID <- sapply(strsplit(casedata$uniqueVarID, ":"), "[", 1)
casedata$gene_filter <- paste(casedata$Blinded.ID, casedata$nearest_gene)
case_gene_filter <- casedata[!duplicated(casedata$gene_filter),]
case_gene <- count(case_gene_filter, nearest_gene)
casedata$fetal_element_gene <- sapply(strsplit((as.character(casedata$fetal_element_ID)), ";"), "[", 4)
casedata$nc_gene <- paste(casedata$Blinded.ID, casedata$fetal_element_gene)
case_nc_gene_filter <- casedata[!duplicated(casedata$nc_gene),]
case_nc_gene <- count(case_nc_gene_filter, fetal_element_gene)
casedata$nc_ID <- paste(casedata$Blinded.ID, casedata$fetal_element_ID)
case_nc_ID_filter <- casedata[!duplicated(casedata$nc_ID),]
case_nc_ID <- count(case_nc_ID_filter, fetal_element_ID)

condata <- con_fetal_enhancer
condata$Blinded.ID <- sapply(strsplit(condata$uniqueVarID, ":"), "[", 1)
condata$gene_filter <- paste(condata$Blinded.ID, condata$nearest_gene)
con_gene_filter <- condata[!duplicated(condata$gene_filter),]
con_gene <- count(con_gene_filter, nearest_gene)
condata$fetal_element_gene <- sapply(strsplit((as.character(condata$fetal_element_ID)), ";"), "[", 4)
condata$nc_gene <- paste(condata$Blinded.ID, condata$fetal_element_gene)
con_nc_gene_filter <- condata[!duplicated(condata$nc_gene),]
con_nc_gene <- count(con_nc_gene_filter, fetal_element_gene)
condata$nc_ID <- paste(condata$Blinded.ID, condata$fetal_element_ID)
con_nc_ID_filter <- condata[!duplicated(condata$nc_ID),]
con_nc_ID <- count(con_nc_ID_filter, fetal_element_ID)

#compare burden by enhancer
compare_nc <- merge(case_nc_gene, con_nc_gene, by="fetal_element_gene", all.x=T, all.y=T) 
compare_nc[is.na(compare_nc)] <- 0

#now look at multihit genes and identify genes with nominal enrichment
load("~/Box Sync/Seidman/Rprojects/GMKF_ALL/annotated_wgs_31Jul18.rda")

for(j in 1:nrow(compare_nc)){
  tabl2 <- compare_nc[j,]
  tabl2[2,2] <- nrow(gmkf_final_annotated)-tabl2[1,2]
  tabl2[2,3] <- nrow(ssc_final_annotated)-tabl2[1,3]
  compare_nc$pval_2side[j] <- fisher.test(tabl2[,2:3])$p.value
}

x <- compare_nc %>% filter(pval_2side<0.05) #still 33 with either total DNV list or just enhancer DNV list
tabl <- cbind(c(sum(x$n.x), (nrow(case_fetal_enhancer)-sum(x$n.x))), c(sum(x$n.y), (nrow(con_fetal_enhancer)-sum(x$n.y))))
nominal_dnv <- fisher.test(tabl)$p.value 
nom_genes <- x$fetal_element_gene