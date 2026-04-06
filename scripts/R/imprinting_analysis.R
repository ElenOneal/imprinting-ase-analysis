suppressMessages(library(edgeR))
suppressMessages(library(dplyr))
suppressMessages(library(limma))

#t <- parent_specific_expression(c,samples,'DHR14','SOD',0.05)
parent_specific_expression <- function(counts,cpm.filter,n.cpm.filter,samples,parent1name,parent2name,fdr_cutoff){
  #create TRUE/FALSE statement for crossing directions
  is_sample <- sapply(colnames(counts), FUN=function(x) {strsplit(x, fixed=T, split="_")[[1]][2]}) %in% samples
  design <- data.frame(row.names=colnames(counts),mother = ifelse(is_sample,parent1name,parent2name),readSource=ifelse(grepl('mat_',colnames(counts)),'maternal','paternal'),cross=ifelse(is_sample,'1','2'))
  #print(design)
 #create edgeR.design that incorporates cross direction and readSource (ie, maternal or paternal)
  edgeR.design <- model.matrix(~design$cross + design$readSource)
  # print(edgeR.design)
  edgeR <- DGEList(counts=counts, genes=row.names(counts))
  # str(edgeR)
  #normalize maternal and paternal count data
  edgeR$samples$lib.size <- colSums(edgeR$counts)
  print(edgeR$samples$lib.size)
  # cpm.filter <-  10/(min(edgeR$samples$lib.size)/1000000)
#   print(cpm.filter)
  keep <- rowSums(cpm(edgeR)>cpm.filter) >= n.cpm.filter
  counts <- counts[keep,]
  edgeR <- edgeR[keep, , keep.lib.sizes=FALSE]
  edgeR <- calcNormFactors(edgeR,method='TMM')
  print(edgeR$samples$norm.factors)
  #estimate dispersion
  edgeR <- estimateGLMCommonDisp(edgeR, edgeR.design)
  edgeR <- estimateGLMTrendedDisp(edgeR, edgeR.design)
  edgeR <- estimateGLMTagwiseDisp(edgeR, edgeR.design)
  print(paste0("edgeR common dispersion: ", edgeR$common.dispersion))
  #perform likelihood ratio test for differential expression by parental read source
  edgeR.fit <- glmFit(edgeR, edgeR.design)
  edgeR.lrt <- glmLRT(edgeR.fit, coef="design$readSourcepaternal")
  edgeR.fit$AveLogCPM
  #pdf(paste0(cross,"_smear_edgeR.pdf"))
   #dev.off()
  #get DEGS with false discovery rate and benjamini-hochsberg correction
  de <- decideTestsDGE(edgeR.lrt, p=fdr_cutoff, adjust="BH")
  print(summary(de <- decideTestsDGE(edgeR.lrt, p=fdr_cutoff, adjust="BH")))
  #get dataframe of genes annotated by significance of expression: -1 = MEG, 1=PEG, 0=biallelic
  significant_genes <- data.frame(gene_name=row.names(counts), Decision = de)
  #dataframe of normalized counts
  norm.nc <- data.frame(cpm(edgeR$counts),gene_name=rownames(edgeR$counts))
  #create dataframe with mean maternal and paternal expression of genes from raw counts
  imp.exp <- data.frame(gene_name=row.names(counts))
  imp.exp$mat1 <- rowSums(counts[, is_sample & grepl("mat_", colnames(counts))])
  imp.exp$pat1 <- rowSums(counts[, is_sample & grepl("pat_", colnames(counts))])
  imp.exp$mat2 <- rowSums(counts[, !is_sample & grepl("mat_", colnames(counts))])
  imp.exp$pat2 <- rowSums(counts[, !is_sample & grepl("pat_", colnames(counts))])
  imp.exp$total1 <- with(imp.exp,mat1+pat1)
  imp.exp$total2 <- with(imp.exp,mat2+pat2)
  #imp.exp <- imp.exp %>% filter(total1>0&total2>0)
  imp.exp$propmat1 <- imp.exp$mat1/imp.exp$total1
  imp.exp$propmat2 <- imp.exp$mat2/imp.exp$total2
  imp.exp$propmat1 <- replace(imp.exp$propmat1,is.na(imp.exp$propmat1),0)
  imp.exp$propmat2 <- replace(imp.exp$propmat2,is.na(imp.exp$propmat2),0)
  imp.exp$mean_maternal <- rowMeans(imp.exp[, c('propmat1','propmat2')])
  return(list(norm.nc,de,edgeR.lrt,edgeR.fit,edgeR,significant_genes,imp.exp,edgeR$samples$lib.size))
}


