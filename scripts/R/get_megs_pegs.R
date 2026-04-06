suppressMessages(library(dplyr))
suppressMessages(library(edgeR))

get_imprinted_list <- function(imp.df,cross,reversecross,degs,edgeR.lrt,fdr_cutoff){
  imp.df <- imp.df %>% select(gene_name,propmat1,propmat2,mean_maternal)
  colnames(imp.df) <- c('gene_name',paste0(cross,'_prop_maternal'),paste0(reversecross,'_prop_maternal'),'mean_maternal')
  megs <- degs %>% filter(design.readSourcepaternal==-1)
  pegs <- degs %>% filter(design.readSourcepaternal==1)
  megs <- left_join(megs,imp.df,by='gene_name')
  pegs <- left_join(pegs,imp.df,by='gene_name')
  # sig.imp <- data.frame(degs,gene_name=rownames(degs))
  # sig.imp <- sig.imp %>% filter(design.readSourcepaternal!=0)
  # megs <- data.frame(imp.df$gene_name[decideTestsDGE(edgeR.lrt, p=fdr_cutoff, adjust="BH")==-1])
  # pegs <- data.frame(imp.df$gene_name[decideTestsDGE(edgeR.lrt, p=fdr_cutoff, adjust="BH")==1])
  # colnames(megs) <- 'gene_name'
  # colnames(pegs) <- 'gene_name'
  # megs <- megs %>% filter(gene_name %in% degs$gene_name)
  # pegs <- pegs %>% filter(gene_name %in% degs$gene_name)
  lrt.results <- edgeR.lrt$table
  lrt.results$P_adj <- p.adjust(lrt.results$PValue,method='BH')
  lrt.results$gene_name <- rownames(lrt.results)
  rownames(lrt.results) <- c()
  return(list(imp.df,megs,pegs,lrt.results))
}

# l <- get_imprinted_list(nud.imp.exp,cross,reversecross,pse.nud.significant_genes,pse.nud.edgeR.lrt,fdr_cutoff)
# megs <- l[[2]]
# pegs <- l[[3]]
# # nud.imp.exp <- nud.imp.exp %>% select(gene_name,propmat1,propmat2,mean_maternal)
# # colnames(nud.imp.exp) <- c('gene_name',paste0(cross,'_prop_maternal'),paste0(reversecross,'_prop_maternal'),'mean_maternal')
# # # 
# sig.nud.imp.genes <- data.frame(pse.nud.significant_genes,gene_name=rownames(pse.nud.significant_genes))
# sig.nud.imp.genes <- sig.nud.imp.genes %>% filter(design.readSourcepaternal!=0)
# megs <- sig.nud.imp.genes %>% filter(design.readSourcepaternal==-1) 
# pegs <- sig.nud.imp.genes %>% filter(design.readSourcepaternal==1)
# megs <- left_join(megs,)
# megs2 <- data.frame(nud.imp.exp$gene_name[decideTestsDGE(pse.nud.edgeR.lrt, p=fdr_cutoff, adjust="BH")==-1])
# pegs2 <- data.frame(nud.imp.exp$gene_name[decideTestsDGE(pse.nud.edgeR.lrt, p=fdr_cutoff, adjust="BH")==1])
# # megs <- megs %>% filter(gene_name %in% (pse.nud.significant_genes$gene_name))
# #                         
# #                         %>% filter(design.readSourcepaternal==-1))$gene_name)
# # 
# colnames(megs2) <- 'gene_name'
# colnames(pegs2) <- 'gene_name'
# 
# intersection <- intersect(unlist(list(megs$gene_name)), unlist(list(megs$gene_name)))

# 
# lrt.results$P_adj <- p.adjust(lrt.results$PValue,method='BH')
# lrt.results
# megs <- list(megs,lrt.results,imp.exp,mim_arab_orthologs,gene.coord) %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by='gene_name'), .)
# megs <- list(megs,arab_syn,arab_defline) %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by='arabidopsis_name'), .) %>% dplyr::select(gene_name,contig,start,DHRO22_GMR2_prop_maternal,GMR2_DHRO22_prop_maternal,mean_maternal,logFC,logCPM,LR,PValue,P_adj,arabidopsis_name,synonym,defline)
# 
# megs <- unique(megs)
# 
# pegs <- list(pegs,lrt.results,imp.exp,mim_arab_orthologs,gene.coord) %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by='gene_name'), .)
# pegs <- list(pegs,arab_syn,arab_defline) %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by='arabidopsis_name'), .) %>% dplyr::select(gene_name,contig,start,DHRO22_GMR2_prop_maternal,GMR2_DHRO22_prop_maternal,mean_maternal,logFC,logCPM,LR,PValue,P_adj,arabidopsis_name,synonym,defline)
# pegs <- pegs %>% filter(mean_maternal<(1-0.67))
# pegs <- unique(pegs)
# 
# write.table(megs, file=paste0('edgeR_List_MEGs_',cross,'_annotated.txt'), row.names=F, quote=F, col.names=T,sep='\t')
# write.table(pegs, file=paste0('edgeR_List_PEGs_',cross,'_annotated.txt'), row.names=F, quote=F, col.names=T,sep='\t')