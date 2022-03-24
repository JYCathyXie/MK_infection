library(MetaNeighbor)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(Matrix)
library(Seurat)

all_mx<-as.matrix(ct_all_seu_rn@assays$RNA@data)
all_se <- SummarizedExperiment(list(gene_matrix=all_mx))
all_list_df<-DataFrame(ct_all_seu_rn@meta.data[,c(1,11)])
all_list_df[,"sample_id"]<-row.names(all_list_df)
colnames(all_list_df)<-c("study_id","cell_type","sample_id")
all_list_df$study_id<-as.character(all_list_df$study_id)
all_list_df$cell_type<-as.character(all_list_df$cell_type)
colData(all_se)<-all_list_df

Baso_mx<-as.matrix(Baso_only_seu@assays$RNA@data)
Baso_se <- SummarizedExperiment(list(gene_matrix=Baso_mx))
Baso_list_df<-DataFrame(Baso_only_seu@meta.data[,c(1,7)])
Baso_list_df[,"sample_id"]<-row.names(Baso_list_df)
colnames(Baso_list_df)<-c("study_id","cell_type","sample_id")
colData(Baso_se)<-Baso_list_df

BAd_MK_mx<-as.matrix(BAd_BM_MK@assays$RNA@data)
BAd_MK_se <- SummarizedExperiment(list(gene_matrix=BAd_MK_mx))
BAd_MK_list_df<-DataFrame(BAd_BM_MK@meta.data[,c(1,14)])
BAd_MK_list_df[,"sample_id"]<-row.names(BAd_MK_list_df)
colnames(BAd_MK_list_df)<-c("study_id","cell_type","sample_id")
colData(BAd_MK_se)<-BAd_MK_list_df

JCI_MK_mx<-as.matrix(JCI_BMMK@assays$RNA@data)
JCI_MK_se <- SummarizedExperiment(list(gene_matrix=JCI_MK_mx))
JCI_MK_list_df<-DataFrame(JCI_BMMK@meta.data[,c(1,10)])
JCI_MK_list_df[,"sample_id"]<-row.names(JCI_MK_list_df)
colnames(JCI_MK_list_df)<-c("study_id","cell_type","sample_id")
JCI_MK_list_df$cell_type<-as.character(JCI_MK_list_df$cell_type)
colData(JCI_MK_se)<-JCI_MK_list_df

NI_ct_mx<-as.matrix(NI_ct_bm_sub@assays$RNA@data)
NI_ct_se <- SummarizedExperiment(list(gene_matrix=NI_ct_mx))
NI_ct_list_df<-DataFrame(NI_ct_bm_sub@meta.data[,c(1,14)])
NI_ct_list_df[,"sample_id"]<-row.names(NI_ct_list_df)
colnames(NI_ct_list_df)<-c("study_id","cell_type","sample_id")
NI_ct_list_df$study_id<-as.character(NI_ct_list_df$study_id)
NI_ct_list_df$cell_type<-as.character(NI_ct_list_df$cell_type)
NI_ct_list_df$study_id<-rep("NI_ct",8141)
colData(NI_ct_se)<-NI_ct_list_df

MF_mx<-as.matrix(Mc_seu@assays$RNA@data)
MF_se <- SummarizedExperiment(list(gene_matrix=MF_mx))
MF_list_df<-DataFrame(Mc_seu@meta.data[,c(1,5)])
MF_list_df[,"sample_id"]<-row.names(MF_list_df)
colnames(MF_list_df)<-c("study_id","cell_type","sample_id")
MF_list_df$study_id<-as.character(MF_list_df$study_id)
MF_list_df$cell_type<-as.character(MF_list_df$cell_type)
colData(MF_se)<-MF_list_df

common_genes1 <- intersect(rownames(all_se), rownames(NI_ct_se))
common_genes2 <- intersect(common_genes1, rownames(BAd_MK_se))
common_genes3 <- intersect(common_genes2, rownames(JCI_MK_se))
common_genes4 <- intersect(common_genes3, rownames(Baso_se))
common_genes5 <- intersect(common_genes4, rownames(MF_se))
JCI_MK_se_c<-JCI_MK_se[common_genes5,]
BAd_MK_se_c<-BAd_MK_se[common_genes5,]
all_se_c<-all_se[common_genes5,]
NI_ct_se_c<-NI_ct_se[common_genes5,]
Baso_se_c<-Baso_se[common_genes5,]
MF_se_c<-MF_se[common_genes5,]

new_colData = DataFrame(
  study_id = rep("cells",14624),
  cell_type = c(as.character(BAd_MK_se_c$cell_type),
                as.character(JCI_MK_se_c$cell_type),
                as.character(NI_ct_se_c$cell_type),
                as.character(Baso_se_c$cell_type),
                as.character(MF_se_c$cell_type),
                as.character(all_se_c$cell_type)))
nw_mx<-cbind(assay(BAd_MK_se_c, 1),assay(JCI_MK_se_c, 1),assay(NI_ct_se_c, 1),
             assay(Baso_se_c, 1),assay(MF_se_c, 1),assay(all_se_c, 1))
MK_mix6 <- SingleCellExperiment(
  list(gene_matrix=nw_mx),
  colData = new_colData)
dim(MK_mix6)
var_genes <- variableGenes(dat = MK_mix6, exp_labels = MK_mix6$study_id)
length(var_genes)
celltype_NV <- MetaNeighborUS(var_genes = common_genes5,
                              dat = MK_mix6,
                              study_id = MK_mix6$study_id,
                              cell_type = MK_mix6$cell_type,
                              fast_version = T)

write.csv(celltype_NV,"/media/user/sdd/output.csv")