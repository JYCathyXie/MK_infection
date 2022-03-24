library(symphony)
library(singlecellmethods)
library(patchwork)

###BloodAdv_ControlMK####
ct_MK_seu_nwid<-readRDS("/media/user/sdd/Ctrl_dim13_nwid_seurat.rds")
BAd_MK<-readRDS("/media/user/sdd/ABM_ALU_MK_forSymphony.rds")
DefaultAssay(ct_MK_seu_nwid_clstid)<-"RNA"
DefaultAssay(BAd_BM_MK)<-"RNA"
test.anchors<-FindIntegrationAnchors(object.list = list(ct_MK_seu_nwid_clstid,BAd_BM_MK),
                                     dims = 1:20,k.filter = 100)
aaa<-IntegrateData(anchorset = test.anchors,dims = 1:20)
DefaultAssay(aaa)<-"integrated"
aaa<-FindVariableFeatures(aaa,selection.method = "vst",nfeatures = 2000)
aaa<-ScaleData(aaa,verbose = FALSE)
aaa<-RunPCA(aaa,npcs = 20,verbose = FALSE)
aaa<-SCTransform(aaa,vars.to.regress = "percent.mt", verbose = FALSE)
aaa<-RunHarmony(aaa,"orig.ident",assay.use="SCT")
aaa<-AddMetaData(aaa,aaa@active.ident,col.name = "cell_type")

exprs_norm<-aaa@assays$RNA@data
metadata<-as.data.frame(aaa@meta.data[,c(1,22)])
idx_query <- which(metadata$orig.ident=="Ctrl_all")
ref_exp_full <- exprs_norm[, -idx_query]
ref_metadata <- metadata[-idx_query, ]
query_exp <- exprs_norm[, idx_query]
query_metadata <- metadata[idx_query, ]

var_genes <- vargenes_vst(ref_exp_full, 
                          groups = as.character(ref_metadata[['cell_type']]), 
                          topn = 2000)
ref_exp <- ref_exp_full[var_genes, ]
dim(ref_exp)

vargenes_means_sds <- tibble(symbol = var_genes, mean = Matrix::rowMeans(ref_exp))
vargenes_means_sds$stddev <- singlecellmethods::rowSDs(ref_exp, 
                                                       vargenes_means_sds$mean)

ref_exp_scaled <- singlecellmethods::scaleDataWithStats(ref_exp, 
                                                        vargenes_means_sds$mean, 
                                                        vargenes_means_sds$stddev, 1)
set.seed(0)
s <- irlba(ref_exp_scaled, nv = 20)
Z_pca_ref <- diag(s$d) %*% t(s$v) # [pcs by cells]
loadings <- s$u

set.seed(0)
ref_harmObj <- harmony::HarmonyMatrix(
  data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
  meta_data = ref_metadata, ## dataframe with cell labels
  theta = c(2),             ## cluster diversity enforcement
  vars_use = c('cell_type'),    ## variable to integrate out
  nclust = 100,             ## number of clusters in Harmony model
  max.iter.harmony = 20,
  return_object = TRUE,     ## return the full Harmony model object
  do_pca = FALSE            ## don't recompute PCs
)
reference <- symphony::buildReferenceFromHarmonyObj(
  ref_harmObj,            # output object from HarmonyMatrix()
  ref_metadata,           # reference cell metadata
  vargenes_means_sds,     # gene names, means, and std devs for scaling
  loadings,               # genes x PCs matrix
  verbose = TRUE,         # verbose output
  do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
  save_uwot_path = '/media/user/sdd/symphony_output/blad_ABMMK_uwot_model_1')

saveRDS(reference, '/media/user/sdd/symphony_output/blad_ABMMK_reference.rds')

umap_labels <- cbind(ref_metadata, reference$umap$embedding)

p0<-ggplot(umap_labels,aes(x=UMAP1, y=UMAP2,color=cell_type))+
  geom_point(size=0.5)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

query <- mapQuery(query_exp,             # query gene expression (genes x cells)
                 query_metadata,        # query metadata (cells x attributes)
                 reference,             # Symphony reference object
                 do_normalize = FALSE,  # perform log(CP10k) normalization on query
                 do_umap = TRUE)
query <- knnPredict(query, reference, reference$meta_data$cell_type, k = 5)
head(query$meta_data)

reference$meta_data$cell_type_pred_knn <- NA
reference$meta_data$ref_query <- 'reference'
query$meta_data$ref_query <- 'query'

meta_data_combined <- rbind(query$meta_data, reference$meta_data)
umap_combined <- rbind(query$umap, reference$umap$embedding)
umap_combined_labels <- cbind(meta_data_combined, umap_combined)

p1<-ggplot(umap_combined_labels,aes(x=UMAP1, y=UMAP2,color=cell_type))+
  geom_point(size=0.5)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  scale_color_manual(values=c(alpha("#e9842c",1),alpha("#00b813",1),alpha("#00b5ee",1),
                              alpha("#f863df",0),alpha("#b34180",1),alpha("#1597bb",1),
                              alpha("#81b214",1),alpha("#da723c",1),alpha("#f05454",1)))
p2<-ggplot(umap_combined_labels,aes(x=UMAP1, y=UMAP2,color=cell_type))+
  geom_point(size=0.5)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  scale_color_manual(values=c(alpha("#e9842c",1),alpha("#00b813",1),alpha("#00b5ee",1),
                              alpha("#f863df",0),alpha("#b34180",0),alpha("#1597bb",0),
                              alpha("#81b214",0),alpha("#da723c",0),alpha("#f05454",0)))
p3<-ggplot(umap_combined_labels,aes(x=UMAP1, y=UMAP2,color=cell_type))+
  geom_point(size=0.5)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  scale_color_manual(values=c(alpha("#e9842c",0),alpha("#00b813",0),alpha("#00b5ee",0),
                              alpha("#f863df",0),alpha("#b34180",1),alpha("#1597bb",1),
                              alpha("#81b214",1),alpha("#da723c",1),alpha("#f05454",1)))
Bld_MK_sym_plot<-p1+p2+p3
ggsave("/media/user/sdd/symphony_output/Bld_BMMK_sym_umap.pdf",
       Bld_MK_sym_plot,width=15,height=4)

###JCI_MK_Ctrl_MK####
ct_MK_seu_nwid<-readRDS("/media/user/sdd/Ctrl_dim13_nwid_seurat.rds")
JCI_newident<-readRDS("/media/user/sdd/GSE158358_lung_JCI/JCI_newident.rds")
JCI_id<-as.data.frame(JCI_newident@active.ident)
LU_MK_cell<-row.names(subset(JCI_id,JCI_id$`JCI_newident@active.ident`=="Lung_MK"))
BM_MK_cell<-row.names(subset(JCI_id,JCI_id$`JCI_newident@active.ident`=="BM_MK"))
JCI_MK<-subset(JCI_newident,cells=c(BM_MK_cell))

DefaultAssay(ct_MK_seu_nwid_clstid)<-"RNA"
DefaultAssay(JCI_MK)<-"RNA"
test.anchors<-FindIntegrationAnchors(object.list = list(ct_MK_seu_nwid_clstid,JCI_MK),
                                     dims = 1:20,k.filter = 100)
aaa<-IntegrateData(anchorset = test.anchors,dims = 1:20)
DefaultAssay(aaa)<-"integrated"
aaa<-FindVariableFeatures(aaa,selection.method = "vst",nfeatures = 2000)
aaa<-ScaleData(aaa,verbose = FALSE)
aaa<-RunPCA(aaa,npcs = 20,verbose = FALSE)
aaa<-SCTransform(aaa,vars.to.regress = "percent.mt", verbose = FALSE)
aaa<-RunHarmony(aaa,"orig.ident",assay.use="SCT")

exprs_norm<-aaa@assays$RNA@data
metadata<-as.data.frame(aaa@meta.data[,c(1,18)])
idx_query <- which(metadata$orig.ident=="Ctrl_all")
ref_exp_full <- exprs_norm[, -idx_query]
ref_metadata <- metadata[-idx_query, ]
query_exp <- exprs_norm[, idx_query]
query_metadata <- metadata[idx_query, ]

var_genes <- vargenes_vst(ref_exp_full, 
                          groups = as.character(ref_metadata[['cell_type']]), 
                          topn = 2000)
ref_exp <- ref_exp_full[var_genes, ]
dim(ref_exp)

vargenes_means_sds <- tibble(symbol = var_genes, mean = Matrix::rowMeans(ref_exp))
vargenes_means_sds$stddev <- singlecellmethods::rowSDs(ref_exp, 
                                                       vargenes_means_sds$mean)

ref_exp_scaled <- singlecellmethods::scaleDataWithStats(ref_exp, 
                                                        vargenes_means_sds$mean, 
                                                        vargenes_means_sds$stddev, 1)
set.seed(0)
s <- irlba(ref_exp_scaled, nv = 20)
Z_pca_ref <- diag(s$d) %*% t(s$v) # [pcs by cells]
loadings <- s$u

set.seed(0)
ref_harmObj <- harmony::HarmonyMatrix(
  data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
  meta_data = ref_metadata, ## dataframe with cell labels
  theta = c(2),             ## cluster diversity enforcement
  vars_use = c('cell_type'),    ## variable to integrate out
  nclust = 100,             ## number of clusters in Harmony model
  max.iter.harmony = 20,
  return_object = TRUE,     ## return the full Harmony model object
  do_pca = FALSE            ## don't recompute PCs
)
reference <- symphony::buildReferenceFromHarmonyObj(
  ref_harmObj,            # output object from HarmonyMatrix()
  ref_metadata,           # reference cell metadata
  vargenes_means_sds,     # gene names, means, and std devs for scaling
  loadings,               # genes x PCs matrix
  verbose = TRUE,         # verbose output
  do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
  save_uwot_path = '/media/user/sdd/symphony_output/JCI_ABMMK_uwot_model_1')

saveRDS(reference, '/media/user/sdd/symphony_output/JCI_ABMMK_reference.rds')

umap_labels <- cbind(ref_metadata, reference$umap$embedding)

p0<-ggplot(umap_labels,aes(x=UMAP1, y=UMAP2,color=cell_type))+
  geom_point(size=0.5)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

query <- mapQuery(query_exp,             # query gene expression (genes x cells)
                  query_metadata,        # query metadata (cells x attributes)
                  reference,             # Symphony reference object
                  do_normalize = FALSE,  # perform log(CP10k) normalization on query
                  do_umap = TRUE)
query <- knnPredict(query, reference, reference$meta_data$cell_type, k = 5)
head(query$meta_data)

reference$meta_data$cell_type_pred_knn <- NA
reference$meta_data$ref_query <- 'reference'
query$meta_data$ref_query <- 'query'

meta_data_combined <- rbind(query$meta_data, reference$meta_data)
umap_combined <- rbind(query$umap, reference$umap$embedding)
umap_combined_labels <- cbind(meta_data_combined, umap_combined)

p1<-ggplot(umap_combined_labels,aes(x=UMAP1, y=UMAP2,color=cell_type))+
  geom_point(size=0.5)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  scale_color_manual(values=c(alpha("#00b813",1),alpha("#b34180",1),
                              alpha("#1597bb",1),alpha("#81b214",1),alpha("#da723c",1),
                              alpha("#f05454",1)))
p2<-ggplot(umap_combined_labels,aes(x=UMAP1, y=UMAP2,color=cell_type))+
  geom_point(size=0.5)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  scale_color_manual(values=c(alpha("#00b813",1),alpha("#b34180",0),
                              alpha("#1597bb",0),alpha("#81b214",0),alpha("#da723c",0),
                              alpha("#f05454",0)))
p3<-ggplot(umap_combined_labels,aes(x=UMAP1, y=UMAP2,color=cell_type))+
  geom_point(size=0.5)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  scale_color_manual(values=c(alpha("#00b813",0),alpha("#b34180",1),
                              alpha("#1597bb",1),alpha("#81b214",1),alpha("#da723c",1),
                              alpha("#f05454",1)))
JCI_MK_sym_plot<-p1+p2+p3
ggsave("/media/user/sdd/symphony_output/JCI_MK_sym_umap.pdf",
       JCI_MK_sym_plot,width=15,height=4)

###NI_immune_cell####
NI_ct_bm_sub
ct_all_seu
ct_id<-as.data.frame(ct_all_seu@active.ident)
imm_cell<-row.names(subset(ct_id,ct_id$`ct_all_seu@active.ident`%in%c("6.MP","7.Neu","8.Neu","9.Neu","10.Neu",
                                                      "11.Neu","12.B","13.B","14.B","15.Mon", "16.Mon",
                                                      "17.Mon","18.DC","19.T")))
ct_imm_seu<-subset(ct_all_seu,cells=imm_cell)

DefaultAssay(NI_ct_bm_sub)<-"RNA"
DefaultAssay(ct_imm_seu)<-"RNA"
test.anchors<-FindIntegrationAnchors(object.list = list(ct_imm_seu,NI_ct_bm_sub),
                                     dims = 1:20,k.filter = 100)
aaa<-IntegrateData(anchorset = test.anchors,dims = 1:20)
DefaultAssay(aaa)<-"integrated"
aaa<-FindVariableFeatures(aaa,selection.method = "vst",nfeatures = 2000)
aaa<-ScaleData(aaa,verbose = FALSE)
aaa<-RunPCA(aaa,npcs = 20,verbose = FALSE)
aaa<-SCTransform(aaa,vars.to.regress = "percent.mt", verbose = FALSE)
aaa<-RunHarmony(aaa,"orig.ident",assay.use="SCT")
aaa<-AddMetaData(aaa,aaa@active.ident,col.name = "cell_type")

exprs_norm<-aaa@assays$RNA@data
metadata<-as.data.frame(aaa@meta.data[,c(1,20)])
idx_query <- which(metadata$orig.ident=="Ctrl_all")
ref_exp_full <- exprs_norm[, -idx_query]
ref_metadata <- metadata[-idx_query, ]
query_exp <- exprs_norm[, idx_query]
query_metadata <- metadata[idx_query, ]

var_genes <- vargenes_vst(ref_exp_full, 
                          groups = as.character(ref_metadata[['cell_type']]), 
                          topn = 2000)
ref_exp <- ref_exp_full[var_genes, ]
dim(ref_exp)

vargenes_means_sds <- tibble(symbol = var_genes, mean = Matrix::rowMeans(ref_exp))
vargenes_means_sds$stddev <- singlecellmethods::rowSDs(ref_exp, 
                                                       vargenes_means_sds$mean)

ref_exp_scaled <- singlecellmethods::scaleDataWithStats(ref_exp, 
                                                        vargenes_means_sds$mean, 
                                                        vargenes_means_sds$stddev, 1)
set.seed(0)
s <- irlba(ref_exp_scaled, nv = 20)
Z_pca_ref <- diag(s$d) %*% t(s$v) # [pcs by cells]
loadings <- s$u

set.seed(0)
ref_harmObj <- harmony::HarmonyMatrix(
  data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
  meta_data = ref_metadata, ## dataframe with cell labels
  theta = c(2),             ## cluster diversity enforcement
  vars_use = c('cell_type'),    ## variable to integrate out
  nclust = 100,             ## number of clusters in Harmony model
  max.iter.harmony = 20,
  return_object = TRUE,     ## return the full Harmony model object
  do_pca = FALSE            ## don't recompute PCs
)
reference <- symphony::buildReferenceFromHarmonyObj(
  ref_harmObj,            # output object from HarmonyMatrix()
  ref_metadata,           # reference cell metadata
  vargenes_means_sds,     # gene names, means, and std devs for scaling
  loadings,               # genes x PCs matrix
  verbose = TRUE,         # verbose output
  do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
  save_uwot_path = '/media/user/sdd/symphony_output/NI_uwot_model_1')

saveRDS(reference, '/media/user/sdd/symphony_output/NI_reference.rds')

umap_labels <- cbind(ref_metadata, reference$umap$embedding)

p0<-ggplot(umap_labels,aes(x=UMAP1, y=UMAP2,color=cell_type))+
  geom_point(size=0.5)+
  xlab("UMAP_1")+
  ylab("UMAP_2")

query <- mapQuery(query_exp,             # query gene expression (genes x cells)
                  query_metadata,        # query metadata (cells x attributes)
                  reference,             # Symphony reference object
                  do_normalize = FALSE,  # perform log(CP10k) normalization on query
                  do_umap = TRUE)
query <- knnPredict(query, reference, reference$meta_data$cell_type, k = 5)

reference$meta_data$cell_type_pred_knn <- NA
reference$meta_data$ref_query <- 'reference'
query$meta_data$ref_query <- 'query'

meta_data_combined <- rbind(query$meta_data, reference$meta_data)
umap_combined <- rbind(query$umap, reference$umap$embedding)
umap_combined_labels <- cbind(meta_data_combined, umap_combined)

p1<-ggplot(umap_combined_labels,aes(x=UMAP1, y=UMAP2,color=cell_type))+
  geom_point(size=0.5)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  scale_color_manual(values=c(alpha("#325288",1),alpha("#325288",1),
                              alpha("#310b0b",1),alpha("#310b0b",1),alpha("#310b0b",1),
                              alpha("#ff8303",1),alpha("#ff8303",1),alpha("#ff8303",1),
                              alpha("#c67ace",1),alpha("#1b2021",1),
                              alpha("#81b214",1),alpha("#325288",1),alpha("#325288",1),
                              alpha("#325288",1),alpha("#310b0b",1),alpha("#c67ace",1),
                              alpha("#ff8303",1),alpha("#81b214",1),alpha("#325288",1),
                              alpha("#1b2021",1)))
p2<-ggplot(umap_combined_labels,aes(x=UMAP1, y=UMAP2,color=cell_type))+
  geom_point(size=0.5)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  scale_color_manual(values=c(alpha("#325288",0),alpha("#325288",0),
                              alpha("#310b0b",0),alpha("#310b0b",0),alpha("#310b0b",0),
                              alpha("#ff8303",0),alpha("#ff8303",0),alpha("#ff8303",0),
                              alpha("#c67ace",0),alpha("#1b2021",0),
                              alpha("#81b214",0),alpha("#325288",0),alpha("#325288",0),
                              alpha("#325288",0),alpha("#310b0b",1),alpha("#c67ace",1),
                              alpha("#ff8303",1),alpha("#81b214",1),alpha("#325288",1),
                              alpha("#1b2021",1)))
p3<-ggplot(umap_combined_labels,aes(x=UMAP1, y=UMAP2,color=cell_type))+
  geom_point(size=0.5)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  scale_color_manual(values=c(alpha("#325288",1),alpha("#325288",1),
                              alpha("#310b0b",1),alpha("#310b0b",1),alpha("#310b0b",1),
                              alpha("#ff8303",1),alpha("#ff8303",1),alpha("#ff8303",1),
                              alpha("#c67ace",1),alpha("#1b2021",1),
                              alpha("#81b214",1),alpha("#325288",1),alpha("#325288",1),
                              alpha("#325288",1),alpha("#310b0b",0),alpha("#c67ace",0),
                              alpha("#ff8303",0),alpha("#81b214",0),alpha("#325288",0),
                              alpha("#1b2021",0)))
NI_imm_sym_plot<-p1+p2+p3
ggsave("/media/user/sdd/symphony_output/NI_imm_sym_umap.pdf",
       NI_imm_sym_plot,width=15,height=4)
