library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(Seurat)
library(DoubletDecon)
library(clusterProfiler)
library(doParallel)
library(parallel)
library(doSNOW)
library(BiocGenerics)
library(mygene)
library(gplots)

Ctrl_all<-Read10X(data.dir = "/media/user/sdd/MKs_10x/2.cellranger.out/Mk_Ctrl/Ctrl/outs/filtered_feature_bc_matrix")
lm_all<-Read10X(data.dir = "/media/user/sdd/MKs_10x/2.cellranger.out/Mk_LM/LM/outs/filtered_feature_bc_matrix")
Ctrl_seu<-CreateSeuratObject(counts = Ctrl_all,project = "1.Ctrl",min.cells = 3,min.features = 200)
Ctrl_seu[["percent.mt"]]<-PercentageFeatureSet(Ctrl_seu,pattern = "^mt-")
Ctrl_seu<-subset(Ctrl_seu,subset = nFeature_RNA>200 & percent.mt<12)
Ctrl_seu<-NormalizeData(Ctrl_seu,verbose = FALSE)
Ctrl_seu<-FindVariableFeatures(Ctrl_seu,selection.method = "vst",nfeatures = 2000)
Vlnplot<-VlnPlot(Ctrl_seu,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3,
                 pt.size = 0.5)
Ctrl_seu<-ScaleData(Ctrl_seu,verbose = FALSE)
Ctrl_seu<-RunPCA(Ctrl_seu,npcs = 20,verbose = FALSE)
Ctrl_seu<-RunTSNE(Ctrl_seu,dims = 1:20)
Ctrl_seu<-RunUMAP(Ctrl_seu,reduction = "pca",dims = 1:20)
Ctrl_seu<-FindNeighbors(Ctrl_seu,reduction = "pca",dims = 1:20)
Ctrl_seu<-FindClusters(Ctrl_seu,resolution = 0.5)
Ctrl_seu_d<-Improved_Seurat_Pre_Process(Ctrl_seu, num_genes=50, write_files=FALSE)
location<-"/media/user/sdd/"
filename<-"ctrl_mt12"
write.table(Ctrl_seu_d$newExpressionFile, paste0(location, filename, "_expression"), sep="\t")
write.table(Ctrl_seu_d$newFullExpressionFile, paste0(location, filename, "_fullExpression"), sep="\t")
write.table(Ctrl_seu_d$newGroupsFile, paste0(location, filename , "_groups"), sep="\t", col.names = F)
ctrl_results<-Main_Doublet_Decon(rawDataFile=Ctrl_seu_d$newExpressionFile, 
                            groupsFile=Ctrl_seu_d$newGroupsFile, 
                            filename=filename, 
                            location=location,
                            fullDataFile=NULL, 
                            removeCC=FALSE, 
                            species="mmu", 
                            rhop=0.25, 
                            write=TRUE, 
                            PMF=TRUE, 
                            useFull=FALSE, 
                            heatmap=FALSE,
                            centroids=TRUE,
                            num_doubs=100, 
                            only50=FALSE,
                            min_uniq=4,nCores = 1)
LIST<-row.names(ctrl_results$Final_nondoublets_groups)
head(LIST)
LIST=gsub('[.]','-',LIST)
ctrl_remove<-subset(x = Ctrl_seu, cells=LIST)
ctrl_remove<-ScaleData(ctrl_remove,verbose = FALSE)
ctrl_remove<-RunPCA(ctrl_remove,npcs = 20,verbose = FALSE)
ctrl_remove<-RunTSNE(ctrl_remove,dims = 1:20)
ctrl_remove<-RunUMAP(ctrl_remove,reduction = "pca",dims = 1:20)
ctrl_remove<-FindNeighbors(ctrl_remove,reduction = "pca",dims = 1:20)
ctrl_remove<-FindClusters(ctrl_remove,resolution = 0.2)
p1<-DimPlot(ctrl_remove,reduction = "umap", label = TRUE)
Ctrl_markers<-FindAllMarkers(object = ctrl_remove,only.pos = TRUE, 
                           min.pct = 0.1, logfc.threshold = 0.1)

lm_seu<-CreateSeuratObject(counts = lm_all,project = "lm_all",min.cells = 3,min.features = 200)
lm_seu[["percent.mt"]]<-PercentageFeatureSet(lm_seu,pattern = "^mt-")
Vlnplot<-VlnPlot(lm_seu,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3,
                 pt.size = 0.5)
lm_seu<-subset(lm_seu,subset = nFeature_RNA>200& percent.mt<12)
lm_seu<-NormalizeData(lm_seu,verbose = FALSE)
lm_seu<-FindVariableFeatures(lm_seu,selection.method = "vst",nfeatures = 2000)
Vlnplot<-VlnPlot(lm_seu,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3,
                 pt.size = 0.5)
lm_seu<-ScaleData(lm_seu,verbose = FALSE)
lm_seu<-RunPCA(lm_seu,npcs = 20,verbose = FALSE)
lm_seu<-RunTSNE(lm_seu,dims = 1:20)
lm_seu<-RunUMAP(lm_seu,reduction = "pca",dims = 1:20)
lm_seu<-FindNeighbors(lm_seu,reduction = "pca",dims = 1:20)
lm_seu<-FindClusters(lm_seu,resolution = 0.5)
lm_seu_d<-Improved_Seurat_Pre_Process(lm_seu, num_genes=50, write_files=FALSE)
location<-"/media/user/sdd/"
filename<-"lm_mt12"
write.table(lm_seu_d$newExpressionFile, paste0(location, filename, "_expression"), sep="\t")
write.table(lm_seu_d$newFullExpressionFile, paste0(location, filename, "_fullExpression"), sep="\t")
write.table(lm_seu_d$newGroupsFile, paste0(location, filename , "_groups"), sep="\t", col.names = F)

lm_results<-Main_Doublet_Decon(rawDataFile=lm_seu_d$newExpressionFile, 
                                 groupsFile=lm_seu_d$newGroupsFile, 
                                 filename=filename, 
                                 location=location,
                                 fullDataFile=NULL, 
                                 removeCC=FALSE, 
                                 species="mmu", 
                                 rhop=1, 
                                 write=TRUE, 
                                 PMF=TRUE, 
                                 useFull=FALSE, 
                                 heatmap=FALSE,
                                 centroids=TRUE,
                                 num_doubs=100, 
                                 only50=FALSE,
                                 min_uniq=4,nCores = 1)
LIST<-row.names(lm_results$Final_nondoublets_groups)
head(LIST)
LIST=gsub('[.]','-',LIST)
lm_remove<-subset(x = lm_seu, cells=LIST)

test.anchors<-FindIntegrationAnchors(object.list = list(ctrl_remove,lm_remove),
                                     dims = 1:20,k.filter = 100)
mt12_dedou_combined<-IntegrateData(anchorset = test.anchors,dims = 1:15)
DefaultAssay(mt12_dedou_combined)<-"integrated"
mt12_dedou_combined<-FindVariableFeatures(mt12_dedou_combined,selection.method = "vst",nfeatures = 2000)
mt12_dedou_combined<-ScaleData(mt12_dedou_combined,verbose = FALSE)
mt12_dedou_combined<-RunPCA(mt12_dedou_combined,npcs = 20,verbose = FALSE)
mt12_dedou_combined<-SCTransform(mt12_dedou_combined, 
                                 vars.to.regress = "percent.mt", verbose = FALSE)
mt12_dedou_combined<-RunHarmony(mt12_dedou_combined,"orig.ident",assay.use="SCT")
mt12_dedou_combined<-RunTSNE(mt12_dedou_combined,dims = 1:20)
mt12_dedou_combined<-RunUMAP(mt12_dedou_combined,reduction = "harmony",dims = 1:15)
mt12_dedou_combined<-FindNeighbors(mt12_dedou_combined,dims = 1:20)
mt12_dedou_combined<-FindClusters(mt12_dedou_combined,resolution = 0.5)

MK_combine<-subset(mt12_dedou_combined,idents=c("0","1","9","12","17"))
DefaultAssay(MK_combine)<-"SCT"
MK_combine<-FindVariableFeatures(MK_combine,selection.method = "vst",nfeatures = 2000)
MK_combine<-ScaleData(MK_combine,verbose = FALSE)
MK_combine<-RunPCA(MK_combine,npcs = 20,verbose = FALSE)
MK_combine<-SCTransform(MK_combine, 
                        vars.to.regress = "percent.mt", verbose = FALSE)
MK_combine<-RunHarmony(MK_combine,"orig.ident",assay.use="SCT")
MK_combine<-RunTSNE(MK_combine,dims = 1:20)
MK_combine<-RunUMAP(MK_combine,reduction = "harmony",dims = 1:5)
MK_combine<-FindNeighbors(MK_combine,dims = 1:20)
MK_combine<-FindClusters(MK_combine,resolution = 0.08)

ct_MK<-row.names(subset(MK_combine@meta.data,
                        MK_combine@meta.data$orig.ident=="Ctrl_all"))
lm_MK<-row.names(subset(MK_combine@meta.data,
                        MK_combine@meta.data$orig.ident=="lm_all"))
ct_MK_seu<-subset(MK_combine,cells=ct_MK)
lm_MK_seu<-subset(MK_combine,cells=lm_MK)

ct_MK_seu<-FindVariableFeatures(ct_MK_seu,selection.method = "vst",nfeatures = 2000)
ct_MK_seu<-ScaleData(ct_MK_seu,verbose = FALSE)
ct_MK_seu<-RunPCA(ct_MK_seu,npcs = 20,verbose = FALSE)
ct_MK_seu<-SCTransform(ct_MK_seu, 
                        vars.to.regress = "percent.mt", verbose = FALSE)
ct_MK_seu<-RunHarmony(ct_MK_seu,"orig.ident",assay.use="SCT")
ct_MK_seu<-RunTSNE(ct_MK_seu,dims = 1:20)
ct_MK_seu<-RunUMAP(ct_MK_seu,reduction = "harmony",dims = 1:13)
ct_MK_seu<-FindNeighbors(ct_MK_seu,dims = 1:20)
ct_MK_seu<-FindClusters(ct_MK_seu,resolution = 0.55)
ct_markers<-FindAllMarkers(object = ct_MK_seu,only.pos = TRUE, 
                           min.pct = 0.1, logfc.threshold = 0.1)

ct_df<-as.data.frame(ct_MK_seu@assays[["RNA"]]@data)
ct_df_sp<-ct_df[c("Cxcr4","Tnf","Mki67","Mcm2","Ccl3","Il6","Cd74",
                  "Ccr2","Cdh1","S100a9","Gp9","Myl9","F13a1",
                  "Fcgr3","Pf4","Itga2b","Spi1","Il1b"),]
ct_df_sp_t<-as.data.frame(t(ct_df_sp))
data4pea_scale<-as.data.frame(scale(ct_df_sp_t))
peas_1<-cor(data4pea_scale, method = "spearman");peas_1
pas<-cor.test(data4pea_scale$Cxcr4,data4pea_scale$Il1b, method = "spearman");pas