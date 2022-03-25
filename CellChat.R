library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(future)
library(Seurat)

ct_all_seu_rn_nwid<-readRDS("/media/user/sdd/ctrl_all_MK.rds")
data.input<-ct_all_seu_rn_nwid@assays$RNA@data
###merge####
data.input1<-as.data.frame(ct_all_seu_rn_nwid@assays$RNA@data)
data.input2<-as.data.frame(MF_mf_seu@assays$RNA@data)
data.input<-merge(data.input1,data.input2,by="row.names",all=T)
data.input[is.na(data.input)==TRUE]<-0
row.names(data.input)<-data.input$Row.names
data.input<-data.input[,-1]
data.input<-as.matrix(data.input)

identity<-as.data.frame(ct_all_seu_rn_nwid@active.ident)
colnames(identity)<-"labels"
###merge####
identity1<-as.data.frame(ct_all_seu_rn_nwid@active.ident)
colnames(identity1)<-"labels"
identity2<-as.data.frame(MF_mf_seu@active.ident)
colnames(identity2)<-"labels"
identity<-rbind(identity1,identity2)

cellchat<-createCellChat(data.input)
cellchat<-addMeta(cellchat,meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
levels(cellchat@idents)
as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.mouse 
#unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
#ligand<-unique(CellChatDB.use[["interaction"]][["ligand"]])
#receptor<-unique(CellChatDB.use[["interaction"]][["receptor"]])
#CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
#CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 1)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat,PPI.mouse)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
df.net <- subsetCommunication(cellchat,slot.name = "net")
df.netP <- subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.net,"/media/user/sdd/cellchat_secreted.csv")

df.dot<-read.csv("/media/user/sdd/cellchat_secreted_select.csv")
for (i in 1:nrow(df.dot)){
  df.dot[i,8]<-paste0(df.dot[i,1],"_",df.dot[i,2])
}
colnames(df.dot)[8]<-"pair"

dotp<-ggplot(df.dot,aes(x=pair,y=interaction_name_2))+
  geom_point(aes(size=pval,colour=(log2(prob))))+
  scale_size(limits = c(1,2))+
  scale_color_gradientn(colors = c("#F43B3F","#F3E544","#3EBEC3","#1E3164")[4:1])
ggsave("/media/user/sdd/cellchat/cellchat_dot.pdf",
       dotp,width = 8,height = 4)