library(DOSE)
library(GOSemSim)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(GO.db)
library(tibble)
library(enrichplot)
library(ggplot2)

#function
get_GO_data <- function(OrgDb, ont, keytype) {
  GO_Env <- get_GO_Env()
  use_cached <- FALSE
  
  if (exists("organism", envir=GO_Env, inherits=FALSE) &&
      exists("keytype", envir=GO_Env, inherits=FALSE)) {
    
    org <- get("organism", envir=GO_Env)
    kt <- get("keytype", envir=GO_Env)
    
    if (org == DOSE:::get_organism(OrgDb) &&
        keytype == kt &&
        exists("goAnno", envir=GO_Env, inherits=FALSE)) {
      use_cached <- TRUE
    }
  }
  
  if (use_cached) {
    goAnno <- get("goAnno", envir=GO_Env)
  } else {
    OrgDb <- GOSemSim:::load_OrgDb(OrgDb)
    kt <- keytypes(OrgDb)
    if (! keytype %in% kt) {
      stop("keytype is not supported...")
    }
    
    kk <- keys(OrgDb, keytype=keytype)
    goAnno <- suppressMessages(
      AnnotationDbi::select(OrgDb, keys=kk, keytype=keytype,
             columns=c("GOALL", "ONTOLOGYALL")))
    
    goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])
    
    assign("goAnno", goAnno, envir=GO_Env)
    assign("keytype", keytype, envir=GO_Env)
    assign("organism", DOSE:::get_organism(OrgDb), envir=GO_Env)
  }
  
  if (ont == "ALL") {
    GO2GENE <- unique(goAnno[, c(2,1)])
  } else {
    GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
  }
  
  GO_DATA <- DOSE:::build_Anno(GO2GENE, get_GO2TERM_table())
  
  goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
  goOnt <- goOnt.df[,2]
  names(goOnt) <- goOnt.df[,1]
  assign("GO2ONT", goOnt, envir=GO_DATA)
  return(GO_DATA)
}

get_GO_Env <- function () {
  if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
  }
  get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}

get_GO2TERM_table <- function() {
  GOTERM.df <- get_GOTERM()
  GOTERM.df[, c("go_id", "Term")] %>% unique
}

get_GOTERM <- function() {
  pos <- 1
  envir <- as.environment(pos)
  if (!exists(".GOTERM_Env", envir=envir)) {
    assign(".GOTERM_Env", new.env(), envir)
  }
  GOTERM_Env <- get(".GOTERM_Env", envir = envir)
  if (exists("GOTERM.df", envir = GOTERM_Env)) {
    GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
  } else {
    GOTERM.df <- toTable(GOTERM)
    assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
  }
  return(GOTERM.df)
}
#get_GO_data_GSEA
GO_DATA <- get_GO_data("org.Mm.eg.db", "ALL", "SYMBOL") 
GO_Gene_1 <- unique(unlist(GO_DATA$PATHID2EXTID["GO:0006954"]));GO_Gene_1

gene_go1<-cbind(rep("xxx",length(GO_Gene_1)),as.data.frame(GO_Gene_1))
colnames(gene_go1)<-c("ont","gene")
gene_fc_df<-read.csv("/media/user/sde/MK.csv")
gene_df<-gene_fc_df
Gene_id<-bitr(gene_df$SYMBOL, fromType="SYMBOL", toType=c("SYMBOL","ENTREZID"),
              OrgDb="org.Mm.eg.db")

gesa_gene_go<-GSEA(gene_df_list, nPerm = 10000, TERM2GENE = Hallmark_geneset,
                   verbose=FALSE, pvalueCutoff = 5)

gsea_result<-data.frame(gesa_gene_go)
gseaplot<-gseaplot2(gesa_gene_go,3,pvalue_table = T,color = "darkgrey",title = "GO:0045669");gseaplot