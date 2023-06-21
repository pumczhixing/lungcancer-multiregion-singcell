
library(ggplot2)
library(reshape2)
library(ggsci)
library(gridExtra)
library(Seurat)
library(monocle)
library(reshape)
library("dplyr")
library("Matrix")
library(grid)


load("Tcell_seurat_cluster.2.RData")

sample_info <- read.table("sample.4.info",sep="\t",stringsAsFactors=F)
matrix <- as.matrix(Tcell_obj@meta.data[,"sample"])
matrix_2 <- matrix
matrix_3 <- matrix
rownames(matrix) <- rownames(Tcell_obj@meta.data)
rownames(matrix_2) <- rownames(Tcell_obj@meta.data)
rownames(matrix_3) <- rownames(Tcell_obj@meta.data)

for( i in 1:23){matrix <- gsub(sample_info[i,"V3"],sample_info[i,"V2"],matrix)}
Tcell_obj <- AddMetaData(object = Tcell_obj, metadata = matrix, col.name = "Infiltration")

for( i in 1:23){matrix_2 <- gsub(sample_info[i,"V3"],sample_info[i,"V5"],matrix_2)}
Tcell_obj <- AddMetaData(object = Tcell_obj, metadata = matrix_2, col.name = "Patients_info")

for( i in 1:23){matrix_3 <- gsub(sample_info[i,"V3"],sample_info[i,"V4"],matrix_3)}
Tcell_obj <- AddMetaData(object = Tcell_obj, metadata = matrix_3, col.name = "Sample_info")

project="All_samples_Tcells"


#tcr_info="QJM/L191024-P1-TCR_filtered_contig_annotations.csv"

tcr_plot <- function(file,lable_tem){
    qjm_tcr_N=""
    qjm_tcr_N <- read.table(file,header=T,sep=",",stringsAsFactors=F)
    
    for( i in 1:nrow(qjm_tcr_N)){if(qjm_tcr_N[i,"raw_clonotype_id"]!="clonotype1" & qjm_tcr_N[i,"raw_clonotype_id"]!="clonotype2" & qjm_tcr_N[i,"raw_clonotype_id"]!="clonotype3" & qjm_tcr_N[i,"raw_clonotype_id"]!="clonotype4" & qjm_tcr_N[i,"raw_clonotype_id"]!="clonotype5"){qjm_tcr_N[i,"raw_clonotype_id"]="TCR"} }
Tcr_info<- (qjm_tcr_N[,c("barcode","raw_clonotype_id")])[!duplicated(qjm_tcr_N[,c("barcode","raw_clonotype_id")], fromLast=TRUE), ]
    change_name=paste(".", lable_tem , sep="")

    for( i in 1:nrow(Tcr_info)){Tcr_info[,"barcode"] <- gsub("-1",change_name ,Tcr_info[,"barcode"])}

    matrt <- as.matrix(Tcr_info$raw_clonotype_id)
    rownames(matrt) <- Tcr_info$barcode
    change_metadata_name= paste("Patients_",lable_tem,"_TCR_info",sep="")
    print(change_metadata_name)
    change_metadata_name=gsub("-","_",change_metadata_name)
    Tcell_obj_tmp <- AddMetaData(object = Tcell_obj, metadata = matrt, col.name = change_metadata_name)
#eemp <- qjm_tcr_N[,c("barcode","raw_clonotype_id")]
pdf(file=paste0(project,"_",lable_tem,"_tcr_Umap_cells.pdf"),width=8,height=6, onefile=FALSE)
        print(DimPlot(Tcell_obj_tmp, reduction = "umap",label=FALSE,group.by= change_metadata_name ,cols=c("#FF0000","#CCCC00","#99CC00","#33CCFF","#0000FF","#77DDFF"),pt.size=1,na.value="#DDDDDD" ,cols.highlight=rownames(matrt) ))
    dev.off()

}


tcr_all_plot <- function(file1,lable_tem){
    qjm_tcr_N=""
    qjm_tcr_N <- read.table(file1,header=T,sep=",",stringsAsFactors=F)

    for( i in 1:nrow(qjm_tcr_N)){if(qjm_tcr_N[i,"raw_clonotype_id"]!="clonotype1" & qjm_tcr_N[i,"raw_clonotype_id"]!="clonotype2" & qjm_tcr_N[i,"raw_clonotype_id"]!="clonotype3" & qjm_tcr_N[i,"raw_clonotype_id"]!="clonotype4" & qjm_tcr_N[i,"raw_clonotype_id"]!="clonotype5"){qjm_tcr_N[i,"raw_clonotype_id"]="TCR"} }
Tcr_info<- (qjm_tcr_N[,c("barcode","raw_clonotype_id")])[!duplicated(qjm_tcr_N[,c("barcode","raw_clonotype_id")], fromLast=TRUE), ]
    change_name=paste(".", lable_tem , sep="")

    for( i in 1:nrow(Tcr_info)){Tcr_info[,"barcode"] <- gsub("-1",change_name ,Tcr_info[,"barcode"])}

    matrt <- as.matrix(Tcr_info$raw_clonotype_id)
    rownames(matrt) <- Tcr_info$barcode
    
    return(matrt)
}

qjm_tcr_N <- tcr_all_plot("QJM/L191024-P1-TCR_filtered_contig_annotations.csv","L191024-P1")
lsd_tcr_N <- tcr_all_plot("LSD/L191111-P-TCR_filtered_contig_annotations.csv","L191111-P")
yxq_tcr_N <- tcr_all_plot("YXQ/L191112-N1-TCR_filtered_contig_annotations.csv","L191112-N1")
zy_tcr_N <- tcr_all_plot("ZY/L191126-P1-TCR_filtered_contig_annotations.csv","L191126-P1")
dcy_tcr_N <- tcr_all_plot("DCY/L191126-P2-TCR_filtered_contig_annotations.csv","L191126-P2")
ll_tcr_N <- tcr_all_plot("LL/L191209-P1-TCR_filtered_contig_annotations.csv","L191209-P1")


all_tcr_N <-do.call("rbind", list(qjm_tcr_N,lsd_tcr_N,yxq_tcr_N,zy_tcr_N,dcy_tcr_N,ll_tcr_N))

Tcell_obj_all_tcr_N <- AddMetaData(object = Tcell_obj, metadata =  all_tcr_N, col.name = "All_Normal_TCR")

pdf(file=paste0(project,"_","All_Normal_TCR","_tcr_Umap_cells.pdf"),width=8,height=6, onefile=FALSE)
    print(DimPlot(Tcell_obj_all_tcr_N, reduction = "umap",label=FALSE,group.by= "All_Normal_TCR",cols=c("#FF0000","#CCCC00","#99CC00","#33CCFF","#0000FF","#77DDFF"),pt.size=1,na.value="#DDDDDD" ,cols.highlight=rownames(matrt) ))
dev.off()


qjm_tcr_L1 <- tcr_all_plot("QJM/L191024-T1-TCR_filtered_contig_annotations.csv","L191024-T1")
qjm_tcr_L2 <- tcr_all_plot("QJM/L200604-T1-TCR_filtered_contig_annotations.csv","L200604-T1")
qjm_tcr_L3 <- tcr_all_plot("QJM/L200604-T2-TCR_filtered_contig_annotations.csv","L200604-T2")

lsd_tcr_L <- tcr_all_plot("LSD/L191111-T2-TCR_filtered_contig_annotations.csv","L191111-T2")

dcy_tcr_L1 <- tcr_all_plot("DCY/L191126-T3-TCR_filtered_contig_annotations.csv","L191126-T3")
dcy_tcr_L2 <- tcr_all_plot("DCY/L191126-T4-TCR_filtered_contig_annotations.csv","L191126-T4")

ll_tcr_L1 <- tcr_all_plot("LL/L191209-T1-TCR_filtered_contig_annotations.csv","L191209-T1")
ll_tcr_L2 <- tcr_all_plot("LL/L200611-T-TCR_filtered_contig_annotations.csv","L200611-T")

all_tcr_L <-do.call("rbind", list(qjm_tcr_L1,qjm_tcr_L2,qjm_tcr_L3,lsd_tcr_L,dcy_tcr_L1,dcy_tcr_L2,ll_tcr_L1,ll_tcr_L2))
Tcell_obj_all_tcr_L <- AddMetaData(object = Tcell_obj, metadata =  all_tcr_L, col.name = "All_Low_infiltration_TCR")

pdf(file=paste0(project,"_","All_Low_infiltration_TCR","_tcr_Umap_cells.pdf"),width=8,height=6, onefile=FALSE)
    print(DimPlot(Tcell_obj_all_tcr_L, reduction = "umap",label=FALSE,group.by= "All_Low_infiltration_TCR",cols=c("#FF0000","#CCCC00","#99CC00","#33CCFF","#0000FF","#77DDFF"),pt.size=1,na.value="#DDDDDD" ,cols.highlight=rownames(matrt) ))
dev.off()




#qjm_tcr_N <- read.table("QJM/L191024-P1-TCR_filtered_contig_annotations.csv",header=T,sep=",",stringsAsFactors=F)
#lsd_tcr_N <- read.table("LSD/L191111-P-TCR_filtered_contig_annotations.csv",header=T,sep=",",stringsAsFactors=F)
#yxq_tcr_N <- read.table("YXQ/L191112-N1-TCR_filtered_contig_annotations.csv",header=T,sep=",",stringsAsFactors=F)
#zy_tcr_N <- read.table("ZY/L191126-P1-TCR_filtered_contig_annotations.csv",header=T,sep=",",stringsAsFactors=F)
#dcy_tcr_N <- read.table("DCY/L191126-P2-TCR_filtered_contig_annotations.csv",header=T,sep=",",stringsAsFactors=F)
#ll_tcr_N <- read.table("LL/L191209-P1-TCR_filtered_contig_annotations.csv",header=T,sep=",",stringsAsFactors=F)
lsd_tcr_D1 <- tcr_all_plot("LSD/L191111-T3-TCR_filtered_contig_annotations.csv","L191111-T3")
lsd_tcr_D2 <- tcr_all_plot("LSD/T-0817-VDJ-TCR_filtered_contig_annotations.csv","T-0817-VDJ")
yxq_tcr_D1 <- tcr_all_plot("YXQ/L191112-T1-TCR_filtered_contig_annotations.csv","L191112-T1")
yxq_tcr_D2 <- tcr_all_plot("YXQ/L191112-T2-TCR_filtered_contig_annotations.csv","L191112-T2")
zy_tcr_D1 <- tcr_all_plot("ZY/L191126-T1-TCR_filtered_contig_annotations.csv","L191126-T1")
zy_tcr_D2 <- tcr_all_plot("ZY/L191126-T2-TCR_filtered_contig_annotations.csv","L191126-T2")
zy_tcr_D3 <- tcr_all_plot("ZY/L200525-T-TCR_filtered_contig_annotations.csv","L200525-T")
dcy_tcr_D1 <- tcr_all_plot("DCY/L200519-T-TCR_filtered_contig_annotations.csv","L200519-T")
ll_tcr_D1 <- tcr_all_plot("LL/L191209-T2-TCR_filtered_contig_annotations.csv","L191209-T2")
all_tcr_D <-do.call("rbind", list(lsd_tcr_D1,lsd_tcr_D2,yxq_tcr_D1,yxq_tcr_D2,zy_tcr_D1,zy_tcr_D2,zy_tcr_D3,dcy_tcr_D1,ll_tcr_D1))
#yxq_tcr_D1 <- tcr_all_plot("yxq/L191112-T1-TCR_filtered_contig_annotations.csv","L191112-T1")

Tcell_obj_all_tcr_D <- AddMetaData(object = Tcell_obj, metadata =  all_tcr_D, col.name = "All_Deep_infiltration_TCR")

pdf(file=paste0(project,"_","All_Deep_infiltration_TCR","_tcr_Umap_cells.pdf"),width=8,height=6, onefile=FALSE)
    print(DimPlot(Tcell_obj_all_tcr_D, reduction = "umap",label=FALSE,group.by= "All_Deep_infiltration_TCR",cols=c("#FF0000","#CCCC00","#99CC00","#33CCFF","#0000FF","#77DDFF"),pt.size=1,na.value="#DDDDDD" ,cols.highlight=rownames(matrt) ))
dev.off()










tcr_plot("QJM/L191024-P1-TCR_filtered_contig_annotations.csv","L191024-P1")
tcr_plot("QJM/L191024-T1-TCR_filtered_contig_annotations.csv","L191024-T1")
tcr_plot("QJM/L200604-T1-TCR_filtered_contig_annotations.csv","L200604-T1")
tcr_plot("QJM/L200604-T2-TCR_filtered_contig_annotations.csv","L200604-T2")
tcr_plot("LSD/L191111-P-TCR_filtered_contig_annotations.csv","L191111-P")
tcr_plot("LSD/L191111-T2-TCR_filtered_contig_annotations.csv","L191111-T2")
tcr_plot("LSD/L191111-T3-TCR_filtered_contig_annotations.csv","L191111-T3")
tcr_plot("LSD/T-0817-VDJ-TCR_filtered_contig_annotations.csv","T-0817-VDJ")
tcr_plot("YXQ/L191112-T1-TCR_filtered_contig_annotations.csv","L191112-T1")
tcr_plot("YXQ/L191112-T2-TCR_filtered_contig_annotations.csv","L191112-T2")
tcr_plot("YXQ/L191112-N1-TCR_filtered_contig_annotations.csv","L191112-N1")
tcr_plot("ZY/L191126-P1-TCR_filtered_contig_annotations.csv","L191126-P1")
tcr_plot("ZY/L191126-T1-TCR_filtered_contig_annotations.csv","L191126-T1")
tcr_plot("ZY/L191126-T2-TCR_filtered_contig_annotations.csv","L191126-T2")
tcr_plot("ZY/L200525-T-TCR_filtered_contig_annotations.csv","L200525-T")
tcr_plot("DCY/L191126-P2-TCR_filtered_contig_annotations.csv","L191126-P2")
tcr_plot("DCY/L191126-T3-TCR_filtered_contig_annotations.csv","L191126-T3")
tcr_plot("DCY/L191126-T4-TCR_filtered_contig_annotations.csv","L191126-T4")
tcr_plot("DCY/L200519-T-TCR_filtered_contig_annotations.csv","L200519-T")
tcr_plot("LL/L191209-P1-TCR_filtered_contig_annotations.csv","L191209-P1")
tcr_plot("LL/L191209-T1-TCR_filtered_contig_annotations.csv","L191209-T1")
tcr_plot("LL/L191209-T2-TCR_filtered_contig_annotations.csv","L191209-T2")
tcr_plot("LL/L200611-T-TCR_filtered_contig_annotations.csv","L200611-T")






#pdf(file=paste0(project,"_qjm_tcr_N_Umap_cells.pdf"),width=8,height=6, onefile=FALSE)
#DimPlot(Tcell_obj, reduction = "umap",label=TRUE,)
#dev.off()

