library(futile.logger)
flog.info("::InferCNV program:Start")
flog.info("::Load function:Start")
suppressMessages(library('Seurat'))
suppressMessages(library(ggplot2))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
library(Seurat)
library(monocle)
library(reshape2)

#source("/glusterfs/home/xia_hao0101/project_202007/scRNA/Guo_wei/New_sample_20200811/scRNA/src/inferCNV.R")
source("/glusterfs/home/local_xia_hao/project_202007/scRNA/Guo_wei/New_sample_20200811/scRNA/src/inferCNV.R")
flog.info("::Load input file:Start")
#gene=read.delim("/titan2/liu_yuhao0101/database/inferCNV/gene.GRCH38.all.gene.txt",row.names=1)
#object.seurat=eval(parse(text = load("../all_cells_include_typeCNV.rData")))
#obj_epi <-subset(obj,idents=c("Alveolar cells","Epithelial cells"))
#load("../Endothelialcells_seurat_cluster.RData")
load("../Myeloid_cells_seurat_cluster.RData")

a <- list()
temp <- table(Idents(Myeloid_cells_obj))
for( i in unique(Idents(Myeloid_cells_obj))){
    samples_cell <- Idents(Myeloid_cells_obj)[ Idents(Myeloid_cells_obj)==i ]
    a[[paste("cell_name",i,sep="_")]] <- names(sample(samples_cell, temp[i] %/% 3,replace=F))
}

choose_cells<-list()
for( i in unique(Idents(Myeloid_cells_obj))){
    choose_cells <- unlist(list(choose_cells,c(a[[paste("cell_name",i,sep="_")]])))
}

#Tcell_obj_tmp <-  subset(Tcell_obj,cells = c(a$cell_name_0,a$cell_name_1,a$cell_name_2,a$cell_name_3,
#                   a$cell_name_4,a$cell_name_5,a$cell_name_6,a$cell_name_7,a$cell_name_8))

Myeloid_cells_obj_tmp <-  subset(Myeloid_cells_obj,cells=choose_cells)

info <- read.table("/glusterfs/home/local_xia_hao/Guo_wei/scRNA/Guo_wei/New_sample_20200811/scRNA/dataDelivery/6_malignant/Pseudotime_epithelial_to_tumor_cells/sample.info",sep="\t",header=T,stringsAsFactors=F)
rownames(info) <- info$Sample
cells_obj <- Myeloid_cells_obj_tmp
temp_matrix <- cells_obj@meta.data[,c("sample","index")]
for( i in rownames(info)){temp_matrix[temp_matrix[,"sample"]==i,"index"]=info[i,"TMB"]}
cells_obj=AddMetaData(object = cells_obj, metadata = temp_matrix[,"index"], col.name = "TMB")
for( i in rownames(info)){temp_matrix[temp_matrix[,"sample"]==i,"index"]=info[i,"TMB_H"]}
cells_obj=AddMetaData(object = cells_obj, metadata = temp_matrix[,"index"], col.name = "TMB_H")
for( i in rownames(info)){temp_matrix[temp_matrix[,"sample"]==i,"index"]=info[i,"MSI"]}
cells_obj=AddMetaData(object = cells_obj, metadata = temp_matrix[,"index"], col.name = "MSI")
for( i in rownames(info)){temp_matrix[temp_matrix[,"sample"]==i,"index"]=info[i,"MSS_e"]}
cells_obj=AddMetaData(object = cells_obj, metadata = temp_matrix[,"index"], col.name = "MSS_e")

Myeloid_cells_obj_tmp <- cells_obj
object.seurat <- cells_obj

save(file="Myeloid_cells_obj.rData",Myeloid_cells_obj_tmp)

pseudoresult <- function(obj_tmp_test,type,projectname){

data <- as(as.matrix(obj_tmp_test@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = obj_tmp_test@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
head(pData(monocle_cds))

clustered_spleen_monocle  <- monocle_cds



names(pData(clustered_spleen_monocle))[names(pData(clustered_spleen_monocle))=="RNA_snn_res.0.5"]="Cluster"
clustered_spleen_monocle <- estimateSizeFactors(clustered_spleen_monocle)
clustered_spleen_monocle <- estimateDispersions(clustered_spleen_monocle)

diff_test_res <- differentialGeneTest(clustered_spleen_monocle,fullModelFormulaStr = "~Cluster",cores=20)

save(diff_test_res, file = paste("diff_test_res.3",type,"Rdata",sep="."))
#ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
#  ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
ordering_genes <- row.names(diff_test_res[diff_test_res$qval<0.01,])[ order((diff_test_res[diff_test_res$qval<0.01,])$qval)][1:2000]

clustered_spleen_monocle <- setOrderingFilter(clustered_spleen_monocle, ordering_genes)
save(clustered_spleen_monocle, file = paste("clustered_spleen_monocle.3",type,"Rdata",sep="."))

pdf(paste(projectname,"_plot_ordering_genes.",type,".3.pdf",sep = ""),width = 8,height = 6)
print(plot_ordering_genes(clustered_spleen_monocle))
dev.off()

clustered_spleen_monocle <- reduceDimension(clustered_spleen_monocle, max_components = 2, reduction_method = "DDRTree",norm_method = 'log',cores=20)
clustered_spleen_monocle <- orderCells(clustered_spleen_monocle)
save(clustered_spleen_monocle, file = paste("clustered_spleen_monocle.2.3",type,"Rdata",sep="."))


pdf(paste(projectname,"_pseudotime_color_by_","cluster.3.",type,".pdf",sep = ""),width = 8,height = 6)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "Cluster",cell_size = 0.2,cell_link_size = 0.5))
dev.off()

pdf(paste(projectname,"_pseudotime_color_by_","cluster_each.3.",type,".pdf",sep = ""),width = 12,height = 6)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "Cluster",cell_size = 0.2,cell_link_size = 0.5) + facet_wrap(~Cluster, nrow = 1))
dev.off()

pdf(paste(projectname,"_pseudotime_color_by_","pseudotime.3.",type,".pdf",sep = ""),width = 8,height = 6)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "Pseudotime",cell_size = 0.2,cell_link_size = 0.5))
dev.off()

pdf(paste(projectname,"_pseudotime_color_by_","Infiltration.3.",type,".pdf",sep = ""),width = 8,height = 6)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "Infiltration",cell_size = 0.2,cell_link_size = 0.5))
dev.off()


pdf(paste(projectname,"_pseudotime_color_by_","Infiltration_each.3.",type,".pdf",sep = ""),width = 12,height = 6)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "Infiltration",cell_size = 0.2,cell_link_size = 0.5)+ facet_wrap(~Infiltration, nrow = 1))
dev.off()

pdf(paste(projectname,"_pseudotime_color_by_","clust_info.3.",type,".pdf",sep = ""),width = 8,height = 6)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "clust_info",cell_size = 0.2,cell_link_size = 0.5))
dev.off()

pdf(paste(projectname,"_pseudotime_color_by_","clust_info_each.3.",type,".pdf",sep = ""),width = 12,height = 6)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "clust_info",cell_size = 0.2,cell_link_size = 0.5)+ facet_wrap(~clust_info, nrow = 1))
dev.off()


#pdf(paste(projectname,"_pseudotime_color_by_","celltype.3.",type,".pdf",sep = ""),width = 8,height = 6)
#print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "Cell_type",cell_size = 0.2,cell_link_size = 0.5))
#dev.off()

#pdf(paste(projectname,"_pseudotime_color_by_","celltype_each.3.",type,".pdf",sep = ""),width = 12,height = 6)
#print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "Cell_type",cell_size = 0.2,cell_link_size = 0.5) + facet_wrap(~Cell_type, nrow = 1))
#dev.off()

#pdf(paste(projectname,"_pseudotime_color_by_","Infiltration.3.",type,".pdf",sep = ""),width =5,height = 5)
#print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "typeCNV",cell_size = 0.2,cell_link_size = 0.5))
#dev.off()

#pdf(paste(projectname,"_pseudotime_color_by_","typeCNV.3.",type,".pdf",sep = ""),width =5,height = 5)
#print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "typeCNV",cell_size = 0.2,cell_link_size = 0.5)+scale_colour_manual(values = c("#00BA38","#F8766D","#619CFF")))
#dev.off()





prefix <- projectname

for (i in c(1)){

branch_point=i
BEAM_res <- BEAM(clustered_spleen_monocle, branch_point = branch_point, cores = 20)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

pdf(paste(prefix,branch_point,type,"pseudotime_heatmap.3.2.pdf",sep='_'), w=8, h=50)
print(plot_genes_branched_heatmap(clustered_spleen_monocle[row.names(subset(BEAM_res,qval < 1e-5)),],branch_point = branch_point,num_clusters = 3,cores = 4,use_gene_short_name = T,show_rownames = T))
dev.off()
}
}

pseudoresult(object.seurat,"Myeloid_cells_obj","All")


exprData <- clustered_spleen_monocle@assayData$exprs
exprData <- LogNormalize(exprData)
genelist <- c("FDFT1","CD52","AQP3","MCEMP1","S100A13","FABP3","GCHFR","TGM2","SCD","FABP4","RBP4","SERPING1","RND3","CES1","GLDN","INHBA","AC026369.3","TSPAN3","PLA2G16","SVIL",
"OASL","ACOT7","HPGD","NMB","IGFBP2","CSF1","PHLDA3","PCOLCE2","ADTRP","PDLIM1","FHL1","MLPH","FAM89A","AC013457.1","MME","GPD1","CPE","LINC02154","SLC19A3")

for(i in genelist){
if(i %in% rownames(BEAM_res)){
clustered_spleen_monocle$TEMP <- exprData[i,]
pdf(paste(prefix,branch_point,i,"pseudotime.pdf",sep='_'), w=8, h=6)
p<-plot_cell_trajectory(clustered_spleen_monocle, color_by = "TEMP" ,cell_size=1,show_backbone=TRUE)
print(p + scale_color_viridis_c())
dev.off()
}
}




genelist <- c("C020656.1","PMP22","ITGAM","CCL20","LYZ","BLVRB","S100A8","CFD","SDC2","NPL","CPM","GNPDA1","OTOA","SPP1","PLA2G7","LGMN","APOE","TREM2","HAMP","CHI3L1","CHIT1","CTSH","HDHD5","GSDME","DNPH1","ZFAND5","THBS1","C5AR1","TNFAIP6","CH25H","WIPI1","SLC1A3","CCL3","DMXL2","IL6","CD36","CXCL2","AQP9")



for( i in rownames(pData(clustered_spleen_monocle))){pData(clustered_spleen_monocle)[i,"Final_Cluster"]=Myeloid_cells_obj_tmp@meta.data[i,"clust_info"]}


pdf(paste(projectname,"_pseudotime_color_by_","Final_Cluster.",type,".pdf",sep = ""),width = 8,height = 6)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "Final_Cluster",cell_size = 0.2,cell_link_size = 0.5))
dev.off()

pdf(paste(projectname,"_pseudotime_color_by_","Final_Cluster_each.",type,".pdf",sep = ""),width = 12,height = 6)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "Final_Cluster",cell_size = 0.1,cell_link_size = 0.5)+ facet_wrap(~Final_Cluster, nrow = 1))
dev.off()

branch2gene <- read.table("branch2_gene.list",header=F)
branch2gene <- branch2gene[,1]
setwd("/glusterfs/home/local_xia_hao/Guo_wei/scRNA/Guo_wei/New_sample_20200811/scRNA/dataDelivery/Myeloid_cells/Pseudotime_Myeloid_cells_3/Branch2_Gene")
for(i in branch2gene){
if(i %in% rownames(BEAM_res)){
    assign(paste(i,"",sep=""),i)
    pData(clustered_spleen_monocle)[get(paste(i,"",sep=""))] <- exprData[i,]
    pdf(paste(prefix,branch_point,i,"pseudotime.pdf",sep='_'), w=8, h=6)
    p<-plot_cell_trajectory(clustered_spleen_monocle, color_by = get(paste(i,"",sep="")) ,cell_size=1,show_backbone=TRUE)
    print(p + scale_color_viridis_c())
    dev.off()
}
}


