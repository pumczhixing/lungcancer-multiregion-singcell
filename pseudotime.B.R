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
load("../Bcells_seurat_cluster.RData")
info <- read.table("/glusterfs/home/local_xia_hao/Guo_wei/scRNA/Guo_wei/New_sample_20200811/scRNA/dataDelivery/6_malignant/Pseudotime_epithelial_to_tumor_cells/sample.info",sep="\t",header=T,stringsAsFactors=F)
rownames(info) <- info$Sample

temp_matrix <- Bcells_obj@meta.data[,c("sample","index")]
for( i in rownames(info)){temp_matrix[temp_matrix[,"sample"]==i,"index"]=info[i,"TMB"]}
Bcells_obj=AddMetaData(object = Bcells_obj, metadata = temp_matrix[,"index"], col.name = "TMB")
for( i in rownames(info)){temp_matrix[temp_matrix[,"sample"]==i,"index"]=info[i,"TMB_H"]}
Bcells_obj=AddMetaData(object = Bcells_obj, metadata = temp_matrix[,"index"], col.name = "TMB_H")
for( i in rownames(info)){temp_matrix[temp_matrix[,"sample"]==i,"index"]=info[i,"MSI"]}
Bcells_obj=AddMetaData(object = Bcells_obj, metadata = temp_matrix[,"index"], col.name = "MSI")
for( i in rownames(info)){temp_matrix[temp_matrix[,"sample"]==i,"index"]=info[i,"MSS_e"]}
Bcells_obj=AddMetaData(object = Bcells_obj, metadata = temp_matrix[,"index"], col.name = "MSS_e")

object.seurat <- Bcells_obj

save(file="Bcells_obj.rData",Bcells_obj)

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
ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]


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

pdf(paste(projectname,"_pseudotime_color_by_","celltype.3.",type,".pdf",sep = ""),width = 8,height = 6)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "Cell_type",cell_size = 0.2,cell_link_size = 0.5))
dev.off()

pdf(paste(projectname,"_pseudotime_color_by_","celltype_each.3.",type,".pdf",sep = ""),width = 12,height = 6)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "Cell_type",cell_size = 0.2,cell_link_size = 0.5) + facet_wrap(~Cell_type, nrow = 1))
dev.off()

pdf(paste(projectname,"_pseudotime_color_by_","Infiltration.3.",type,".pdf",sep = ""),width =5,height = 5)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "typeCNV",cell_size = 0.2,cell_link_size = 0.5))
dev.off()

pdf(paste(projectname,"_pseudotime_color_by_","typeCNV.3.",type,".pdf",sep = ""),width =5,height = 5)
print(plot_cell_trajectory(clustered_spleen_monocle, color_by = "typeCNV",cell_size = 0.2,cell_link_size = 0.5)+scale_colour_manual(values = c("#00BA38","#F8766D","#619CFF")))
dev.off()





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

pseudoresult(object.seurat,"Bcells_cells","All")


Bcells_obj <- RenameIdents(Bcells_obj, `0` = "0", `1` = "0" , `5` = "0", `6` = "0", `7`="0", `3`="3", `8` = "1", `2` = "2",`3` = "3", `4` = "4")
Bcells_obj <- AddMetaData(object = Bcells_obj, metadata = Idents(Bcells_obj), col.name = "Cluster")
for( i in rownames(pData(clustered_spleen_monocle))){pData(clustered_spleen_monocle)[i,"Final_Cluster"]=Bcells_obj@meta.data[i,"Cluster"]}

markers_C3_vs_C0 <- FindMarkers(Bcells_obj, ident.1 = "3" ,ident.2 = "0")
markers_C2_vs_C0 <- FindMarkers(Bcells_obj, ident.1 = "2" ,ident.2 = "0")
markers_C3_vs_C2 <- FindMarkers(Bcells_obj, ident.1 = "3" ,ident.2 = "2")
write.table(markers_C3_vs_C0,"Bcells_markers_C3_vs_C0_markers.txt",sep="\t",quote=F)
write.table(markers_C2_vs_C0,"Bcells_markers_C2_vs_C0_markers.txt",sep="\t",quote=F)
write.table(markers_C3_vs_C2,"Bcells_markers_C3_vs_C2_markers.txt",sep="\t",quote=F)

genelist <- read.table("genelist.xls",header=F)
genelist <- genelist[,1]
clustered_spleen_monocle_gene <- row.names(subset(fData(clustered_spleen_monocle),gene_short_name %in% genelist))
pdf(paste("Bcell","_pseudotime_","Pseudotime",".pdf",sep = ""),width =6,height = 10)
plot_genes_branched_pseudotime(clustered_spleen_monocle[clustered_spleen_monocle_gene,],branch_point = 1,color_by = "Infiltration",ncol = 1)
dev.off()





genelist <- read.table("genelist.1.xls",header=F)
genelist <- genelist[,1]
pdf(paste("Bcell","_pseudotime_","Pseudotime",".pdf",sep = ""),width =6,height = 30)
clustered_spleen_monocle_gene <- row.names(subset(fData(clustered_spleen_monocle),gene_short_name %in% genelist))
plot_genes_branched_pseudotime(clustered_spleen_monocle[clustered_spleen_monocle_gene,],branch_point = 1,color_by = "Infiltration",ncol = 1)
genelist <- read.table("genelist.2.xls",header=F)
genelist <- genelist[,1]
clustered_spleen_monocle_gene <- row.names(subset(fData(clustered_spleen_monocle),gene_short_name %in% genelist))
plot_genes_branched_pseudotime(clustered_spleen_monocle[clustered_spleen_monocle_gene,],branch_point = 1,color_by = "Infiltration",ncol = 1)
genelist <- read.table("genelist.3.xls",header=F)
genelist <- genelist[,1]
clustered_spleen_monocle_gene <- row.names(subset(fData(clustered_spleen_monocle),gene_short_name %in% genelist))
plot_genes_branched_pseudotime(clustered_spleen_monocle[clustered_spleen_monocle_gene,],branch_point = 1,color_by = "Infiltration",ncol = 1)
genelist <- read.table("genelist.4.xls",header=F)
genelist <- genelist[,1]
clustered_spleen_monocle_gene <- row.names(subset(fData(clustered_spleen_monocle),gene_short_name %in% genelist))
plot_genes_branched_pseudotime(clustered_spleen_monocle[clustered_spleen_monocle_gene,],branch_point = 1,color_by = "Infiltration",ncol = 1)
genelist <- read.table("genelist.5.xls",header=F)
genelist <- genelist[,1]
clustered_spleen_monocle_gene <- row.names(subset(fData(clustered_spleen_monocle),gene_short_name %in% genelist))
plot_genes_branched_pseudotime(clustered_spleen_monocle[clustered_spleen_monocle_gene,],branch_point = 1,color_by = "Infiltration",ncol = 1)
dev.off()


genelist <- read.table("genelist.1.xls",header=F)
genelist <- genelist[,1]
pdf(paste("Bcell","_pseudotime_","Pseudotime_Final_Cluster",".pdf",sep = ""),width =6,height = 25)
clustered_spleen_monocle_gene <- row.names(subset(fData(clustered_spleen_monocle),gene_short_name %in% genelist))
plot_genes_branched_pseudotime(clustered_spleen_monocle[clustered_spleen_monocle_gene,],branch_point = 1,color_by = "Final_Cluster",ncol = 1)
genelist <- read.table("genelist.2.xls",header=F)
genelist <- genelist[,1]
clustered_spleen_monocle_gene <- row.names(subset(fData(clustered_spleen_monocle),gene_short_name %in% genelist))
plot_genes_branched_pseudotime(clustered_spleen_monocle[clustered_spleen_monocle_gene,],branch_point = 1,color_by = "Final_Cluster",ncol = 1)
genelist <- read.table("genelist.3.xls",header=F)
genelist <- genelist[,1]
clustered_spleen_monocle_gene <- row.names(subset(fData(clustered_spleen_monocle),gene_short_name %in% genelist))
plot_genes_branched_pseudotime(clustered_spleen_monocle[clustered_spleen_monocle_gene,],branch_point = 1,color_by = "Final_Cluster",ncol = 1)
genelist <- read.table("genelist.4.xls",header=F)
genelist <- genelist[,1]
clustered_spleen_monocle_gene <- row.names(subset(fData(clustered_spleen_monocle),gene_short_name %in% genelist))
plot_genes_branched_pseudotime(clustered_spleen_monocle[clustered_spleen_monocle_gene,],branch_point = 1,color_by = "Final_Cluster",ncol = 1)
genelist <- read.table("genelist.5.xls",header=F)
genelist <- genelist[,1]
clustered_spleen_monocle_gene <- row.names(subset(fData(clustered_spleen_monocle),gene_short_name %in% genelist))
plot_genes_branched_pseudotime(clustered_spleen_monocle[clustered_spleen_monocle_gene,],branch_point = 1,color_by = "Final_Cluster",ncol = 1)
dev.off()



