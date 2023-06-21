library(Seurat)
library(monocle)
library(reshape)
library(celldex)
library(SingleR)
library(Matrix)
library(ggpubr)


projectname="T_cell"
project <- projectname

load("../../../../Tcell_seurat_cluster.2.RData")
temp <- table(Idents(Tcell_obj))

a <- list()

for( i in unique(Idents(Tcell_obj))){
    samples_Tcell <- Idents(Tcell_obj)[ Idents(Tcell_obj)==i ]
    a[[paste("cell_name",i,sep="_")]] <- names(sample(samples_Tcell, temp[i] %/% 10,replace=F))
}

Tcell_obj_tmp <-  subset(Tcell_obj,cells = c(a$cell_name_0,a$cell_name_1,a$cell_name_2,a$cell_name_3,
                   a$cell_name_4,a$cell_name_5,a$cell_name_6,a$cell_name_7,a$cell_name_8))

#细胞数减倍后的CD4
cd4_tcell_obj <- subset(Tcell_obj_tmp,idents=c('0','1','6','7'))
cd4_tcell_obj <- RenameIdents(cd4_tcell_obj, `0` = "Trm cells", `1` = "Naïve T cells", `2` = "NK cells/Exhausted T cells",`3`="Effector T/Trm cells",`4`="Activated/exhaustion T cells",`5`="Memory/effector T cells",`6`="Treg cells",`7`="Naïve /Memory T cells",`8` = "Cytotoxic T cells /NK cells")
cd4_tcell_obj  <- AddMetaData(cd4_tcell_obj, Idents(cd4_tcell_obj), col.name = "Cell_type")


pdf(file=paste0(project,"_Umap_cd4_T_start_cells.type.3.pdf"),width=12,height=7, onefile=FALSE)
DimPlot(cd4_tcell_obj, reduction = "umap",label = TRUE,cols =c("#3B499299","#EE000099","#008B4599","#63187999","#00828099","#BB002199","#5F559B99","#A2005699", "#80818099"))
dev.off()

cd4_tcell_obj <- NormalizeData(cd4_tcell_obj, normalization.method = "LogNormalize", scale.factor = 10000)
cd4_tcell_obj <- FindVariableFeatures(cd4_tcell_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(cd4_tcell_obj), 10)
all.genes <- rownames(cd4_tcell_obj)
cd4_tcell_obj <- ScaleData(cd4_tcell_obj, features = all.genes)
cd4_tcell_obj <- RunPCA(cd4_tcell_obj, features = VariableFeatures(object = cd4_tcell_obj))
ElbowPlot(cd4_tcell_obj)

project="cd4_tcell"
cd4_tcell_obj <- FindNeighbors(cd4_tcell_obj, dims = 1:20)
cd4_tcell_obj <- FindClusters(cd4_tcell_obj, resolution = 0.5)
cd4_tcell_obj <- RunUMAP(cd4_tcell_obj, dims = 1:20)

pdf(file=paste0(project,"_Umap_cd4_T_cells.type.pdf"),width=12,height=7, onefile=FALSE)
DimPlot(cd4_tcell_obj, reduction = "umap",label = TRUE)
dev.off()


data <- as(as.matrix(cd4_tcell_obj@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = cd4_tcell_obj@meta.data)
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
  
save(diff_test_res, file = paste("diff_test_res.3","CD4","Rdata",sep="."))
#ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
  
  
clustered_spleen_monocle <- setOrderingFilter(clustered_spleen_monocle, ordering_genes)
save(clustered_spleen_monocle, file = paste("clustered_spleen_monocle.3",type,"Rdata",sep="."))


###自此开始做拟时序分析

load("clustered_spleen_monocle.3.CD4.Rdata")#需要加载夏老师选来保存的数据，保证画出的图一致
pdf(paste(projectname,"_plot_ordering_genes.","CD4",".3.pdf",sep = ""),width = 8,height = 6)
plot_ordering_genes(clustered_spleen_monocle)
dev.off()
  
clustered_spleen_monocle <- reduceDimension(clustered_spleen_monocle, max_components = 2, reduction_method = "DDRTree",norm_method = 'log',cores=20)
clustered_spleen_monocle <- orderCells(clustered_spleen_monocle)
save(clustered_spleen_monocle, file = paste("clustered_spleen_monocle.2.3","CD4","Rdata",sep="."))

  
pdf(paste(projectname,"_pseudotime_color_by_","cluster.3.","CD4",".pdf",sep = ""),width = 8,height = 6)
plot_cell_trajectory(clustered_spleen_monocle, color_by = "Cluster",cell_size = 0.2,cell_link_size = 0.5)
dev.off()
  
pdf(paste(projectname,"_pseudotime_color_by_","cluster_each.3.","CD4",".pdf",sep = ""),width = 12,height = 6)
plot_cell_trajectory(clustered_spleen_monocle, color_by = "Cluster",cell_size = 0.2,cell_link_size = 0.5) + facet_wrap(~Cluster, nrow = 1)
dev.off()

pdf(paste(projectname,"_pseudotime_color_by_","pseudotime.3.","CD4",".pdf",sep = ""),width = 8,height = 6)
plot_cell_trajectory(clustered_spleen_monocle, color_by = "Pseudotime",cell_size = 0.2,cell_link_size = 0.5)
dev.off()
  
cols =c("#3B499299","#EE000099","#5F559B99","#A2005699")

pdf(paste(projectname,"_pseudotime_color_by_","celltype.3.","CD4",".pdf",sep = ""),width = 8,height = 6)
plot_cell_trajectory(clustered_spleen_monocle, color_by = "Cell_type",cell_size = 1,cell_link_size = 0.5)+
  scale_color_manual(values = cols)+theme(legend.position = "none")
dev.off()
  
pdf(paste(projectname,"_pseudotime_color_by_","celltype_each.3.","CD4",".pdf",sep = ""),width = 12,height = 6)
plot_cell_trajectory(clustered_spleen_monocle, color_by = "Cell_type",cell_size = 0.2,cell_link_size = 0.5) + facet_wrap(~Cell_type, nrow = 1)
dev.off()

#拟时序匹配到原始数据
load("clustered_spleen_monocle.3.CD4.Rdata")#需要加载夏老师选来保存的数据，保证画出的图一致

clustered_spleen_monocle <- reduceDimension(clustered_spleen_monocle, max_components = 2, reduction_method = "DDRTree",norm_method = 'log',cores=20)
clustered_spleen_monocle <- orderCells(clustered_spleen_monocle)
plot_cell_trajectory(clustered_spleen_monocle, color_by = "State")

cd4_cell_names <- rownames(pData(clustered_spleen_monocle))
cd4_tcell_obj_tmp1 <- subset(x = Tcell_obj, cells = cd4_cell_names)

FeaturePlot(cd4_tcell_obj_tmp1,features = c("TIGIT","ENTPD1","FOXP3","CXCR6","CCR7","IL7R"))
# pdf(file=paste0(project,"_Umap_cd4_tmp1_T_cells.type.3.pdf"),width=12,height=7, onefile=FALSE)
DimPlot(cd4_tcell_obj_tmp1, reduction = "umap",label = TRUE,label.size = 6,
        cols =c("#3B499299","#EE000099","#5F559B99","#A2005699"))
# dev.off()

f <- cbind(cd4_tcell_obj_tmp1@meta.data,pData(clustered_spleen_monocle))
cd4_tcell_obj_tmp1@meta.data$State=factor(f[match(rownames(f),
                                                  rownames(cd4_tcell_obj_tmp1@meta.data)),'State'])

cd4_tcell_obj_tmp1 <- RenameIdents(cd4_tcell_obj_tmp1, `0` = "CD4+Trm cells", `1` = "CD4+Naïve T cells", `2` = "NK cells/Exhausted T cells",`3`="Effector T/Trm cells",`4`="Activated/exhaustion T cells",`5`="Memory/effector T cells",`6`="CD4+Treg cells",`7`="CD4+Naïve /Memory T cells",`8` = "Cytotoxic T cells /NK cells")
cd4_tcell_obj_tmp1  <- AddMetaData(cd4_tcell_obj_tmp1, Idents(cd4_tcell_obj_tmp1), col.name = "Cell_type")
# DimPlot(cd4_tcell_obj_tmp1, reduction = "umap",label = TRUE,label.size = 3,
#         cols =c("#3B499299","#EE000099","#5F559B99","#A2005699"))
cd4_tcell_obj_tmp <- NormalizeData(cd4_tcell_obj_tmp1, normalization.method = "LogNormalize", scale.factor = 10000)
cd4_tcell_obj_tmp <- FindVariableFeatures(cd4_tcell_obj_tmp, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(cd4_tcell_obj_tmp), 10)
all.genes <- rownames(cd4_tcell_obj_tmp)
cd4_tcell_obj_tmp <- ScaleData(cd4_tcell_obj_tmp, features = all.genes)
cd4_tcell_obj_tmp <- RunPCA(cd4_tcell_obj_tmp, features = VariableFeatures(object = cd4_tcell_obj_tmp))
ElbowPlot(cd4_tcell_obj_tmp)

project="cd4_tcell"
cd4_tcell_obj_tmp <- FindNeighbors(cd4_tcell_obj_tmp, dims = 1:20)
cd4_tcell_obj_tmp <- FindClusters(cd4_tcell_obj_tmp, resolution = 0.4)
cd4_tcell_obj_tmp <- RunUMAP(cd4_tcell_obj_tmp, dims = 1:20)

pdf(file=paste0(project,"_Umap_cd4_tmp_T_cells.type.3.pdf"),width=6,height=4, onefile=FALSE)
DimPlot(cd4_tcell_obj_tmp, reduction = "umap",label.size = 6,label = TRUE)
dev.off()


m <- t(clustered_spleen_monocle@reducedDimS) %>% as.data.frame() %>%
  cbind(Cluster=cd4_tcell_obj_tmp@meta.data$RNA_snn_res.0.4) %>% cbind(Pseudotime=clustered_spleen_monocle$Pseudotime)
colnames(m) <- c("Component1","Component2","Cluster","Pseudotime")
# 
cd4_tcell_type <- m %>% ggplot(aes(Component1,Component2,color=Cluster))+geom_point()+
  facet_wrap(~Cluster, nrow = 2)+
  # facet_grid(~Cluster,space = "free_x",shrink = TRUE,scales = "free_x",switch ="x")+
  theme(strip.text.x = element_blank())+
  # theme_classic()+
  # theme(legend.position="none")+
  theme(legend.title = element_blank(),legend.text = element_text(size=8))+
  theme(axis.text.x = element_text(angle=0, hjust=0.8, size=8),
        axis.text.y = element_text(size=8))+
  theme(panel.grid =element_blank()) +
  # theme(axis.ticks.x = element_blank())+
  # theme(axis.ticks.y = element_blank())+
  # theme(axis.line = element_blank())+
  theme(axis.ticks.length=unit(0.1,'cm'))
# scale_x_discrete(expand = c(0,0))+
# scale_y_continuous(expand = c(0,0))
ggsave("CD4+Tcell_types.pdf",plot = cd4_tcell_type,width = 8,height = 4)


ff <- cd4_tcell_obj_tmp1@meta.data[cd4_tcell_obj_tmp1@meta.data$State=="6",]
ff$State <- 0
table(ff$Cell_type)
ff$State[ff$Cell_type=="CD4+Trm cells"] <- "6-1"
ff$State[ff$Cell_type=="CD4+Treg cells"] <- "6-1"
ff$State[ff$Cell_type=="CD4+Naïve T cells"] <- "6-2"
ff$State[ff$Cell_type=="CD4+Naïve /Memory T cells"] <- "6-2"

fg <- cd4_tcell_obj_tmp1@meta.data[cd4_tcell_obj_tmp1@meta.data$State!="6",]
sum(table(fg$State))
ffg <- rbind(ff,fg)
# ffg <- ffg[order(row.names(ffg)), ]
# ffg[row.names(cd4_tcell_obj_tmp1@meta.data)]

cd4_tcell_obj_tmp2 <- subset(x = Tcell_obj, cells = rownames(ffg))
cd4_tcell_obj_tmp2@meta.data$State_type=factor(ffg[match(rownames(ffg),rownames(cd4_tcell_obj_tmp2@meta.data)),'State'])
cd4_tcell_obj_tmp2 <- RenameIdents(cd4_tcell_obj_tmp2, `0` = "CD4+Trm cells", `1` = "CD4+Naïve T cells", `2` = "NK cells/Exhausted T cells",`3`="Effector T/Trm cells",`4`="Activated/exhaustion T cells",`5`="Memory/effector T cells",`6`="CD4+Treg cells",`7`="CD4+Naïve /Memory T cells",`8` = "Cytotoxic T cells /NK cells")
cd4_tcell_obj_tmp2  <- AddMetaData(cd4_tcell_obj_tmp2, Idents(cd4_tcell_obj_tmp2), col.name = "Cell_type")

Idents(cd4_tcell_obj_tmp2)="State_type"
cd4_tcell_obj_tmp2 <- RenameIdents(cd4_tcell_obj_tmp2, `1` = "Treg C1", `3` = "Naïve C2", 
                                   `4` = "Treg C2",`5`="Naïve C3",`6-1`="Treg C4",`6-2`="Naïve C4",
                                   `7`="Naïve C3",`8`="Treg C3",`9`="Naïve C2",`10` = "Naïve C2",`11`="Naïve C1")
cd4_tcell_obj_tmp2 <- AddMetaData(cd4_tcell_obj_tmp2, Idents(cd4_tcell_obj_tmp2), col.name = "cd4_type")

cd4_markers <- FindAllMarkers(cd4_tcell_obj_tmp2,only.pos = TRUE,min.pct = 0.1,logfc.threshold = 0.25)
all.markers <- cd4_markers %>% subset(p_val < 0.05) 
top10 <- all.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
# write.table(top10,"CD4_tcell_deg_top10.txt",sep="\t",quote=F,row.names=F)
pdf("top10_DoHeatmap.pdf",width = 10,height = 8,onefile = FALSE)
DoHeatmap(cd4_tcell_obj_tmp2,features = top10$gene)+NoLegend()
dev.off()

dot_top10_p <- DotPlot(cd4_tcell_obj_tmp2,features = unique(top10$gene))+ RotatedAxis()+
  scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)
dot_top10_p$data$id <- factor(dot_top10_p$data$id,levels = c("Treg C1","Treg C2","Treg C3","Treg C4","Naïve C1", "Naïve C2","Naïve C3","Naïve C4"))

ggsave("cd4_top10_Dotplot.pdf",plot=dot_top10_p,width = 20,height = 8)

pdf("cd4_cluster.pdf",width = 12,height = 8,onefile = FALSE)
# DimPlot(cd4_tcell_obj_tmp2,reduction = "umap",label = TRUE,label.size = 4,
#         cols =c("#FF0000","#FF8247","#FFAEB9","#EE00EE","#40E0D0","#2E8B57","#FFFF00","#8B4513"))
cd4_cluster <- cd4_tcell_obj_tmp2@reductions$umap@cell.embeddings %>% as.data.frame() %>% 
  cbind(cd4_type = cd4_tcell_obj_tmp2@meta.data$cd4_type)

cols =c("#FF0000","#FF8247","#FFAEB9","#EE00EE","#40E0D0","#2E8B57","#FFFF00","#8B4513")
cd4_cluster$cd4_type=factor(cd4_cluster$cd4_type,levels = c("Treg C1","Treg C2","Treg C3","Treg C4","Naïve C1", "Naïve C2","Naïve C3","Naïve C4"))
ggplot(cd4_cluster,aes(x=UMAP_1,y=UMAP_2,color = cd4_type))+
  geom_point(size=2)+
  # stat_summary(data = dplyr::filter(umap_low_tcr,grepl('TCR',All_Low_infiltration_TCR)), geom = "point")+
  # stat_summary(data = dplyr::filter(umap_low_tcr,grepl('clonotype1',All_Low_infiltration_TCR)), geom = "point")+
  # stat_summary(data = dplyr::filter(umap_low_tcr,grepl('clonotype2',All_Low_infiltration_TCR)), geom = "point")+
  # stat_summary(data = dplyr::filter(umap_low_tcr,grepl('clonotype3',All_Low_infiltration_TCR)), geom = "point")+
  # stat_summary(data = dplyr::filter(umap_low_tcr,grepl('clonotype4',All_Low_infiltration_TCR)), geom = "point")+
  # stat_summary(data = dplyr::filter(umap_low_tcr,grepl('clonotype5',All_Low_infiltration_TCR)), geom = "point")+
  scale_color_manual(values = cols)+
  theme_classic()+
  theme(panel.grid =element_blank()) +
  theme(axis.ticks.x = element_blank())+
  theme(axis.ticks.y = element_blank())+
  # theme(axis.ticks.length=unit(0.1,'cm'))+
  theme(legend.title = element_blank(),legend.text = element_text(size=20))+
  theme(axis.title.x=element_text(vjust=1, size=20,color = "black"),  # X axis title
        axis.title.y=element_text(size=20,color = "black"),
        axis.text.x = element_text(angle=0, hjust=1, vjust=0.5, size=20),
        axis.text.y = element_text(size=20),
        axis.line = element_line(color = "black",linetype = 1,size = 1))
dev.off()

###treg
tregc1 = c("DUSP4","FOXP3","TNFRSF18")
tregc1_p1 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "DUSP4",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 1.5,max.cutoff = 2)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
tregc1_p2 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "FOXP3", cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 1.5,max.cutoff = 2)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
tregc1_p3 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "TNFRSF18",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 1.5,max.cutoff = 2)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()

tregc2 = c("MVD","MPC1")
tregc2_p1 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "MVD",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 1.5,max.cutoff = 2)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
tregc2_p2 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "MPC1",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 1.5,max.cutoff = 2)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()

tregc3 = c("CD52","MT1E","MT2A")
tregc3_p1 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "CD52",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 3,max.cutoff = 4)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
tregc3_p2 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "MT1E",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 3,max.cutoff = 4)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
tregc3_p3 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "MT2A",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 3,max.cutoff = 4)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()


tregc4 = c("CXCR4","ZFP36")
tregc4_p1 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "CXCR4",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 4,max.cutoff = 5)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
tregc4_p2 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "ZFP36",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 4,max.cutoff = 5)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()


###naive
naivec1 = c("RSRP1","AHNAK","NEAT1")
naivec1_p1 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "RSRP1",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 2.5,max.cutoff = 3)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend() 
naivec1_p2 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "AHNAK",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 2.5,max.cutoff = 3)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
naivec1_p3 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "NEAT1",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 2.5,max.cutoff = 3)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
  
naivec2 = c("CCL5","NKG7")
naivec2_p1 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "CCL5",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 2,max.cutoff = 3)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
naivec2_p2 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "NKG7",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff =1,max.cutoff = 3)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()

naivec3 = c("DUSP1","DNAJB1")
naivec3_p1 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "DUSP1",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 3.5,max.cutoff = 4)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
naivec3_p2 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "DNAJB1",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 3.5,max.cutoff = 4)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()

naivec4 = c("RPS3A","RPS8","CCR7")
naivec4_p1 <- FeaturePlot(cd4_tcell_obj_tmp2,features = c("RPS3A"),cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 4,max.cutoff = 5)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
naivec4_p2 <- FeaturePlot(cd4_tcell_obj_tmp2,features = c("RPS8"),cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 4,max.cutoff = 5)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
naivec4_p3 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "CCR7",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 1,max.cutoff = 3)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()

marker_umap <- ggarrange(tregc1_p1,tregc1_p2,tregc1_p3,tregc2_p1,tregc2_p2,tregc3_p1,tregc3_p2,tregc3_p3,tregc4_p1,tregc4_p2,
          naivec1_p1,naivec1_p2,naivec1_p3,naivec2_p1,naivec2_p2,naivec3_p1,naivec3_p2,naivec4_p1,naivec4_p2,
          naivec4_p3,ncol = 5,nrow = 4)
ggsave("cd4_marker_umap.pdf",plot=marker_umap,width = 15,height = 10)

gene = c("DUSP4","FOXP3","TNFRSF18","MVD","MPC1","CD52","MT1E","MT2A","CXCR4","ZFP36","RSRP1","AHNAK",
         "NEAT1","CCL5","NKG7","DUSP1","DNAJB1","RPS3A","RPS8","CCR7")
dot_p <- DotPlot(cd4_tcell_obj_tmp2, features = gene)+ RotatedAxis()+scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0) 
dot_p$data$id <- factor(dot_p$data$id,levels = c("Treg C1","Treg C2","Treg C3","Treg C4","Naïve C1", "Naïve C2","Naïve C3","Naïve C4"))
ggsave("cd4_Dotplot.pdf",plot=dot_p,width = 10,height = 5)


DUSP4 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "DUSP4",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 1.5,max.cutoff = 2)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
MPC1 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "MPC1",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 1.5,max.cutoff = 2)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
CD52 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "CD52",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 3,max.cutoff = 4)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
HPGD <- FeaturePlot(cd4_tcell_obj_tmp2,features = "HPGD",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 2,max.cutoff = 3)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()

AHNAK <- FeaturePlot(cd4_tcell_obj_tmp2,features = "AHNAK",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 2.5,max.cutoff = 3)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
GZMK <- FeaturePlot(cd4_tcell_obj_tmp2,features = "GZMK",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 2.5,max.cutoff = 3)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
DUSP1 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "DUSP1",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 3.5,max.cutoff = 4)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
FOS <- FeaturePlot(cd4_tcell_obj_tmp2,features = "FOS",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 3.5,max.cutoff = 4)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()

CCR7 <- FeaturePlot(cd4_tcell_obj_tmp2,features = "CCR7",cols = c("lightgrey" ,"#DE1F1F"),min.cutoff = 1,max.cutoff = 3)+
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()

marker_umap <- ggarrange(DUSP4,MPC1,CD52,HPGD,AHNAK,GZMK,FOS,CCR7,ncol = 4,nrow = 2)
ggsave("cd4_marker_umap_1.pdf",plot=marker_umap,width = 10,height = 4)
marker_umap <- ggarrange(DUSP4,MPC1,CD52,HPGD,AHNAK,GZMK,DUSP1,CCR7,ncol = 4,nrow = 2)
ggsave("cd4_marker_umap_2.pdf",plot=marker_umap,width = 10,height = 4)

marker_umap <- ggarrange(DUSP4,MPC1,CD52,HPGD,AHNAK,GZMK,FOS,DUSP1,CCR7,ncol = 3,nrow = 3)
ggsave("cd4_marker_umap_3.pdf",plot=marker_umap,width = 8,height = 6)


gene = c("DUSP4","MPC1","CD52","HPGD","AHNAK","GZMK","FOS","DUSP1","CCR7")
dot_p <- DotPlot(cd4_tcell_obj_tmp2, features = gene)+ RotatedAxis()+scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0) 
dot_p$data$id <- factor(dot_p$data$id,levels = c("Treg C1","Treg C2","Treg C3","Treg C4","Naïve C1", "Naïve C2","Naïve C3","Naïve C4"))
ggsave("cd4_Dotplot_1.pdf",plot=dot_p,width = 10,height = 5)

