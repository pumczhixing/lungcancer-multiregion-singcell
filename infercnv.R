
library(futile.logger)
flog.info("::InferCNV program:Start")
flog.info("::Load function:Start")
suppressMessages(library('Seurat'))
suppressMessages(library(ggplot2))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
source("/glusterfs/home/local_xia_hao/project_202007/scRNA/Guo_wei/New_sample_20200811/scRNA/src/inferCNV.R")
flog.info("::Load input file:Start")
gene=read.delim("/titan2/liu_yuhao0101/database/inferCNV/gene.GRCH38.all.gene.txt",row.names=1)
object.seurat=eval(parse(text = load("./All_samples_obs_for_cnv.RData")))
obj_not_immune <-subset(obj,idents=c("Alveolar cells","Myeloid cells","Endothelial cells","Epithelial cells","Mast cells","Fibroblast cells"))
object.seurat <- obj_not_immune
if(length(setdiff(row.names(object.seurat),row.names(gene)))!=0){
  stop("Not all gene was in gene order files")
}
#########filter CHRX & CHRY################################################
gene=subset(gene,is.element(gene[,1],paste("chr",1:22,sep=""))==TRUE)
gene[,1]=factor(gene[,1],levels=paste("chr",1:22,sep=""))
counts=GetAssayData(object =object.seurat, assay = "RNA",slot = 'counts')
counts_matrix = as.matrix(counts[(rownames(counts) %in% row.names(gene)),])

keep.genes=intersect(rownames(counts_matrix),row.names(counts_matrix)[rowSums(counts_matrix!=0) >0.05*ncol(counts_matrix)])
#keep.genes=intersect(rownames(counts_matrix),names(rowSums(counts_matrix==0) >0.1*ncol(counts_matrix)))
if(length(keep.genes)<5000){
  warning("genes used for infer CNV in %d",length(keep.genes))
}
print(length(keep.genes))
counts_matrix=counts_matrix[keep.genes,]
cpm=sweep(counts_matrix, STATS=colSums(counts_matrix), MARGIN=2, FUN="/")*1e6
cp100k = cpm/10
flog.info("::Prepare input file for infercnv:Start")
expression_data = as.data.frame(log2(cp100k+1))
input_gene_order=read.delim("/titan2/liu_yuhao0101/database/inferCNV/gene.GRCH38.all.gene.txt",row.names=1)
CD45.cell=colnames(counts[,counts["PTPRC",]!=0])
input_reference_samples=CD45.cell
order_ret <- order_reduce(data=expression_data,
                                    genomic_position=input_gene_order)
expression_data <- order_ret$expr
input_gene_order <- order_ret$order
flog.info("::infercnv:Start")
ret_list = infer_cnv(data=expression_data,
                               gene_order=input_gene_order,
                               cutoff=0,
                               reference_obs=input_reference_samples,
                               transform_data=FALSE,
                               window_length=101,
                               max_centered_threshold=3,
                               noise_threshold=0,
                               num_ref_groups=1,
                               out_path=".",
                               k_obs_groups=1,
                               plot_steps=FALSE,
                               contig_tail=0,
                               method_bound_vis="average_bound",
                               lower_bound_vis=-1,
                               upper_bound_vis=1,
                               ref_subtract_method="by_mean",
                               hclust_method='complete')
flog.info("::infercnv:End")
save(ret_list,file="ret_list.rData")

reference_idx=ret_list[["REF_OBS_IDX"]]
obs_data0=data.frame(ret_list[["VIZ"]])
colnames(obs_data0)=gsub("\\.","-",colnames(obs_data0))
ref_data=obs_data0[,colnames(obs_data0) %in% reference_idx]

obs_data=obs_data0[,!(colnames(obs_data0) %in% reference_idx)]
cnv.single       = apply(obs_data,2,function(x){mean(x*x)})
cnv.top5.cell    = names(head(sort(cnv.single,decreasing=TRUE),ceiling(0.05*ncol(obs_data))))
cnv.top5.avgcnv  = apply(obs_data[,cnv.top5.cell],1,mean)
cor.matrix_top5=cor(obs_data0,cnv.top5.avgcnv )
scatter.m=data.frame("CNV.signal"=apply(obs_data0,2,function(x){mean(x*x)}),"CNV.correlation"=cor.matrix_top5[,1])

scatter.m$state=sapply(row.names(scatter.m),function(x) {
  if(scatter.m[x,1]>=0.03 & scatter.m[x,2]>=0.4){
    "malignant cells"
  }else if(scatter.m[x,1]<0.03 & scatter.m[x,2]<0.4){
    "non-malignant cells"
  }else{
    "unresolved cells"
  }
})
#write.table(scatter.m,snakemake@output[[1]],sep="\t")
#write.table(scatter.m,snakemake@output[[1]],sep="\t")

write.table(scatter.m,"scatter.m.table.xls",sep="\t",quote=F)
#outdir=gsub("_si.*","",])
outdir="Temp"

if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
setwd(outdir)



pdf("cnv_scatter.pdf")
ggplot(data=scatter.m, aes(x=CNV.signal, y=CNV.correlation, color=state))+geom_point(size=0.2)+
  scale_y_continuous(limits=c(-0.6,1)) +scale_x_continuous(limits=c(0,0.2)) +#坐标轴翻转由coord_flip()实现
  theme_bw()+labs(title="CNV inference",y="CNV correlation",x="CNV signal")#+
dev.off()

malignant.cell=row.names(scatter.m[scatter.m$stat %in% "malignant cells",])
non.malignant.cell=row.names(scatter.m[scatter.m$stat %in% "non-malignant cells",])
gene=read.delim("/titan2/liu_yuhao0101/database/inferCNV/gene.GRCH38.all.gene.txt",row.names=1)
gene=gene[row.names(gene) %in% row.names(obs_data0),]
as.numeric(gsub("chr","",gene[,1])) ->gene[,1]
gene=gene[order(gene[,1]),]
obs_data0=obs_data0[row.names(gene),]
non.malignant.cnv=obs_data0[,non.malignant.cell]
malignant.cnv=obs_data0[,malignant.cell]
save(obs_data0,file="All_samples_obs_data0.CNV.RData")
save(ref_data,file="ref_data.CNV.RData")
save(obs_data,file="obs_data.CNV.RData")

ha = HeatmapAnnotation(df=factor(gene$geneID,levels = 1:22),
                       col = list(df=c(
                         "1" = "black","2" = "green", "3" = "grey","4" = "black",
                         "5" = "green","6" = "grey","7" = "black", "8" = "green",
                         "9" = "grey","10" = "black", "11" = "green","12" = "grey",
                         "13" = "black","14" = "green","15" = "grey",  "16" = "black",
                         "17" = "green","18" = "grey","19" = "black", "20" = "green",
                         "21" = "grey","22" = "black")))
p1=Heatmap(t(malignant.cnv),cluster_rows = T,cluster_columns = F,
           show_column_names = F,show_row_names = FALSE,
           col=colorRamp2(c(-1, 0, 1), c("midnightblue", "white", "darkred")),
           top_annotation = ha)
p2=Heatmap(t(non.malignant.cnv),cluster_rows = T,cluster_columns = F,
           show_column_names = F,show_row_names = FALSE,
           col=colorRamp2(c(-1, 0, 1), c("midnightblue", "white", "darkred")),
           top_annotation = ha)
pdf("malignant.cnv.pdf")
print(p1)
dev.off()
pdf("non.malignant.cnv.pdf")
print(p2)
dev.off()
#save(obs_data0,file="All_samples_obs_data0.CNV.RData")
#save(ref_data,file="ref_data.CNV.RData")
#save(obs_data,file="obs_data.CNV.RData")
flog.info("::InferCNV program:complete")


#flist = list.files(".",pattern="scatter.m.table.xls",recursive =T)
mm.cell.index={}
unresolved.cell.index={}
#summary=read.delim(flist[1])
summary=read.table("../scatter.m.table.2.xls",sep="\t",header=T)

#metadata=obj@meta.data
mm.cell.index[[dirname( )]]=row.names(summary)[row.names(metadata[metadata$Sample %in% dirname(flist[1]),]) %in% summary$state %in% "malignant cells"]
for(index in 2:length(flist)){
    fname=read.delim(flist[index])
     summary=rbind(summary,fname)
 }

write.table(summary,"summary.inferCNV.xls",sep="\t")
obj = eval(parse(text = load("../3_option_anno/obj_cellanno.RData")))
suppressMessages(library('Seurat'))
obj=obj[,row.names(obj@meta.data)]
obj <- AddMetaData(object = obj, metadata = Idents(obj), col.name = "cellType")
cellIndex = row.names(obj@meta.data)
test=factor(
    unlist(sapply(cellIndex,function(x) {
        ifelse(x %in% row.names(summary),as.character(summary[gsub("\\.","-",x),"state"]),"non-malignant cells")})),
        levels=c("malignant cells","non-malignant cells","unresolved cells"))
obj <- AddMetaData(object = obj, metadata = test, col.name = "typeCNV")
head(obj@meta.data)
save(obj,file="obj.cnv.anno.RData")
pdf("tSNE_cnv_1.pdf", w = 8, h = 6)
print(DimPlot(object = obj, group.by="typeCNV",dims=c(1,2),reduction = "tsne",pt.size=0.1))
dev.off()






