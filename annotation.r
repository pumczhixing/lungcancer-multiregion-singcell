suppressMessages(library('Seurat'))
args=commandArgs(T)
obj=args[1]
obj = eval(parse(text = load(obj)))
outdir=args[2]
if(!dir.exists(outdir)){dir.create(outdir)}
setwd(outdir)

new.cluster.ids <- c(rep("Alveolar",2),"B cell","Endothelial","Fibroblast",
rep("T cell",6),rep("Myeloid",8),rep("Epithelial",3),"Mast cell")

names(new.cluster.ids) <- c(17,19,10,15,14,3,2,0,18,8,16,1,20,11,5,13,4,22,21,6,7,9,12)
obj <- RenameIdents(obj, new.cluster.ids)

pdf("tsne_anno.pdf",width=10,height=8)
p1=DimPlot(obj, reduction = "tsne",label = FALSE,pt.size=0.1)
print(p1)
dev.off()

marker=c("CLDN18","FOLR1","AQP4","PEBP4",#Alveolar
"CD79A","IGKC","IGLC3","IGHG3",#B cell
"CLDN5","FLT1","CDH5","RAMP2",#Endothelial
"COL1A1","DCN","COL1A2","C1R",#Fibroblast
"CD3D","CD3E","CD3G","CD7","CD2","TRBC1","TRBC2",# T cell
"LYZ","CD68","FCGR3A","MARCO",#Myeloid:
"KIT",#Mast cell
"EPCAM","TMEM190","CAPS","PIFO","SNTN")#Epithelial

obj <- ScaleData(obj, features=marker,vars.to.regress = c("nCount","percent.mt"))

pdf("marker_2.pdf",height=8,w=10)
DoHeatmap(obj, features = marker) + NoLegend()
dev.off()

jpeg("vlnplot.jpeg",width = 480*6, height = 480*3)
VlnPlot(obj, marker,ncol=4)
dev.off()

pdf("feature_tSNE.pdf", w = 8, h = 6)
for(x  in marker) {
print(FeaturePlot(obj, features = x,pt.size=0.2))
}
dev.off()

save(obj,file="obj_cellanno.RData")
