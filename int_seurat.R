log <- file(snakemake@log[[1]], open="wt")#打开读写
sink(log) #.sink函数将输出结果重定向到文件。
sink(log, type="message")#.sink函数将输出结果重定向到文件。
suppressMessages(library('Seurat'))
addInfo = function(obj,data){
    #obj: seurat
    #data:dataframe sample flowed samples.csv
    cellIndex = row.names(obj@meta.data)
    obj= obj[,cellIndex]
    if(length(colnames(data))>=2){
        for(type in colnames(data) ){
            test=unlist(sapply(cellIndex,function(x) data[gsub(".*\\.","",x),type]))
            obj <- AddMetaData(object = obj, metadata = test, col.name = type)
        }
    }
    return(obj)
}
html=snakemake@input[[1]]
outdir=dirname(snakemake@output[[1]])
fname=paste(dirname(html),"filtered_feature_bc_matrix.h5",sep="/")
obj=Read10X_h5(fname)
colnames(obj)=paste(gsub("-.*","",colnames(obj)),".",snakemake@params[[1]],sep="")
obj <- CreateSeuratObject(obj, min.cells = 3, min.features = 200)
sprintf("raw cell number is %s",ncol(obj))
#QC
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-|^MT-")

#3.filter data
#high.nUMI <- isOutlier(sce$pct_counts_Mito, nmads=4, type="higher")
## method 1
## obj <- subset(obj,subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)
## method 2
max_mt   = attr(scater::isOutlier(obj@meta.data$percent.mt,nmads=4,type="higher"),"thresholds")[2]
max_gene = attr(scater::isOutlier(obj@meta.data$nFeature_RNA,nmads=4,type="higher"),"thresholds")[2]
obj <- subset(obj,subset = nFeature_RNA > 200 & nFeature_RNA < max_gene & percent.mt < max_mt)
sprintf("filtered cell number is %s",ncol(obj))

obj[['max_gene']] = max_gene
obj[['max_mt']] = max_mt
Info = read.delim("samples.csv",sep=",",header=T)
rownames(Info)=Info$sample
Info = Info[,!(colnames(Info) %in% "path")]
sampleInfo = Info[snakemake@params[[1]],]
obj=addInfo(obj,sampleInfo)
head(obj@meta.data)
save(obj,file=snakemake@output[[1]])  # 4_seurat_int/*_seurat_filtered.RData
