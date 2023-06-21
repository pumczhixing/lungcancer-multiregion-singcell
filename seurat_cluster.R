log <- file(snakemake@log[[1]], open="wt")#打开读写
sink(log) #.sink函数将输出结果重定向到文件。
sink(log, type="message")#.sink函数将输出结果重定向到文件。
suppressMessages(library(Seurat))

Cal_npc = function(scRNA.PCA){
    findKneePoint <- function(pcs)
    {
        npts <- length(pcs)
        if(npts<=3){
            return(npts)
        }else{
            P1 <- c(1,pcs[1])
            P2 <- c(npts,pcs[npts])
            v1 <- P1 - P2
            dd <- sapply(2:(npts-1),function(i){
                Pi <- c(i, pcs[i])
                v2 <- Pi - P1
                m <- cbind(v1,v2)
                d <- abs(det(m))/sqrt(sum(v1*v1))
            })
            return(which.max(dd))
        }
    }
    ### find elbow point and get number of components to be used
    eigenv.prop <- (scRNA.PCA@reductions$pca@stdev^2)/sum(scRNA.PCA@reductions$pca@stdev^2)
    eigengap <- eigenv.prop[-length(eigenv.prop)]-eigenv.prop[-1]
    
    ### method 1
    ######pca.res$kneePts <- which(pca.res$eigengap<1e-4)[1]
    
    ### method 2 (max distance to the line defined by the first and last point in the scree plot)
    kneePts <- findKneePoint(head(eigenv.prop, n=100)) + 5
    
    if(!is.na(kneePts)){ pca.npc <- kneePts}
    pca.npc <- min(pca.npc,ncol(scRNA.PCA@reductions$pca@cell.embeddings))
    sprintf("set pca.npc to %d while kneePts is at %d (ssc.reduceDim)\n",pca.npc,
    
    if(!is.na(kneePts)) kneePts else -1)
    eigenv <- eigenv.prop * 100
    dat.plot.eigenv <- data.frame(PC=seq_along(eigenv),
                                eigenv=eigenv,
                                isKneePts=as.character(seq_along(eigenv)==kneePts),
                                stringsAsFactors = F)

    # p1=ggplot2::ggplot(head(dat.plot.eigenv,n=100),mapping = aes(PC,eigenv)) +
    # geom_point(aes(colour=isKneePts),show.legend=F) + ylab("Variation explained (%)") +
    # scale_colour_manual(values = c("TRUE"="#E41A1C","FALSE"="#377EB8")) +
    # theme_bw()
    return(list("pca.npc"=pca.npc,"dat.plot.eigenv"=dat.plot.eigenv))
}

obj <- eval(parse(text = load(snakmake@input[[1]])))
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj,selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj,vars.to.regress = c("nCount_RNA", "percent.mt"))
obj <- RunPCA(obj,features = VariableFeatures(object = obj))

npc = Cal_npc(obj)
jpeg("ElbowPlot.jpeg")
ggplot2::ggplot(head(npc[["dat.plot.eigenv"]],n=100),mapping = aes(PC,eigenv)) +
    geom_point(aes(colour=isKneePts),show.legend=F) + ylab("Variation explained (%)") +
    scale_colour_manual(values = c("TRUE"="#E41A1C","FALSE"="#377EB8")) +theme_bw()
dev.off()

ndim = npc[["pca.npc"]]  # ndim = 10

obj <- FindNeighbors(obj, dims = 1:ndim)  # obj@commands$FindNeighbors.RNA.pca
obj <- FindClusters(obj, resolution = 0.5)

obj <- RunTSNE(obj, dims = 1:ndim)

# embed_tsne <- Embeddings(obj, 'tsne')   #提取tsne图坐标
# write.csv(embed_tsne,'embed_tsne.csv')

obj <- RunUMAP(obj, dims = 1:20)
# embed_umap <- Embeddings(obj, 'umap')   #提取umap图坐标
# write.csv(embed_umap, 'embed_umap.csv')

save(obj, file = snakemake@output[[1]])  # 5_cluster/seurat_cluster.RData
                                         # 5_final/seurat_cluster.3.RData with UMAP


