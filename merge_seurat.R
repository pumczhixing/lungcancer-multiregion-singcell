log <- file(snakemake@log[[1]], open="wt")#打开读写
sink(log) #.sink函数将输出结果重定向到文件。
sink(log, type="message")#.sink函数将输出结果重定向到文件。
#obj_seurat=eval(parse(text = load(obj)))
print(snakemake@input)
objlst=unlist(snakemake@input)
suppressMessages(library('Seurat'))
num=length(objlst)
if(length(objlst)<2){
    cat(" obj should be offered at least 2")
}else{
    obj0 = eval(parse(text = load(objlst[1])))
    for(x in 2:num){
        temp=objlst[x]
        objy=eval(parse(text = load(temp)))
        obj0=merge(x=obj0,y=objy)
    }
    save(obj0,file=snakemake@output[[1]])
}
nrow(obj0@meta.data)   # 4_seurat_int/merge.RData
