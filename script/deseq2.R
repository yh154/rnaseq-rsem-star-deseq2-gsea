log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

rnk <- snakemake@params[["rnk"]]
adjusted_pvalue <- snakemake@params[["adjusted_pvalue"]]
gtf <- snakemake@params[["gtf"]]

dir.create("diffexp")

################
# functions
################
fun_pca <- function(dds, title=""){
  require("ggplot2")
  require("DESeq2")
  vsd <- vst(dds,blind=TRUE);
  pcaData <- plotPCA(vsd, intgroup = c( "condition"),returnData=TRUE);
  percentVar <- round(100 * attr(pcaData, "percentVar"));
  ggplot(pcaData,aes(PC1, PC2, color=condition)) + geom_point(size=3) +
    guides(fill=FALSE) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text(aes(label=name), vjust= -0.5, hjust= -0.1) + 
    geom_vline(xintercept = 0, linetype="dashed",color="gray50") + 
    geom_hline(yintercept = 0, linetype="dashed",color="gray50") +
    ggtitle(title) + theme_bw() + coord_fixed() +
    theme(plot.title = element_text(size = 12, face = "bold",hjust=0.5)
          ,panel.grid.minor = element_blank()
          #,panel.grid.major = element_blank()
          ,aspect.ratio = 1)
}
fun_plot_padj_cnt <- function(result, title=''){
  par(cex=.8) # control bar top label size
  hist(result[["padj"]],
       main = title,
       xaxt='n',
       cex.lab=1.2, # control xlab/ylab size
       breaks = 0:20/20,col = "steelblue",
       border = "white", labels = T,
       xlab = "Adjusted-pvalue",ylab="Gene count")
       axis(side=1, at=0:10/10, labels=0:10/10)
}

fun_write_result_csv <- function(result, file_id="output"){
  if(nrow(result)>0){
    fout <- sprintf("diffexp/Deseq2_res_%s.csv",file_id)
    write.csv(result, fout, quote=F, row.names=FALSE)
  }else{
    message(sprintf("%s result table is empty. There might be an issue with GTF file.", file_id))
  }
}

fun_process_gtf <- function(gtf){
  gtf <- read.table(gtf, sep="\t");
  if(ncol(gtf) >= 9){
    gtf <- gtf[gtf$V3=="gene",]
    V9 <- t(sapply(as.character(gtf$V9),function(x){
            y=as.character(unlist(strsplit(x,"; ")))
            idx_id = grepl("gene_id", y)
            idx_type = grepl("gene_type|gene_biotype", y)
            idx_name = grepl("gene_name",y)
            c(y[idx_id],y[idx_type],y[idx_name])
          }))
    if(ncol(V9)==3){
        V9 <- apply(V9,2,function(x){gsub("gene_id\\s+|gene_biotype\\s+|gene_type\\s+|gene_name\\s+","",x)})
        colnames(V9)=c("gene_id","gene_type","gene_name")
        return(V9)
    }else{stop("GTF has fewer columns as expected.")}
  }else{
    stop("GTF has fewer columns as expected.")
  }
}

# add deseq normalized count to result
fun_add_ncnt <- function(result, file_name, dds){
  ncnt <- DESeq2::counts(dds, normalized=TRUE)
  sel <- as.character(unlist(strsplit(file_name, '_vs_')))
  coldata <- SummarizedExperiment::colData(dds)
  within <- coldata[["condition"]] %in% sel
  if(sum(within)>1){
    result <- merge(ncnt[,within], result, by='row.names')
  }else{
    message(sprintf("Error in selecting columns from normalized count table for %s", file_name))
  }
}

# add gene symbol to de result table
# gtf - ensembl format, containing 'gene_id','gene_type' & 'gene_name' in column 9.
fun_add_symbol <- function(result, gtf){
   merged <- merge(gtf, result, by=1, all.y=TRUE)
}

fun_write_rnk <- function(rnk, file_id="output"){
  fout <- sprintf("diffexp/gsea_%s.rnk", file_id)
  write.table(rnk, fout, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# colData and countData must have the same sample order
contrasts = NULL
cnt<- read.table(snakemake@input[["counts"]], sep="\t", header=TRUE, row.names=1, check.names=FALSE)
coldata <- read.table(snakemake@params[["coldata"]], sep="\t", header=TRUE, row.names=1, check.names=FALSE)
#if(nrow(coldata)==0 | ncol(coldata)<4){
#    coldata <- read.csv(snakemake@params[["coldata"]], header=TRUE, row.names=1, check.names=FALSE)
#}

colnames(cnt) = gsub("_expected_count$","",colnames(cnt))

print(head(cnt))
print(coldata)

if(!all(rownames(coldata)==colnames(cnt))){
  od=match(rownames(coldata), colnames(cnt))
  cnt = cnt[,od]
}

# remove low expressed genes
#cnt <- cnt[rowSums(sign(cnt))>1,]
cnt <- cnt[rowMeans(cnt)>1,]
# make sure counts in integer
rns <- rownames(cnt)
cnt <- apply(cnt,2,as.integer)
rownames(cnt) <- rns

dds <- DESeq2::DESeqDataSetFromMatrix(countData=cnt, colData=coldata, design=~ condition)

# normalization and preprocessing
dds <- DESeq2::DESeq(dds)
saveRDS(dds, file=snakemake@output[[1]])
if(is.null(contrasts)){
    contrasts <- levels(SummarizedExperiment::colData(dds)[["condition"]])
    contrasts <- combn(contrasts, 2, simplify=FALSE)
}else{
    contrasts <- as.list(contrasts)
    contrasts <- lapply(contrasts, function(x){as.character(unlist(strsplit(x,"_vs_")))}
}
res <- lapply(contrasts, function(x){
  DESeq2::results(dds, contrast = c("condition",x[1],x[2]), alpha = adjusted_pvalue, independentFiltering =T)
})
names(res) <- lapply(contrasts, paste, collapse="_vs_")

# plotting, before any table merging.
pdf("diffexp/diffexp.pdf", height=8.5, width=11)
fun_pca(dds)
mapply(DESeq2::plotMA, res, main=names(res))
mapply(fun_plot_padj_cnt, res, names(res))
dev.off()

# write de result table
# add deseq2 normalized counts
res <- mapply(fun_add_ncnt, res, names(res), MoreArgs = list(dds=dds), SIMPLIFY = FALSE)
# check if gene identify is ensembl, add symbol to result

# add symbol if rowname is ensembl id.
if(all(grepl("^ENS", rns))){
     gtf <- fun_process_gtf(gtf)
     res <- lapply(res, fun_add_symbol, gtf=gtf)
 }else{
     res <- lapply(res,function(x){colnames(x)[1]="gene_name";x})
 }

mapply(fun_write_result_csv, res, file_id=names(res))

if(toupper(rnk)){
  rnk <- lapply(res, function(x){
    x <- x[,c("gene_name","log2FoldChange")]
    x <- x[!is.na(x[["log2FoldChange"]]),]
    x <- x[with(x, order(-abs(x["log2FoldChange"]))), ]
    x <- x[!duplicated(x[["gene_name"]]), ]
    x[["gene_name"]] <- toupper(x[["gene_name"]])
    x
  })
  mapply(fun_write_rnk, rnk, names(rnk))
}
