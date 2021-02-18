#DESEQ2 code 
R
library('DESeq2')
coldata<-read.csv('final_counts_combined_high_quality.csv', row.names=1)
rownames(coldata)<-sub(".v5.1", "", rownames(coldata))
names(coldata)<-sub("LIB.", "", names(coldata))
names(coldata)<-sub(".1_counts.tsv", "", names(coldata))
names(coldata)<-sub("_merged_samples.counts.tsv", "", names(coldata))
cts<-coldata
coldata<-read.csv('coldata_high_quality.csv')
dds<-DESeqDataSetFromMatrix(countData=cts, colData=coldata, design= ~Sex)
dds<-DESeq(dds)
res<-results(dds)
plotMA(res, ylim=c(-15,15))
res1<-merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by='row.names', sort=FALSE)
write.csv(res1, 'DESEQ2 M vs F high quality with counts.csv')

ntd<-normTransform(dds)
library('pheatmap')
select<-order(rowMeans(counts(dds)), decreasing=TRUE)[1:5000]
df<-as.data.frame(colData(dds))
pdf('top_5000_DE_genes_heatmap.pdf', width=24, height=12)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE)
dev.off()

select<-order(res$log2FoldChange), decreasing=TRUE)[1:1000]
pdf('top_1000_DE_genes_heatmap.pdf', width=24, height=12)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE)
dev.off()

select<-order(res$log2FoldChange), decreasing=FALSE)[1:1000]
pdf(width=24, height=12)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE)
dev.off()


Auxin<-read.csv('Auxin.csv', header=TRUE)
Cytokinin<-read.csv('Cytokinin.csv', header=TRUE)
GA<-read.csv('GA.csv', header=TRUE)
Jasmonate<-read.csv('Jasmonate.csv', header=TRUE)
Ethylene<-read.csv('Ethylene.csv', header=TRUE)
RDDM<-read.csv('RDDM.csv', header=TRUE)
MADSbox<-read.csv('MADS-box.csv', header=TRUE)

select<-which(rownames(assay(dds)) %in% Auxin$transcriptID)
pdf('auxin_DE.pdf', width=24, height=12)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE)
dev.off()

select<-which(rownames(assay(dds)) %in% Cytokinin$transcriptID)
pdf('cytokinin_DE.pdf', width=24, height=12)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE)
dev.off()

select<-which(rownames(assay(dds)) %in% GA$transcriptID)
pdf('gibberellic_acid_DE.pdf', width=24, height=12)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE)
dev.off()

select<-which(rownames(assay(dds)) %in% Ethylene$transcriptID)
pdf('ethylene_DE.pdf', width=24, height=12)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE)
dev.off()

select<-which(rownames(assay(dds)) %in% Jasmonate$transcriptID)
pdf('Jasmonate_DE.pdf', width=24, height=12)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE)
dev.off()

select<-which(rownames(assay(dds)) %in% RDDM$transcriptID)
pdf('RDDM_DE.pdf', width=24, height=12)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE)
dev.off()

select<-which(rownames(assay(dds)) %in% MADSbox$transcriptID)
pdf('MADSbox_DE.pdf', width=24, height=12)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE)
dev.off()












