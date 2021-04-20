library(monocle)
library(umap)
dataATAC <- read.table("GSE74912_ATACseq_All_Counts.txt")
samplesATAC <- read.table("sampleInfoATAC.txt")

rownames(samplesATAC) <- samplesATAC[,1]
samplesATAC <- samplesATAC[ ,-1]

colnames(samplesATAC) <- samplesATAC[1, ]
colnames(dataATAC) <- dataATAC[1, ]
dataATAC <- dataATAC[-1, ]

rownames(dataATAC) <- c(1:590650)
featuresATAC <- data.frame(gene_short_name=paste("Peak[", rownames(data), "]" ))
rownames(samplesATAC) <- c("SampleID", "CellType")

dataATAC <- dataATAC[, c(-1,-2,-3)]
samplesATAC <- samplesATAC[, !grepl("pHSC", samplesATAC[2, ])]
data1ATAC <- dataATAC[ , c(colnames(dataATAC) %in% samplesATAC[1, ])]
data1ATAC <- apply(data1ATAC, 2, as.numeric)

unmapes_data_cells <- umap(t(apply(data1ATAC, 2, as.numeric)), n_neighbors = 5L, cores=4, knn = TRUE)
clustered <- kmeans(unmapes_data_cells$layout, centers=15)
clusters <- data.frame(clusters=clustered$cluster)

result <- data.frame(x=unmapes_data_cells[["layout"]][,1], y=unmapes_data_cells[["layout"]][,2])
clusters1 <- clusters[match(rownames(result), rownames(clusters)),] 
clusters$SampleID <- rownames(clusters)

result$SampleID <- rownames(result)
result <- merge(result, clusters, by="SampleID")
result <- merge(result, samples1, by="SampleID")

ggplot(data=result, aes(x=x, y=y, color=as.factor(CellType)))+geom_point()+labs(title = "Clustering for ATAC-seq", x = "First component", y = "Second component", color="CellType" )
ggsave(plot=last_plot(), device="pdf", filename="ATAC_umap_Cell_types.pdf", width=7, height=5)

samplesATAC <- samplesATAC[ , grepl("HSC", samplesATAC[ 2, ]) | grepl("CD4", samplesATAC[2, ])]

data1ATAC <- dataATAC[ , c(colnames(dataATAC) %in% samplesATAC[1, ])]
data1ATAC <- apply(data1ATAC, 2, as.numeric)

counts <- c()
len <- nrow(data1ATAC)
for(i in 1:len){
	counts <- c(counts, sum(data1ATAC[i, ] > 5 ))
	}

featuresATAC$num_cells_expr <- counts
featATAC <- AnnotatedDataFrame(data=featuresATAC)

df1 <- data.frame(x=colnames(data1ATAC))
samples1ATAC <- samplesATAC[ , match(df1$x, samplesATAC[ 1, ])] 
rownames(samples1ATAC) <- c("SampleID", "CellType")
samples1ATAC <- data.frame(t(samples1ATAC))
sampATAC <-AnnotatedDataFrame(data=samples1ATAC)


ATAC <- newCellDataSet(as.matrix(data1ATAC), phenoData = sampATAC, featureData = featATAC, expressionFamily=negbinomial.size())
ATAC <- estimateSizeFactors(ATAC)
ATAC <- estimateDispersions(ATAC)

disp_table <- dispersionTable(ATAC)
expressed_genesATAC <- row.names(subset(fData(ATAC),
    num_cells_expr >= 4))
    
diff_test_resATAC <- differentialGeneTest(ATAC[expressed_genesATAC,], fullModelFormulaStr = "~CellType", cores=5)
ATAC_ordering_genes <-
    row.names(diff_test_resATAC)[order(diff_test_resATAC$qval)][1:1000]
    
ATAC <- setOrderingFilter(ATAC,
        ordering_genes = ATAC_ordering_genes)
ATAC <- reduceDimension(ATAC, method = 'DDRTree')
ATAC <- orderCells(ATAC)
plot_cell_trajectory(ATAC, color_by="CellType")

# SAVE PLOT(OPTIONAL) 

diff_test_res_ptATAC <- differentialGeneTest(ATAC[ATAC_ordering_genes, ],
					   fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_namesATAC <- row.names(diff_test_res_ptATAC)[order(diff_test_res_ptATAC$qval)][1:100]

plot_pseudotime_heatmap(ATAC[sig_gene_namesATAC,],
                num_clusters = 2,
                cores = 1,
                show_rownames = T)
                
                
#SAVE PLOT
save(ATAC, file="ATACseq.Rdata")

#FOR DIFFERENTIAL GENE EXPRESSIONS
dataRNA <- read.table("GSE74246_RNAseq_All_Counts.txt")
samplesRNA <- read.table("sampleInfoRNA.txt")

rownames(samplesRNA) <- samplesRNA[,1]
samplesRNA <- samplesRNA[ ,-1]

colnames(samplesRNA) <- samplesRNA[1, ]
colnames(dataRNA) <- dataRNA[1, ]
dataRNA <- dataRNA[-1, ]

rownames(dataRNA) <- dataRNA[ , 1] 
featuresRNA <- data.frame(gene_short_name=dataRNA[, 1])
rownames(samplesRNA) <- c("SampleID", "CellType")

dataRNA <- dataRNA[, -1]
samplesRNA <- samplesRNA[, !grepl("pHSC", samplesRNA[2, ])]
data1RNA <- dataRNA[ , c(colnames(dataRNA) %in% samplesRNA[1, ])]
data1RNA <- apply(data1RNA, 2, as.numeric)

unmapes_data_cells <- umap(t(data1RNA), n_neighbors = 5L, cores=4, knn = TRUE)
clustered <- kmeans(unmapes_data_cells$layout, centers=15)
clusters <- data.frame(clusters=clustered$cluster)
result <- data.frame(x=unmapes_data_cells[["layout"]][,1], y=unmapes_data_cells[["layout"]][,2])
clusters1 <- clusters[match(rownames(result), rownames(clusters)),] 
clusters$SampleID <- rownames(clusters)



result$SampleID <- rownames(result)
result <- merge(result, clusters, by="SampleID")
samples1RNA <- data.frame(t(samplesRNA))
result <- merge(result,samples1RNA , by="SampleID")

ggplot(data=result, aes(x=x, y=y, color=as.factor(CellType)))+geom_point()+labs(title = "Clustering for RNA-seq", x = "First component", y = "Second component", color="CellType" )
ggsave(plot=last_plot(), device="pdf", filename="RNA_umap_Cell_types.pdf", width=7, height=5)



samplesRNA <- samplesRNA[ , grepl("HSC", samplesRNA[ 2, ]) | grepl("CD4", samplesRNA[2, ])]

data1RNA <- data1RNA[ , c(colnames(data1RNA) %in% samplesRNA[1, ])]
data1RNA <- apply(data1RNA, 2, as.numeric)

counts <- c()
len <- nrow(data1RNA)
for(i in 1:len){
  counts <- c(counts, sum(data1RNA[i, ] > 5 ))
}

featuresRNA$num_cells_expr <- counts

featRNA <- AnnotatedDataFrame(data=featuresRNA)

df1 <- data.frame(x=colnames(data1RNA))
samples1RNA <- samplesRNA[ , match(df1$x, samplesRNA[ 1, ])] 
rownames(samples1RNA) <- c("SampleID", "CellType")
samples1RNA <- data.frame(t(samples1RNA))
sampRNA <-AnnotatedDataFrame(data=samples1RNA)

RNA <- newCellDataSet(as.matrix(data1RNA), phenoData = sampRNA, featureData = featRNA, expressionFamily=negbinomial.size())
RNA <- estimateSizeFactors(RNA)
RNA <- estimateDispersions(RNA)

disp_table <- dispersionTable(RNA)
expressed_genesRNA <- row.names(subset(fData(RNA),
    num_cells_expr >= 4))
    
diff_test_resRNA <- differentialGeneTest(RNA[expressed_genesRNA,], fullModelFormulaStr = "~CellType", cores=5)
RNA_ordering_genes <-
    row.names(diff_test_resRNA)[order(diff_test_resRNA$qval)][1:1000]
    
RNA <- setOrderingFilter(RNA,
        ordering_genes = RNA_ordering_genes)
RNA <- reduceDimension(RNA, method = 'DDRTree')
RNA <- orderCells(RNA)
plot_cell_trajectory(RNA, color_by="CellType")
#SAVE PLOT(OPTIONAL)



diff_test_res_ptRNA <- differentialGeneTest(RNA[RNA_ordering_genes, ],
                                            fullModelFormulaStr = "~sm.ns(Pseudotime)", cores=5)

sig_gene_namesRNA <- row.names(diff_test_res_ptRNA)[order(diff_test_res_ptRNA$qval)][1:100]


plot_pseudotime_heatmap(RNA[sig_gene_namesRNA, ],
                        num_clusters = 2,
                        cores = 1,
                        show_rownames = T)






