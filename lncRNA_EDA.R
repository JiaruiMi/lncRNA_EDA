# RNA-seq analysis of zebrafish pancreatic cells from Bernard Peer's lab

#=======================================================================================
#
#                     DESeq object creation and Normalization
#
#=======================================================================================

## Set working directory
setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/lncRNA_EDA')
### getwd()
### dir()


# Load packages
library('reshape2')
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("amap") ## 'hcluster' function
library("ggplot2")
library("ggcorrplot")
library("biomaRt")
library('pheatmap')
library("ggbiplot")
library("ggrepel")
library("devtools")

## packages for enrichment analysis
library("org.Dr.eg.db")
readable = T
library("DOSE")
library("clusterProfiler")
library("enrichplot")
library('viridis')
library("BiocParallel")
library("stringr")
# library("pathview")
library("topGO")
library("ReactomePA")


#####################################################################################

# Read in data (counts/reads table and sampleGroup)
## Read in counts/reads table
reads <- read.csv('counts_total.csv', header = T, row.names = 1)
reads[1:3,1:3];dim(reads)
## In total, 23 samples with 32266 genes detected

## Initial Gene filtering
reads <- reads[rowSums(reads)>2,];dim(reads)

## Read in Coldata(Phenodata)
sample <- read.csv('SampleGroup_total.csv', header = T, row.names = 1, colClasses = 'factor')
head(sample);dim(sample);levels(sample$celltype);table(sample$celltype)

## For visualization of normalization (ggplot2 tidy data preparation)
group_List <- sample$celltype
exprSet_L <- row.names(reads)
exprSet_L <- cbind(exprSet_L, reads); head(exprSet_L)
exprSet_L <- melt(data = exprSet_L, id.vars = 'exprSet_L'); head(exprSet_L)
exprSet_L$group <- rep(group_List, each = nrow(reads)); head(exprSet_L)
## Sequencing depth investigation ##
# boxplot
p <- ggplot(data = exprSet_L, aes(x = variable, y = value, fill = group))+ geom_boxplot()
p
# violinplot
p <- ggplot(data = exprSet_L,aes(x = variable, y = value, fill = group))+geom_violin()
p

############### Make DESeq object and Do the normalization (size factor) #####################

# Make DESeq object and Do the normalization (size factor)
## Create DESeqDataSet
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = reads, colData = sample, design  = ~ celltype)

## Calculation and Normalization: get normalized_counts
dds <- DESeq(ddsFullCountTable)
normalized_counts <- counts(dds, normalized= T)
class(normalized_counts); head(normalized_counts)

## After Size factor normalization, check the sequencing depth across samples
exprSet_L <- row.names(normalized_counts)    
exprSet_L <- cbind(exprSet_L, as.data.frame(normalized_counts)) # Must convert normalized_count into data.frame before cbind
exprSet_L <- melt(data = exprSet_L, id.vars = 'exprSet_L')
exprSet_L$group <- rep(group_List, each = nrow(normalized_counts))
## Sequencing depth investigation ##
### boxplot
p <- ggplot(data = exprSet_L, aes(x = variable, y = value, fill = group))+ geom_boxplot()
print(p)
#### violinplot
p <- ggplot(data = exprSet_L,aes(x = variable, y = value, fill = group))+geom_violin()
p

#==============================================================================================
#
#     Scater object creation, Gene filtering and Exploratory data analysis (SC3, M3Drop)
#                     Get the genelist of highly variable genes
#
#==============================================================================================

# Load packages
library('SingleCellExperiment')
library('scater')
library('SC3')
options(stringsAsFactors = F)
data <- SingleCellExperiment(assays = list(counts = as.matrix(reads)), colData = sample)
data
data <- calculateQCMetrics(data)


## We don't do cell filtering here, but we do need to check the gene expression (Highly Expressed Genes)
plotQC(data, type = 'highest-expression', n = 20)
# ENSDARG00000099970(MALAT1), ENSDARG00000035350(ins), ENSDARG00000042993(try), ENSDARG00000073742(prss59.2),
# ENSDARG00000014190(sst2), ENSDARG00000079274(prss59.1), ENSDARG00000021339(cpa5), ENSDARG00000079296(gcga),
# ENSDARG00000017314(zgc:92041, elastase),ENSDARG00000080337(mtRNA), ENSDARG00000033161(sst1.2), ENSDARG00000007276(ela3l)
# ENSDARG00000045442(cpb1), ENSDARG00000090428(ctrb1), ENSDARG00000039730(zgc:112160, hydrolase)
# ENSDARG00000074547(si:ch211-240l19.8, protein coding), ENSDARG00000056744(ela2), ENSDARG00000029822(cel.2)
# ENSDARG00000056765(ela2l), ENSDARG00000063905(mt-co1)


## Data visualization -- PCA and tSNE (After sample QC)
### Before normalization (read counts)
plotPCA(data, exprs_values = 'counts', colour_by = 'celltype', size_by = 'total_features')
plotPCA(data, exprs_values = 'counts', colour_by = 'celltype', size_by = 'total_features', ntop = 500)
plotTSNE(data, exprs_values = 'counts', perplexity = 3, check_duplicates = F, 
         colour_by = 'celltype', size_by = 'total_features', 
         rand_seed = 123456)
plotTSNE(data, exprs_values = 'counts', perplexity = 3, check_duplicates = F, ntop = 500,
         colour_by = 'celltype', size_by = 'total_features', 
         rand_seed = 123456)

### Normalization (size factor/RLE)
data  # Before Normalization, the 'data' object only contains 'counts' expression dataset
data <- normaliseExprs(data, method = "RLE", return_log = T, return_norm_as_exprs = T)
data # After Normalization, the 'data' object contains three expression dataset: counts, normcounts, logcounts
plotPCA(data, exprs_values = 'normcounts', colour_by = 'celltype', size_by = 'total_features', ntop = 500)
plotPCA(data, exprs_values = 'logcounts', colour_by = 'celltype', size_by = 'total_features', ntop = 500)
plotTSNE(data, exprs_values = 'normcounts', perplexity = 3, check_duplicates = F, ntop = 500,
         colour_by = 'celltype', size_by = 'total_features', 
         rand_seed = 123456)
plotTSNE(data, exprs_values = 'logcounts', perplexity = 3, check_duplicates = F, ntop = 500,
         colour_by = 'celltype', size_by = 'total_features', 
         rand_seed = 123456)
plotRLE(data, exprs_mats = list(Raw_read_counts = 'counts', Normalized_read_counts = 'normcounts', log_normalized_counts = 'logcounts'),
        exprs_logged = c(T,T,T), colour_by = 'celltype')



## Biological Analysis -- Clustering 
### load packages
library('SC3')

### Use SC3 to estimate k value (the number of clusters)
data <- sc3_estimate_k(data)
data 
metadata(data)$sc3$k_estimation
#### Before running sc3,you should add feature_symbol to the rowData slot which contains the rownames of the genes
rowData(data)$feature_symbol <- row.names(counts(data)) 
#### Running SC3 (we also ask it to calculate biological properties of the clusters)
data <- sc3(data, ks = 5, biology = T)

#### Consensus matrix
sc3_plot_consensus(data, k = 5, show_pdata = 'celltype')
#### Silhouette plot
sc3_plot_silhouette(data, k = 5)
#### Heatmap of expression matrix
sc3_plot_expression(data, k = 5, show_pdata = 'celltype')
#### Identify marker genes
# sc3_plot_markers(data, k = 5, show_pdata = 'celltype') ops! Have problem here, maybe due to the genesymbol




## Biological Analysis -- Feature Selection
### In most situations only a portion of those will show a response to the biological conditions of interest. Most genes
### detected will only be detected at different levels due to technical noise. One consequence of this is that technical
### noise and batch effect can obscure the biological signal of interest. Thus, it is often advantageous to perform feature
### selection to remove those genes which noly exhibit technical noise from downstream analysis. Not only does this 
### generally increase the signal:noise ratio in the data; it also reduces the computational complexity of analyses, by
### reducing the total amount of data to be processed.
### Load packages
library('M3Drop')
library('matrixStats')
### Highly variable genes: The method proposed to identify features was to identify highly variable genes (HVG). HVG assumes
### that if genes have large differences in expression across cells, some of those differences are due to biological
### differences between the cells rather than technical noise. However, because of the nature of count data, there is 
### a positive relationship between the mean expression of a gene and the variance in the read counts across cells.
### This relationship must be corrected for to properly identify HVGs.
plot(rowMeans(counts(data)), rowVars(counts(data)), xlab = 'Mean Expression', ylab = 'Variance')
plot(rowMeans(counts(data)), rowVars(counts(data)), log = 'xy', xlab = 'Mean Expression', ylab = 'Variance')
### A popular method to correct for the relationship between variance and mean expression was proposed by Brennecke.
### To use the Brennecke method, we first normlize for library size then calculate the mean and square coeffcient of 
### variation (variation divided by the squared mean expression). A quadratic curve is fit to the relationship between 
### these two variables, and then a chi-square test is used to find genes significantly above the curve. This method is
### included in the M3Drop package as the BrenneckeGetVariableGenes function. We will use the entire dataset to estimate
### the technical noise.

#### In the figure, below the red curve is the fitted technical noise model and the dashed line is the 95% CI. Pink dots
#### are the genes with significant biological variability after multiple-testing correction.
Brennecke_HVG <- BrenneckeGetVariableGenes(logcounts(data), fdr = 0.01, minBiolDisp = 0.5)
str(Brennecke_HVG)
class(Brennecke_HVG)

#### Using these Highly Variable Genes, we can do clustering and PCA plot
##### Clustering with normalized_counts derived from DESeq2
normalized_counts_HVG <- normalized_counts[rownames(normalized_counts) %in% Brennecke_HVG,]
pearson_cor <- as.matrix(cor(normalized_counts_HVG, method = 'pearson'))
head(pearson_cor)
hc <- hcluster(t(normalized_counts_HVG), method="pearson")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(pearson_cor, Rowv = as.dendrogram(hc), trace = 'none',symm = T, col = hmcol, main = 'The pearson correlation of each')
pheatmap(pearson_cor)

##### Clusterig with Log-transformed normalized counts derived from DESeq2
rld <- rlog(dds, blind = F) # This step takes some time
rlogMat <- assay(rld)
rlogMat_HVG <- rlogMat[rownames(rlogMat) %in% Brennecke_HVG,]
pearson_cor <- as.matrix(cor(rlogMat_HVG, method = 'pearson'))
hc <- hcluster(t(rlogMat_HVG), method="pearson")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(pearson_cor, Rowv = as.dendrogram(hc), trace = 'none',symm = T, col = hmcol, main = 'The pearson correlation of each')
pheatmap(pearson_cor)





#==============================================================================================
#
#          Seurat object creation, Gene filtering and Exploratory data analysis
#                    Get the genelist of highly variable genes
#
#==============================================================================================

## Loading packages and Seurat object generation
library(Seurat)
library(mclust)
library(dplyr)

### Seurat introduces its own object class -- seurat. All calculations in this chapter are performed
### on an object of this class. To begin the analysis we first need to initialize the object with 
### the raw (non-normalized) data. We can determine the threshold for gene expression (e.g. >= 3 cells)
### with at least (e.g. 200) detected genes.

#### The Seurat object generation step is a little tricky, make sure you load the package. You can checked
#### packages to see whether seurat has been ticked.

#### The raw data slot represents the original expression matrix, input when creating the Seurat object,
#### and prior to any preprocessing by Seurat. This could represent the UMI matrix generated by DropSeqTools,
#### or 10*CellRanger, a count matrix from featureCounts, an FPKM matrix produced by Cufflinks, or a TPM
#### matrix produced by RSEM. Either raw counts or normalized values are fine, but the input expression
#### should not be log-transformed.
seuset <- CreateSeuratObject(raw.data = as.matrix(normalized_counts)) 
seuset; str(seuset)
seuset@raw.data
#### We can see in the meta.data, there are nGene and nUMI. Seurat will detect the ensembl_gene_id 
#### automatically and convert them into entrezgene. The nUMI here is the number of reads in total.
#### "The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat,
#### For non-UMI data, nUMI represents the sum of the non-normalized/normalized values within a cell.
VlnPlot(object = seuset, features.plot = c('nGene','nUMI'), nCol = 2)
GenePlot(object = seuset, gene1 = 'nUMI', gene2 = 'nGene')

### Normalization
#### By default, a global-scaling normalization method 'logNormalize" that normalizes the gene expression
#### measurements for each cell by the total expression, multiplies this by a scale factor (10000 by default)
#### and log transformed the result. It seems like log-transformed CPM. More methods to be added shortly.
seuset <- NormalizeData(object = seuset, normalization.method = 'LogNormalize', scale.factor = 10000)
seuset@data
#### Check the raw.data and data slot, they are different. The data slot (object@data) stores normalized
#### and log-transformed single cell expression. This maintains the relative abundance levels of all genes,
#### and contains only zero or positive values. This data is used for visualization, such as violin and
#### feature plots, most differential expression tests, finding high variable genes, and as input to 
#### ScaleData.


### Highly Variable Genes
#### Seurat calculates highly variable genes and focuses on these for downstream analysis. FindVariableGenes
#### calculates the average expression and dispersion for each gene, places these genes into bins, and then
#### calculates a z-score for dispersion within each bin. This helps control for the relationship between
#### variability and average expression. The object is stored in object@var.genes. This function is unchanged
#### from (Macosko et al.), but new methods for variable gene expression identification are coming soon. We 
#### suggest that users set these parameters to mark visual outliers on the dispersion plot, but the exact
#### parameter settings may vary based on the data type, heterogeneity in the sample, and normalization
#### strategy. This step takes some time. In the plot, the highly variable genes would be labeled with
#### ensembl_gene_id
seuset <- FindVariableGenes(object = seuset, mean.function = ExpMean, dispersion.function = LogVMR,
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(seuset@var.genes)

### Scaling the data and removing unwanted sources of variation
#### Your single cell dataset likely contains 'uninteresting' sources of variation. This could include not
#### only technical noise, but batch effects, or even biological sources of variation (cell cycle stage), As
#### suggested in Buettner et al, NBT, 2015, regressing these signals out of the analysis can improve downstream
#### dimensionality reduction and clustering. To mitigate the effect of these signals, Seurat constructs linear
#### models to predict gene expression based on user-defined variables. The scaled z-scored residuals of these 
#### models are stored in the scale.data slot, and are used for dimensionality reduction and clustering.
seuset <- ScaleData(object = seuset)


### Linear dimensionality reduction (PCA)
#### Seurat will perform on the ScaleData slot , what we need to do is
#### to input the variable genes and run PCA. By default, the genes in object@var.genes are used as input,
#### but can be defined using pc.genes. We have typically found that running dimensionality reduction on
#### highly variable genes can improve performance.
seuset <- RunPCA(object = seuset, pc.genes = seuset@var.genes,
                 do.print = T, pcs.print = 1:5,
                 genes.print = 5)
#### Seurat provides several useful ways of visualizing both cells and genes that define the PCA, including
#### PrintPCA, VizPCA, PCAPlot and PCHeatmap
PrintPCA(object = seuset, pcs.print = 1:5, genes.print = 5, use.full = F)
VizPCA(object = seuset, pcs.use = 1:2)
PCAPlot(object = seuset, dim.1 = 1, dim.2 =2)
#### In particular, PCHeatmap allows for easy exploration of the primary sources of heterogeneity in a dataset,
#### and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells
#### and genes are ordered according to their PCA scores. Setting cell.use to a number plots the extreme cells
#### on both ends of the spectrum, which dramatically speeds plotting for large datasets. Though clearly a 
#### supervised analysis, we find this to be a valuable tool for exploring correlated gene sets:
PCHeatmap(object = seuset, pc.use = 1, cells.use = 23, do.balanced = T, 
          label.columns = F, use.full = F)
PCHeatmap(object = seuset, pc.use = 1:6, cells.use = 23, do.balanced = T,
          label.columns = F, use.full = F)


#### To overcome the extensive technical noise in any single gene, Seurat clusters cells based on their PCA
#### score, with each PC essentially representing a 'metagene' that combines information across a correlated
#### gene set. Determining how many PCs to include downstream is therefore an important step. In Macosko et al,
#### we implemented a resampling rest inspired by the jackStraw procedure. We randomly permute a subset of the
#### data (1% by default) and rerun PCA, constructing a 'null distribution' of gene scores, and repeat this 
#### procedure. We identify 'significant' PCs as those who have a strong enrichment of low p-value genes.
seuset <- JackStraw(object = seuset, num.replicate = 100, do.print = F)
#### Here, I got a lot of warning message!

#### The JackStrawPlot function provides a visualization tool for comparing the distribution of p-value for
#### each PC with a uniform distribution (dashed line). 'Significant' PCs will show a strong enrichment of 
#### genes with low p-values (solid curve above the dashed line); In this case it appears that PCs 1-5 
#### are significant.
JackStrawPlot(object = seuset, PCs = 1:9)

#### A more ad hoc method for determining which PCs to use is to look at a plot of the standard deviation of 
#### the principal components and draw your cutoff where there is a clear elbow in the graph. This can be done
#### with PCElbowPlot.
PCElbowPlot(object = seuset)



### Cluster the cell
#### Seurat now includes an graph-based clustering approach compared to Macosko et al. Importantly, the
#### distance metric which drives the clustering analysis (based on previously identified PCs) remains the
#### same. However, our approach to partioning the cellular distance matrix into clusters has dramatically 
#### improved. Briefly, these methods embed cells in a graph structure -- for example a K-nearest neighbour
#### (KNN) graph, with edges drawn between cells with similar gene expression pattern, and then attempt to
#### partition this graph into highly interconnected 'quasi-cliques' or 'communities'. As in PhenoGraph, we 
#### first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights 
#### between any two cells, we apply modularity optimization techniques, to iteratively group cells together,
#### with the goal of optimizing the standard modularity function.
#### RunPCA and ScaleData shoule be performed in prior!!!!!
seuset <- FindClusters(object = seuset, reduction.type = 'pca', dims.use = 1:6,
                       resolution = 1.0, print.output = 0, save.SNN = T, k.param = 4)
PrintFindClustersParams(object = seuset)
table(seuset@ident)
#### Here we can see there are five clusters, with 6,5,5,4,3 samples respectively.


### Run non-linear dimensional reduction (tSNE)
#### Seurat continues to use tSNE as a powerful tool to visualize and explore these datasets, while we no 
#### longer advise clustering directly on tSNE component, cells within the graph-based clusters determied
#### above should co-localize on the tSNE plot This is because the tSNE aims to place cells with similar
#### local neighborhoods in high-dimensional space together in low-dimensional space. As input to the tSNE,
#### we suggest using the same PCs as input to the clustering analysis, although computing the tSNE based 
#### on scaled gene expression is also supported using the genes.use argument.
seuset <- RunTSNE(object = seuset, dims.use = 1:6, do.fast = T, perplexity = 4, check_duplicates = F)
TSNEPlot(object = seuset, do.label = T, pt.size = 5)
#### The tSNE clusteing looks good, but there is no label on it.



### Find differentially expressed genes (cluster biomarkers)
#### Seurat can help you find markers that define clusters via differential expression. By default, it 
#### identifies positive and negative markers of a single cluster (specified in inden.1), compared to all 
#### other cells. FindAllMarkers automates this process for all clusters, but you can also test groups of 
#### clusters vs each other, or against all cells.

#### Find all markers of cluster 0,1,2,3,4; we can check avg_logFC to see whether the gene
#### is highly expressed (> 0) or lowly expressed (< 0) in one cluster of cell compared with other
#### cell clusters.
markers0 <- FindMarkers(seuset,ident.1 = 0) ## cluster0 is "acinal cell"
markers1 <- FindMarkers(seuset, 1) ## cluster1 is "alpha cell"
markers2 <- FindMarkers(seuset, 2) ## cluster2 is "beta cell"
markers3 <- FindMarkers(seuset, 3) ## cluster3 is "delta cell"
markers4 <- FindMarkers(seuset, 4) ## cluster4 is "ductal cell
print(x = head(x = markers0, n = 10)) ## print cluster0 marker genes with p_val, avg_logFC, pct.1, pct.2, p_val_adj
print(x = head(x = markers1, n = 10)) ## print cluster1 marker genes with p_val, avg_logFC, pct.1, pct.2, p_val_adj
print(x = head(x = markers2, n = 10)) ## print cluster2 marker genes with p_val, avg_logFC, pct.1, pct.2, p_val_adj
print(x = head(x = markers3, n = 10)) ## print cluster3 marker genes with p_val, avg_logFC, pct.1, pct.2, p_val_adj
print(x = head(x = markers4, n = 10)) ## print cluster4 marker genes with p_val, avg_logFC, pct.1, pct.2, p_val_adj
VlnPlot(object = seuset, features.plot = rownames(markers0)[1:6]) 
VlnPlot(object = seuset, features.plot = rownames(markers1)[1:6])
VlnPlot(object = seuset, features.plot = rownames(markers2)[1:6])
VlnPlot(object = seuset, features.plot = rownames(markers3)[1:6])
VlnPlot(object = seuset, features.plot = rownames(markers4)[1:6])


#### Find all markers distinguishing cluster 2 from cluster 1 and 3 (markers in beta cell distinguishing from alpha
#### cell and delta cell)
#### Adding speed by excluding tests:
#### min.pct - controls for sparsity
#### min percentage in a group
#### thresh.test - must have this difference in average
#### Many choices for DE: bimod(tests differences in mean and proportions), roc(uses AUC like definition),
#### t(student t test), tobit(tobit regression on a smoothed data)
cluster2.markers <- FindMarkers(object = seuset, ident.1 = 2, ident.2 = c(1,3), min.pct = 0.25)
print(x = head(x = cluster2.markers, n = 10))
#### If you want to check the marker genes, you can type in e.g. rownames(markers0)
#### You can do enrichment analysis of the marker genes for each cluster

#### Find all markers for every cluster compared to all remaining cells, report only the positive ones; FindAllMarkers
#### automates this process and find markers for all clusters:
marker <- FindAllMarkers(object = seuset, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
top10 <- marker %>% group_by(cluster) %>% top_n(10, avg_logFC)

#### We include several tools for visualizing marker expression, VlnPlot (shows expression probability
#### distributions across clusters), and FeaturePlot (visualizes gene expression on a tSNE or PCA plot) are
#### our most commonly used visualizations. We also suggest exploring RidgePlot, CellPlot, and DotPlot as
#### additional methods to view your datasets.
FeaturePlot(seuset, features.plot = head(rownames(markers2), cols.use = c("gray", "blue"), nCol = 3))

#### Do Heatmap generates an expression heatmap for given cells and genes. In this case, we are plotting the
#### top 10 markers (or all markers if less than 10) for each cluster.
DoHeatmap(object = seuset, genes.use = top10$gene, slim.col.label = T, remove.key = T)



### Assigning cell type identity to clusters
#### Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased
#### clustering to known cell type:
current.cluster.ids <- c(0,1,2,3,4)
new.cluster.ids <- c('acinal cells','alpha-cells', 'beta-cells', 'delta-cells','ductal cells')
seuset@ident <- plyr::mapvalues(x = seuset@ident, from = current.cluster.ids, to = new.cluster.ids)
#### RunTSNE should be performed in prior!!!!!
TSNEPlot(object = seuset, do.label = T, pt.size = 5)



#==================================================================================
#
#        Regular DESeq2 Differential Gene Expression Analysis pipeline
#
#==================================================================================

############# Continue performing normalizaton #################
## Sort according to mad
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts < normalized_counts[order(normalized_counts_mad, decreasing = T),]
## Log transformation
rld <- rlog(dds, blind = F) # This step takes some time, we do not put '+1' here
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing = T),]

### check one gene expression profiles across different samples types, take yap1 "ENSDARG00000068401" as an example
singleGeneExprs <- as.vector(as.matrix(rlogMat[row.names(rlogMat) == "ENSDARG00000036456",]))
singleGeneSample <- as.vector(as.matrix(sample[,1]))
singleGeneTable <- data.frame(exprs = singleGeneExprs, cellType = singleGeneSample)
p <- ggplot(data = singleGeneTable, aes(x = cellType, y = exprs, fill = cellType)) + 
     geom_boxplot() + geom_point() + theme_bw(base_size = 10) +  
     labs(title = "anxa4 expression across different cell types", x = 'cell type', y = 'Mean of Normalized counts (log2-transformed)') + 
     theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16)) + 
     theme(plot.title = element_text(hjust = 0.5, size = 22))
p

### After Size factor normalization and log transformation, we can check the expression matrix again and compare the 
### sequencing depth among different samples
exprSet_L <- row.names(rlogMat)    
exprSet_L <- cbind(exprSet_L, as.data.frame(rlogMat)) # Must convert normalized_count into data.frame before cbind
exprSet_L <- melt(data = exprSet_L, id.vars = 'exprSet_L')
exprSet_L$group <- rep(group_List, each = nrow(rlogMat))


############# Sequencing depth investigation -- After Normalization and log transformation #############
### Load the package
library(ggplot2)
#### boxplot
p <- ggplot(data = exprSet_L, aes(x = variable, y = value, fill = group))+ geom_boxplot()
print(p)
p <- p +stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
p
p <- p + theme_set(theme_set(theme_bw(base_size=20)))
p
p <- p + theme(text=element_text(face='bold'),axis.text.x=element_text(angle=90,hjust=1),axis.title=element_blank())
p
#### violinplot
p <- ggplot(data = exprSet_L,aes(x = variable, y = value, fill = group))+geom_violin()+
  stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")+
  theme_set(theme_set(theme_bw(base_size=20)))+
  theme(text=element_text(face='bold'),axis.text.x=element_text(angle=90,hjust=1),axis.title=element_blank())
p
#### histogram (Here we can see the negative binomial distribution of read counts of all genes)
p <- ggplot(exprSet_L, aes(value, fill = group))+geom_histogram(bins = 200)+facet_wrap(~variable, nrow = 4)
p
#### density plot
p <- ggplot(exprSet_L, aes(value, fill = group, col = group))+geom_density()+facet_wrap(~variable, nrow = 4)
p
p <- ggplot(exprSet_L, aes(value, col = group))+geom_density()
p


############# Hierarchical Clustering and PCA ###############

## Hierarchical clustering
### Selecting colors
display.brewer.all(type = 'all') 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
### Calculating the pearson correlation matrix
#### Clustering is based on the relationship among each individual sample. Both distance and correlation
#### are accepted but correlation is preferred. The distance methods and correlation methods are chosen 
#### based on your data.
pearson_cor <- as.matrix(cor(rlogMat, method = 'pearson'))
head(pearson_cor); dim(pearson_cor)
hc <- hcluster(t(rlogMat), method="pearson")
#### Here I use 'heatmap.2' function from 'gplots' package.
heatmap.2(pearson_cor, Rowv = as.dendrogram(hc), trace = 'none',symm = T, col = hmcol, main = 'The pearson correlation of each')
pheatmap(pearson_cor)


## Principal Component Analysis -- PCA
### How to examine the loading scores to determine what variables have the largest effect on the graph
### Prcomp will give us three things: x, sdev, rotation
### x contains the prcincipal component for drawing a graph
### The prcomp() function calls the loading scores rotation, there are loading scores for each PC
pca <- prcomp(t(rlogMat), scale. = T)   # remember to transpose the matrix
plot(pca$x[,1],pca$x[,2])
ggbiplot(pca ,choices = 1:2, obs.scale = T,labels = NULL, var.scale = T,groups = sample$celltype, ellipse = T, circle = T, var.axes = F, alpha = 0.5)+
  theme(legend.direction = 'horizontal', legend.position = 'top')+theme_bw()
### loadings for PC1, sorting based on the contribution from large to small
pca$rotation; dim(pca$rotation)  ## have a look at the loading scores of each principal component
loading_score1 <- pca$rotation[,1]
gene_score1 <- abs(loading_score1)
gene_score_ranked1 <- sort(gene_score1, decreasing = T)
top_25_genes1 <- names(gene_score_ranked1[1:25])
pca$rotation[top_25_genes1,1]
### loadings for PC2, sorting based on the contribution from large to small
loading_score2 <- pca$rotation[,2]
gene_score2 <- abs(loading_score2)
gene_score_ranked2 <- sort(gene_score2, decreasing = T)
top_25_genes2 <- names(gene_score_ranked2[1:25])
top_25_genes2
pca$rotation[top_25_genes2,2]


### Using plotPCA from DESeq2 package
#### Use top 5000 genes (selected variable genes to do PCA is better)
class(rld)
pca_data <- plotPCA(rld, intgroup = c('celltype'), returnData = T, ntop = 5000)    # loading for each principal component (compared with novel analysis with the same important loadings with novel transcripts)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(PC1, PC2, color=celltype))
p + geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))
#### Without gene selection step, the PCA plot seems not be able to distinguish among endocrine cells
plotPCA(rld, intgroup=c("celltype"))
ggcorrplot(pearson_cor, hc.order = T,lab = T)     



############### Differential Gene Expression #################
sampleA <- 'beta'
sampleB <- 'acinal'
sampleC <- 'alpha'
sampleD <- 'delta'
sampleE <- 'ductal'

### Differential gene expression between sampleA and sampleB
#### Comparison between sampleA and sampleB
contrastV <- c('celltype', sampleA, sampleB)
res <- results(dds, contrast = contrastV)
head(res)

#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA))
colnames(baseMeanA) <- sampleA

#### Calculate the mean value of gene expression in sampleB 
baseB <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleB]
baseMeanB <- as.data.frame(rowMeans(baseB))
colnames(baseMeanB) <- sampleB

#### Prepare the res table of differential gene expression between sampleA and sampleB
res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)
##### Please have a look at the 'res' object and to see the mean value of sampleA, sampleB 
##### and logfoldchange, here I can draw the conclusion that the logfoldchange is based on the comparison
##### sampleA vs sampleB



#### Draw MA plot
plotMA(res[,c(4,5,10)], main = "MAplot")

#### Save the comparison results
com <- paste(sampleA,'_vs_', sampleB)
file_base <- paste('DESeq2',com,sep = '.')
write.table(as.data.frame(res),file = file_base,quote = F, row.names = F)


#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleB,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


### Volcano plot
#### Just a volcano plot
logFC <- res$log2FoldChange
FDR <- res$padj
head(res)
plot(logFC,log10(FDR)*(-1),col = ifelse(FDR<=0.01, ifelse(abs(logFC)>=1, 'red', 'black'), 'black'), xlab = 'logFC', ylab = '-1*log10(FDR)', main = "Volcano plot", pch = '.', ylim=c(1,100))

#### Volcano plot showing upregulated and downregulated genes
res$change <- as.factor(ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 2, 
                               ifelse(res$log2FoldChange > 1, 'UP', 'DOWN'), 'NOT'))



##### Use biomaRt to convert gene ID of differential expressed genes with very significant p-value
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart)
listFilters(ensembl)
signifiant_genes <- getBM(attribute=c('ensembl_gene_id', 'entrezgene','zfin_id_symbol'),filters = 'ensembl_gene_id', 
                          values= row.names(subset(res, -log10(padj) > 100)),mart = ensembl)

##### Plot VolcanoPlot
this_title <- paste('VolcanoPlot of',sampleA, 'vs', sampleB,'\nCutoff for logFC is 1', 
                    '\nThe number of upregulated genes is', 
                    nrow(res[res$change == 'UP',]),'\nThe number of downregulated gene is', 
                    nrow(res[res$change == 'DOWN',]))

ggplot(data = res, aes(x=log2FoldChange, y = -log10(padj), color = change, alpha = 0.5)) +
  geom_point()+
  scale_color_manual(values = c('blue','black', 'red')) +
  geom_hline(yintercept = -log10(0.05),lty=4, lwd=0.6,alpha=0.5)+
  geom_vline(xintercept = c(-1,1),lty=4, lwd=0.6,alpha=0.5)+
  theme_bw(base_size = 10)+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title= this_title, x= 'log2(fold change)', y = '-log10(padj)')+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text_repel(data=subset(res, -log10(padj) > 100), aes(label= ID), col= 'black')



#### Heatmap of beta cell and acinal cell
res_de_up_top25_id <- as.vector(head(res_de_up$ID,25))
res_de_down_top25_id <- as.vector(head(res_de_down$ID,25))
res_de_top50 <- c(res_de_up_top25_id,res_de_down_top25_id)
res_de_top25_expr_up <- normalized_counts[rownames(normalized_counts) %in% res_de_up_top25_id,c(2:6,1,16:20)]
res_de_top50_expr_down <- normalized_counts[rownames(normalized_counts) %in% res_de_down_top25_id,c(2:6,1,16:20)]
res_de_top50_expr <- rbind(res_de_top25_expr_up, res_de_top50_expr_down)
##### Still need to do ID conversion
zfin_gene_symbol <- getBM(attribute=c('ensembl_gene_id', 'entrezgene','zfin_id_symbol'),
                          filters = 'ensembl_gene_id', values= row.names(res_de_top50_expr),mart = ensembl)
row.names(res_de_top50_expr) <- zfin_gene_symbol$zfin_id_symbol
annotation_row <- data.frame(GeneClass = factor(c(rep('beta-specific','acinal-specific', 25), rep('acinal-specific', 25))))
rownames(annotation_row) = row.names(res_de_top50_expr)
pheatmap(res_de_top50_expr, color = colorRampPalette(c('navy', 'white', 'firebrick3'))(50), 
         cluster_rows = T, scale = 'row', cluster_cols = T, annotation_col = sample, 
         annotation_row = annotation_row)




######### Enrichment analysis based on Differential gene expression analysis between two cell types #########
## Use biomaRt to concert gene ID (for enrichment analysis, entrezid should be used)
mart <- useMart('ensembl')
### select dataset "drerio_gene_ensembl"
ensembl <- useDataset('drerio_gene_ensembl', mart)
### beta_vs_acinal
entrezid_beta_vs_acinal<- getBM(attribute=c('ensembl_gene_id', 'entrezgene'),
                                filters = 'ensembl_gene_id', values= row.names(res_de),mart = ensembl)
entrezid_beta_vs_acinal <- entrezid_beta_vs_acinal[!is.na(entrezid_beta_vs_acinal$entrezgene),]
geneList_foldchange_beta_vs_acinal <- res_de[res_de$ID %in% entrezid_beta_vs_acinal$ensembl_gene_id,c('ID','log2FoldChange')]
head(geneList_foldchange_beta_vs_acinal)
#### The input for enrichment analysis is entrezgene id 
entrezid_beta_vs_acinal <- unique(entrezid_beta_vs_acinal$entrezgene)

#### Biological Process
BP <- enrichGO(entrezid_beta_vs_acinal,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
result_BP <- simplify(BP,cutoff = 0.7, by = 'p.adjust', select_fun = min)
dotplot(result_BP,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.04)+scale_color_continuous(low = 'purple', high = 'green')  # Use scale_size to change the size of the bubble, if the bubble is too large and be cut by the edge, you can used xlim to extend the x axis
dotplot(result_BP,showCategory = 10,x = 'count')+scale_size(range = c(2,12))+xlim(NA,120)
barplot(result_BP,drop = T, showCategory = 10)+scale_x_discrete(labels = function(x)str_wrap(x, width = 25))
barplot(result_BP,x = 'count', showCategory = 10)
plotGOgraph(result_BP)
goplot(result_BP)
#### Molecular Function
MF <- enrichGO(entrezid_beta_vs_acinal,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'MF', readable = readable)
result_MF <- simplify(MF,cutoff = 0.7, by = 'p.adjust', select_fun = min)
dotplot(result_BP,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.04)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_MF,drop = T, showCategory = 10)
plotGOgraph(result_MF)
goplot(result_MF)
#### Cellular Component
CC <- enrichGO(entrezid_beta_vs_acinal,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'CC', readable = readable)
result_CC <- simplify(CC,cutoff = 0.7, by = 'p.adjust', select_fun = min)
dotplot(result_CC,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.045)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_CC,drop = T, showCategory = 10)
plotGOgraph(result_CC)
goplot(result_CC)
#### KEGG pathway
kk <- enrichKEGG(entrezid_beta_vs_acinal,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)




# To get the pathway image: pv.out <- pathview(gene.data = geneList_foldchange_beta_higherthan_acinal$log2FoldChange,pathway.id = result_kk$ID, species = 'dre', kegg.native = T)
## Preparation of geneList for GSEA analysis
entrezid_DGE_beta_with_acinal<- getBM(attribute=c('ensembl_gene_id', 'entrezgene'),filters = 'ensembl_gene_id', values= row.names(res_de),mart = ensembl)
head(entrezid_DGE_beta_with_acinal)
dim(entrezid_DGE_beta_with_acinal)
dim(res_de)
entrezid_DGE_beta_with_acinal <- entrezid_DGE_beta_with_acinal[is.na(entrezid_DGE_beta_with_acinal$entrezgene) != T,]
summary(entrezid_DGE_beta_with_acinal)
head(entrezid_DGE_beta_with_acinal)
geneList_foldchange_beta_with_acinal<- res_de[res_de$ID %in% entrezid_DGE_beta_with_acinal$ensembl_gene_id,c('ID','log2FoldChange')]
head(geneList_foldchange_beta_with_acinal)
dim(geneList_foldchange_beta_with_acinal)
summary(entrezid_DGE_beta_with_acinal)
summary(geneList_foldchange_beta_with_acinal)
names(geneList_foldchange_beta_with_acinal)= c("ensembl_gene_id", "log2FoldChange")
geneList_foldchange_beta_with_acinal <- merge(entrezid_DGE_beta_with_acinal,geneList_foldchange_beta_with_acinal,by = "ensembl_gene_id")
summary(geneList_foldchange_beta_with_acinal)
head(geneList_foldchange_beta_with_acinal)
geneList_foldchange_beta_with_acinal <- geneList_foldchange_beta_with_acinal[,-1]
geneList_foldchange_beta_with_acinal_2 <- geneList_foldchange_beta_with_acinal[,2]
names(geneList_foldchange_beta_with_acinal_2)<- geneList_foldchange_beta_with_acinal[,1]
geneList_foldchange_beta_with_acinal_2 <- sort(geneList_foldchange_beta_with_acinal_2,decreasing = T)
head(geneList_foldchange_beta_with_acinal_2)
## gsea analysis of GO--biological process
gsea_analysis_GO_BP <- gseGO(geneList_foldchange_beta_with_acinal_2,ont = "BP", OrgDb = org.Dr.eg.db, verbose = F, pvalueCutoff = 0.1)
head(gsea_analysis_GO_BP)
dim(gsea_analysis_GO_BP)
gseaplot(gsea_analysis_GO_BP,geneSetID = "GO:0006508")
dotplot(gsea_analysis_GO_BP, showCategory = 30, split='.sign')+facet_grid(.~.sign)
ridgeplot(gsea_analysis_GO_BP, showCategory = 30,fill = 'pvalue')
## gsea analysis of GO--molecular function
gsea_analysis_GO_MF <- gseGO(geneList_foldchange_beta_with_acinal_2, ont = 'MF', OrgDb = org.Dr.eg.db, verbose = F, pvalueCutoff = 0.1)
head(gsea_analysis_GO_MF)
dim(gsea_analysis_GO_MF)
gseaplot(gsea_analysis_GO_MF,geneSetID = 'GO:0046914')
dotplot(gsea_analysis_GO_MF, showCategory = 30, split='.sign')+facet_grid(.~.sign)
ridgeplot(gsea_analysis_GO_MF, showCategory = 30, fill = 'pvalue')
## gsea analysis of GO--cellular component
gsea_analysis_GO_CC <- gseGO(geneList_foldchange_beta_with_acinal_2, ont = 'CC', OrgDb = org.Dr.eg.db, verbose = F, pvalueCutoff = 0.1)
head(gsea_analysis_GO_CC)
gseaplot(gsea_analysis_GO_CC, geneSetID = 'GO:0005576')
dotplot(gsea_analysis_GO_CC, showCategory = 30, split='.sign')+facet_grid(.~.sign)
ridgeplot(gsea_analysis_GO_CC, showCategory = 30, fill = 'pvalue')
## gsea analysis of KEGG
gsea_analysis_KEGG <- gseKEGG(geneList_foldchange_beta_with_acinal_2,organism = 'dre')
dotplot(gsea_analysis_KEGG)
ridgeplot(gsea_analysis_KEGG)
pv.out <- pathview(gene.data = geneList_foldchange_beta_with_acinal_2,pathway.id = result_kk$ID, species = 'dre', kegg.native = T)









## GENE ID conversion
### Use biomaRt to concert gene ID (For enrichment analysis, entrezid should be used; For text and visualization, zfin_id_symbol should be used)
mart <- useMart('ensembl')
# View(listDatasets(mart))  # select dataset "drerio_gene_ensembl"
ensembl <- useDataset('drerio_gene_ensembl', mart)
listFilters(ensembl)
reads_geneName<- getBM(attribute=c('ensembl_gene_id', 'entrezgene','zfin_id_symbol'),
                       filters = 'ensembl_gene_id', values= row.names(reads),mart = ensembl )
class(reads_geneName); dim(reads_geneName); View(reads_geneName)
reads_geneName[reads_geneName == "" | is.na(reads_geneName) | is.null(reads_geneName)] <- "NA"
?match
length(reads_geneName$zfin_id_symbol[match(reads_geneName$ensembl_gene_id, row.names(reads))])
sort(table(reads_geneName$ensembl_gene_id), decreasing = T)[1:100]       





#==============================================================================================
#
#          Performing enrichment analysis with the cluster marker genes from Seurat
#
#==============================================================================================
## Here we can take beta-cells as an example, the cluster number is 0; the genes (highly-expressed and lowly-expressed
## are stored in row.names(markers2))

class(markers2); dim(markers2); markers2[1:6,]; summary(markers2$avg_logFC); table(markers2$avg_logFC > 0)


### ID conversion (from ensembl_gene_id to entrezgene)
library('biomaRt')
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart)   #### select dataset "drerio_gene_ensembl"
#### Get entrezgene id from biomaRt
acinal_cell_markerGenes <- getBM(attribute=c('ensembl_gene_id', 'entrezgene'),filters = 'ensembl_gene_id', values= row.names(markers0),mart = ensembl)
alpha_cell_markerGenes <- getBM(attribute=c('ensembl_gene_id', 'entrezgene'),filters = 'ensembl_gene_id', values= row.names(markers1),mart = ensembl)
beta_cell_markerGenes<- getBM(attribute=c('ensembl_gene_id', 'entrezgene'),filters = 'ensembl_gene_id', values= row.names(markers2),mart = ensembl)
delta_cell_markerGenes<- getBM(attribute=c('ensembl_gene_id', 'entrezgene'),filters = 'ensembl_gene_id', values= row.names(markers3),mart = ensembl)
ductal_cell_markerGenes<- getBM(attribute=c('ensembl_gene_id', 'entrezgene'),filters = 'ensembl_gene_id', values= row.names(markers4),mart = ensembl)
beta_vs_alpha_delta_markerGenes <- getBM(attribute=c('ensembl_gene_id', 'entrezgene'),filters = 'ensembl_gene_id', 
                                         values= row.names(cluster2.markers),mart = ensembl)

#### Get rid of NA genes
acinal_cell_markerGenes_entrez <- acinal_cell_markerGenes$entrezgene[!is.na(acinal_cell_markerGenes$entrezgene)]
alpha_cell_markerGenes_entrez <- alpha_cell_markerGenes$entrezgene[!is.na(alpha_cell_markerGenes$entrezgene)]
beta_cell_markerGenes_entrez <- beta_cell_markerGenes$entrezgene[!is.na(beta_cell_markerGenes$entrezgene)]
delta_cell_markerGenes_entrez <- delta_cell_markerGenes$entrezgene[!is.na(delta_cell_markerGenes$entrezgene)]
ductal_cell_markerGenes_entrez <- ductal_cell_markerGenes$entrezgene[!is.na(ductal_cell_markerGenes$entrezgene)]
beta_vs_alpha_delta_markerGenes <- beta_vs_alpha_delta_markerGenes$entrezgene[!is.na(beta_vs_alpha_delta_markerGenes$entrezgene)]

#### Export entrezgene id of each category
write.table(acinal_cell_markerGenes_entrez, file = 'acinal_cell_markerGenes_entrez.txt')
write.table(alpha_cell_markerGenes_entrez, file = 'alpha_cell_markerGenes_entrez.txt')
write.table(beta_cell_markerGenes_entrez, file = 'beta_cell_markerGenes_entrez.txt')
write.table(delta_cell_markerGenes_entrez, file = 'delta_cell_markerGenes_entrez.txt')
write.table(ductal_cell_markerGenes_entrez, file = 'ductal_cell_markerGenes_entrez.txt')
write.table(beta_vs_alpha_delta_markerGenes, file = 'beta_vs_alpha_delta_markerGenes.txt')




##################### Enrichment analysis #####################

### Restart rstudio
### Reset working directory
setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/lncRNA_EDA')

### Load Packages
readable = T
library("DOSE")
library("clusterProfiler")
library("enrichplot")
library('viridis')

### Enrichment Analysis of each category using 'clusterProfiler', 'DOSE', and 'enrichplot' packages
beta_vs_alpha_delta_markerGenes <- read.table('beta_vs_alpha_delta_markerGenes.txt')
beta_vs_alpha_delta_markerGenes <- as.vector(as.matrix(beta_vs_alpha_delta_markerGenes)); 
length(beta_vs_alpha_delta_markerGenes); class(beta_vs_alpha_delta_markerGenes)

### Using 'clusterprofiler', 'DOSE', 'enrichplot' packages to do enrichment analysis
#### GO DAG graph (有向无环图) -- clusterprofiler, enrichplot
ego_beta <- enrichGO(beta_vs_alpha_delta_markerGenes, OrgDb = 'org.Dr.eg.db', ont = 'BP', readable = T)
goplot(ego_beta)
#### Barplot
barplot(ego_beta, showCategory = 20)
#### Dotplot
dotplot(ego_beta, showCategory = 20) + scale_color_viridis(option = 'magma')
#### Dotplot facet -- ont = 'all
go_beta <- enrichGO(beta_vs_alpha_delta_markerGenes, OrgDb = 'org.Dr.eg.db', ont = 'all', readable = T)
dotplot(go_beta, split = 'ONTOLOGY') + facet_grid(ONTOLOGY ~., scale = 'free')
#### Gene-concept network (cnetplot): sometimes very complicated and not easy to read
ego2 <- simplify(ego_beta); head(ego2)
cnetplot(ego2)
cnetplot(ego2, circular = T)
#### Upset Plot: The upsetplot is an alternative to cnetplot for visualizing the complex association between genes
#### and gene sets. It emphasizes the gene overlapping among different gene sets.
upsetplot(ego_beta)
#### Heatmap-like functional classification: The heatplot is similar to cnetplot, while displaying the relationships
#### as a heatmap. The gene-concept network may become too complicated if user wants to show a large number of 
#### significant terms. The heatplot can simplify the result and more easy to identify expression patterns. The function
#### is from enrichplot package.
heatplot(ego_beta)
#### Enrichment map: Enrichment map organize enriched terms into a network with edges connecting overlapping gene sets.
#### In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional 
#### modules. This function is from enrichplot package.
emapplot(ego_beta)


#### Here we can see, in beta-cell "BP" GO term, peptide biosynthetic process(GO:0043043), 
#### intracellular protein transport(GO:0006886), ribosome assembly(GO:0042255), translation(GO:0006412),
#### ribonucleoprotein complex biogenesis(GO:0022618), ribosome biogenesis(GO:0042254),
#### ribonucleoprotein complex subunit organization(GO:0071826) are highly associated with each other
str(ego_beta); View(ego_beta@result)
View(listFilters(ensembl))
beta_houseKeeping <- getBM(attributes = 'entrezgene', filters = 'go', values = c('GO:0043043','GO:0006886',
                                                                                 'GO:0042255', 'GO:0006412',
                                                                                 'GO:0022618','GO:0042254',
                                                                                 'GO:0071826','GO:0009126',
                                                                                 'GO:0046034', 'GO:0009259',
                                                                                 'GO:0072521','GO:0006753',
                                                                                 'GO:0055086', 'GO:0006754',
), mart = ensembl)
beta_houseKeeping <- getBM(attributes = 'entrezgene', filters = 'go', values = ego_beta@result$ID[1:36], mart = ensembl)
beta_vs_alpha_delta_markerGenes <- beta_vs_alpha_delta_markerGenes[!beta_vs_alpha_delta_markerGenes %in% beta_houseKeeping$entrezgene]


#========================================================================================================
#
#                                         WGCNA for known genes
#
#========================================================================================================
library(WGCNA)
norm_counts <- normalized_counts[,1:20] # only consider four types of cell (omit ductal cells since with no replicates)
datExpr0 <- as.data.frame(t(norm_counts)); dim(datExpr0)

gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ',')))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Remove samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ',')))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

head(datExpr0); dim(datExpr0)
n <- nrow(datExpr0); n
datExpr0[n+1,] <- apply(datExpr0[c(1:n),],2, function(x){log2(mean(x)+1)}); dim(datExpr0)
datExpr0[n+2,] <- apply(datExpr0[c(1:n),], 2, function(x){log2(sd(x)/mean(x)+1)}) ; dim(datExpr0)# 使用变异系数（coefficient of variance, CV)较正表达量高低对变异度的影响
datExpr1 <- as.data.frame(t(datExpr0))
names(datExpr1)[21] <- 'log2_mean';names(datExpr1)[22] <- 'log2_CV'; names(datExpr1)
head(datExpr1)[21]; head(datExpr1[22]); colnames(datExpr1)
summary(datExpr1)
head(datExpr1)


model_xlog2mean_ylog2CV <- loess(datExpr1$log2_CV ~ datExpr1$log2_mean, span = 0.8, method = 'loess')
prediction <- predict(object = model_xlog2mean_ylog2CV, data.frame(datExpr1$log2_mean), se = T)
datExpr0 <- datExpr1[datExpr1$log2_CV > (prediction$fit + 2*prediction$se.fit) & datExpr1$log2_mean > 5,1:20]; dim(datExpr0) 


p <- ggplot(datExpr1, aes(x = log2_mean, y = log2_CV), col = 'black')+ geom_point() + 
  geom_smooth(span = 0.8, method = 'loess', na.rm = T) + 
  geom_smooth(method = lm, col = 'red', na.rm = T) + 
  ylim(c(0,2.8)) +
  xlim(0,10) +
  geom_vline(xintercept = 5, col = 'darkgreen', lty = 2) +
  theme_classic() +
  labs(x = 'Mean of Normalized counts (log2-transformed)', y = 'Coefficient of variance (log2-transformed)') +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16)) +
  geom_point(data = datExpr1[datExpr1$log2_CV > (prediction$fit + 2*prediction$se.fit) & datExpr1$log2_mean > 5,21:22],col = 'red');  p


filtered_TPM_normalized_counts <- datExpr0
head(filtered_TPM_normalized_counts)
dim(filtered_TPM_normalized_counts)


library('gplots')
library('pheatmap')
library('amap')
library('RColorBrewer')
pearson_cor <- as.matrix(cor(x = datExpr0, method = 'pearson'))
head(pearson_cor); dim(pearson_cor)
hc <- hcluster(t(datExpr0), method="pearson")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(pearson_cor, Rowv = as.dendrogram(hc), trace = 'none',symm = T, col = hmcol, main = 'The pearson correlation of each')
pheatmap(pearson_cor, cutree_rows = 4, cutree_cols = 4, annotation_row = sample, annotation_col = sample)



sampleTree <- hclust(dist(t(datExpr0)), method = 'average')
par(mfrow = c(1,1))
plot(sampleTree, main = "Sample clustering to detect outlier")

###### load trait data ######
traitData <- read.csv(file = 'trait_D.csv', header = T, row.names = 1, check.names = F)
head(traitData)
## convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(traitData, signed = F); traitColors;names(traitData)
plotDendroAndColors(sampleTree, traitColors, groupLabels = names(traitData), main = 'Sample dendrogram and trait heatmap')


################ STEP 2: Network Construction #############
######## Select the best soft-thresholding power #########
## Choose a set of soft-thresholding powers
powers <- seq(1,30,2)
## Call the Network Topological Analysis function
sft <- pickSoftThreshold(t(datExpr0), powerVector = powers, verbose = 5)  # This step takes some time.
par(mfrow = c(1,2), mar = c(6.5,8,3,3))
cex1 = 0.9
str(sft)  
# The output of 'sft' is a list object. powerEstimate is the estimated best power
# The second output of 'sft' is fitIndices, which is a matrix. The fifth column, 'mean.k' denote average connectivity.

## Scale-free topological fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]*sft$fitIndices[,2]), xlab = 'Soft Threshold (power)', 
     ylab = 'Scale free Topological Model Fit, signed R^2', type = 'n', main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]*sft$fitIndices[,2]), labels = powers, cex = cex1, col = 'red')
abline(h = 0.9, col = 'blue') 
# The blue line corresponds to using a R^2 cut-off of h
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = 'Soft Threshold (power)', ylab = 'Mean Connectivity', type = 'n', main = paste('Mean Connectivity'))
text(sft$fitIndices[,1], sft$fitIndices[,5],labels = powers, cex = cex1, col = 'red')
## Choose the softPower
softPower <- sft$powerEstimate
softPower
Adjacency <- adjacency(t(datExpr0), power = softPower)
## Turn adjacency matrix into Topological matrix, this step takes some time ##
TOM <- TOMsimilarity(Adjacency)  
dissTOM <- 1-TOM
## Call the hierarchincal clustering function
geneTree <- hclust(as.dist(dissTOM), method = 'average')   # This step takes some time, to calculate the distance in gene pairs.
par(mfrow = c(1,1))
plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity', labels = F, hang = 0.04)
# We like large modules, so we set the minimum module size relatively high
minModuleSize <- 300
# set cutHeight, otherwise the function will set a cutHeight automatically
# The next step may take some time
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = F, minClusterSize = minModuleSize)
table(dynamicMods)
## Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors,'Dynamic Tree Cut', dendroLabels = F, 
                    hang = 0.03, addGuide = T, guideHang = 0.05, main = 'Gene dendrogram and module colors')


## Calculate eigengene
MEList <- moduleEigengenes(t(datExpr0), dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = 'average')
par(mfrow = c(1,1))
plot(METree, main = 'Clustering of module eigengene', xlab = '', sub = '')
MEDissThres = 0.45 # set the threshold to make some branches together
abline(h = MEDissThres, col = 'red')
Merge <- mergeCloseModules(t(datExpr0), dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- Merge$colors
table(mergedColors)
mergedMEs <- Merge$newMEs
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c('Dynamic Tree Cut', 'Merged dynamic'), 
                    dendroLabels = F, hang = 0.03, addGuide = T, guideHang = 0.05) # It takes time!!!
## Rename to module Colors
moduleColors <- mergedColors
colorOrder <- c('grey', standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
moduleLabels
MEs <- mergedMEs



nGene <- ncol(t(datExpr0))
nSample <- nrow(t(datExpr0))
moduleTraitCor <- cor(MEs, traitData, use = 'p')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSample)
textMatrix = paste(signif(moduleTraitCor,2),"\n(", signif(moduleTraitPvalue,1),")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mfrow = c(1,1))
par(mar = c(3,8.5,2,2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 1.25,
               cex.lab = 1.25,
               zlim <- c(-1,1),
               main = paste('Module-trait relationships'))
?labeledHeatmap


nSelect <- 800
set.seed(10)
select <- sample(nGene, size = nSelect)
selectTOM <- dissTOM[select, select]
selectTree <- hclust(as.dist(selectTOM), method = 'average')
selectColors <- moduleColors[select]
plotDiss <- selectTOM^8
diag(plotDiss) <- NA
TOMplot(plotDiss, selectTree, selectColors)
?TOMplot

###################### Visualizing the gene network of eigengene ###################
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8)

##################### Module membership (MM) and Gene significance ######################
# Gene Significance, GS: 基因显著性参数，为非负数字信息，比如基于样本的临床信息(clinical traits)和基于每个基因的-log(p-value)等
# names (colors) of each module
modNames <- substring(names(MEs),3)
modNames
# 首先计算模块与基因的相关性矩阵
# MEs表示每个模块在样本里的值
geneModuleMembership <- as.data.frame(cor(t(datExpr0), MEs, use = 'p'))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSample))
head(MMPvalue)
names(MMPvalue) <- paste("p.MM", modNames, sep = "")
# names of those traits
traitNames <- names(traitData)
geneTraitSignificance <- as.data.frame(cor(t(datExpr0), traitData, use = 'p'))
head(geneTraitSignificance)
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSample))
head(GSPvalue)
names(geneTraitSignificance) <- paste('GS.',traitNames,sep = "")
names(GSPvalue) <- paste("p.GS.", traitNames, sep = "")


# 也可以指定感兴趣的模块进行分析，每一个module都分配了一个color
# 比如对module = 'greenyellow’ 的模块进行分析
# 'green' module gene
module <- c('yellow')
column <- match(module, modNames)
moduleGenes <- moduleColors == module
head(moduleGenes)
head(moduleColors)
test_module_index <- which(moduleColors == module) # the index of genes which belong to 'blue' module
length(colnames(datExpr0)[test_module_index])
length(rownames(filtered_TPM_normalized_counts)[test_module_index])
# 注意datExpr0和filtered_normalized_counts是转置的关系，所以datExpr0的colnames和filtered_normalized_counts的rownames是一致的
# 都是基因名，相当于后面的probes
test_module_transcriptName <- rownames(filtered_TPM_normalized_counts)[test_module_index]
length(test_module_transcriptName); test_module_transcriptName


#===========================================================================
#
#                          Enrichment analysis
#
#===========================================================================

############################### beta cell ##################################
### Only consider 'greenyellow' module

##### Use biomaRt to convert gene ID of differential expressed genes with very significant p-value
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart)
listFilters(ensembl)
beta_genes <- getBM(attribute=c('ensembl_gene_id', 'entrezgene','zfin_id_symbol'),
                          filters = 'ensembl_gene_id', values= test_module_transcriptName, mart = ensembl)
beta_genes <- beta_genes[!is.na(beta_genes$entrezgene),]

#### The input for enrichment analysis is entrezgene id 
beta_genes <- unique(beta_genes$entrezgene)
length(beta_genes);

#### Biological Process
##### Calculate 'BP' GO 
BP <- enrichGO(beta_genes,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
result_BP <- simplify(BP,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_BP
##### dotplot
dotplot(result_BP,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.1)+scale_color_continuous(low = 'purple', high = 'green')  # Use scale_size to change the size of the bubble, if the bubble is too large and be cut by the edge, you can used xlim to extend the x axis
dotplot(result_BP,showCategory = 10,x = 'count')+scale_size(range = c(2,12))+xlim(NA,25) + scale_color_viridis()
##### barplot
barplot(result_BP,drop = T, showCategory = 10) + 
  scale_x_discrete(labels = function(x)str_wrap(x, width = 25)) 
barplot(result_BP,x = 'count', showCategory = 10) + 
  scale_x_discrete(labels = function(x)str_wrap(x, width = 25)) 
goplot(result_BP)
##### Gene-concept network (cnetplot): sometimes very complicated and not easy to read
cnetplot(result_BP)
cnetplot(result_BP, circular = T)

#### Molecular Function
MF <- enrichGO(beta_genes,'org.Dr.eg.db',pvalueCutoff = 0.2, pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'MF', readable = readable)
result_MF <- simplify(MF,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_MF
dotplot(result_MF,showCategory = 10)+scale_size(range = c(2,10))+ggplot2::xlim(NA, 0.04)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_MF,drop = T, showCategory = 10)
plotGOgraph(result_MF)
goplot(result_MF)
#### Cellular Component
CC <- enrichGO(beta_genes,'org.Dr.eg.db',pvalueCutoff = 0.4,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'CC', readable = readable)
result_CC <- simplify(CC,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_CC
dotplot(result_CC,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.05)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_CC,drop = T, showCategory = 10)
plotGOgraph(result_CC)
goplot(result_CC)
#### KEGG pathway
kk <- enrichKEGG(beta_genes,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.4,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
result_kk
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.085)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 5)
#### ReactomePA
pa <- enrichPathway(beta_genes,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
plot(pa)
?enrichPathway

############################### beta cell with "greenyellow" and "brown" two modules ################################

# 也可以指定感兴趣的模块进行分析，每一个module都分配了一个color
# 比如对module = 'greenyellow’ 的模块进行分析
# 'green' module gene
module <- c('greenyellow')
column <- match(module, modNames)
moduleGenes <- moduleColors == module
head(moduleGenes)
head(moduleColors)
test_module_index <- which(moduleColors == module) # the index of genes which belong to 'blue' module
length(colnames(datExpr0)[test_module_index])
length(rownames(filtered_TPM_normalized_counts)[test_module_index])
# 注意datExpr0和filtered_normalized_counts是转置的关系，所以datExpr0的colnames和filtered_normalized_counts的rownames是一致的
# 都是基因名，相当于后面的probes
test_module_transcriptName <- rownames(filtered_TPM_normalized_counts)[test_module_index]
length(test_module_transcriptName); test_module_transcriptName



# 也可以指定感兴趣的模块进行分析，每一个module都分配了一个color
# 比如对module = 'pink’ 的模块进行分析
# 'pink' module gene
module <- c('pink')
column <- match(module, modNames)
moduleGenes <- moduleColors == module
head(moduleGenes)
head(moduleColors)
test_module_index <- which(moduleColors == module) # the index of genes which belong to 'blue' module
length(colnames(datExpr0)[test_module_index])
length(rownames(filtered_TPM_normalized_counts)[test_module_index])
# 注意datExpr0和filtered_normalized_counts是转置的关系，所以datExpr0的colnames和filtered_normalized_counts的rownames是一致的
# 都是基因名，相当于后面的probes
test_module_transcriptName1 <- rownames(filtered_TPM_normalized_counts)[test_module_index]
length(test_module_transcriptName1); test_module_transcriptName1

test_module_transcriptName <- c(test_module_transcriptName1, test_module_transcriptName)
length(test_module_transcriptName)

### consider 'greenyellow' and 'purple' two modules together

beta_genes <- getBM(attribute=c('ensembl_gene_id', 'entrezgene','zfin_id_symbol'),
                    filters = 'ensembl_gene_id', values= test_module_transcriptName, mart = ensembl)
beta_genes <- beta_genes[!is.na(beta_genes$entrezgene),]

#### The input for enrichment analysis is entrezgene id 
beta_genes <- unique(beta_genes$entrezgene)


#### Biological Process
##### Calculate 'BP' GO 
BP <- enrichGO(beta_genes,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
result_BP <- simplify(BP,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_BP
##### dotplot
dotplot(result_BP,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.06)+scale_color_continuous(low = 'purple', high = 'green')  # Use scale_size to change the size of the bubble, if the bubble is too large and be cut by the edge, you can used xlim to extend the x axis
dotplot(result_BP,showCategory = 10,x = 'count')+scale_size(range = c(2,12))+xlim(NA,30) + scale_color_viridis()
##### barplot
barplot(result_BP,drop = T, showCategory = 10) + 
  scale_x_discrete(labels = function(x)str_wrap(x, width = 25)) 
barplot(result_BP,x = 'count', showCategory = 10) + 
  scale_x_discrete(labels = function(x)str_wrap(x, width = 25)) 
goplot(result_BP)
##### Gene-concept network (cnetplot): sometimes very complicated and not easy to read
cnetplot(result_BP)
cnetplot(result_BP, circular = T)

#### Molecular Function
MF <- enrichGO(beta_genes,'org.Dr.eg.db',pvalueCutoff = 0.05, pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'MF', readable = readable)
result_MF <- simplify(MF,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_MF
dotplot(result_MF,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.04)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_MF,drop = T, showCategory = 10)
plotGOgraph(result_MF)
goplot(result_MF)
#### Cellular Component
CC <- enrichGO(beta_genes,'org.Dr.eg.db',pvalueCutoff = 0.2,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'CC', readable = readable)
result_CC <- simplify(CC,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_CC
dotplot(result_CC,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.05)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_CC,drop = T, showCategory = 10)
plotGOgraph(result_CC)
goplot(result_CC)
#### KEGG pathway
kk <- enrichKEGG(beta_genes,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.4,pAdjustMethod = 'BH',qvalueCutoff = 0.5)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
result_kk
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.085)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 5)
#### ReactomePA
pa <- enrichPathway(beta_genes,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)






############################### alpha cell ##################################
# 也可以指定感兴趣的模块进行分析，每一个module都分配了一个color
# 比如对module = 'turquoise’ 的模块进行分析
# 'turquoise' module gene
module <- c('blue')
column <- match(module, modNames)
moduleGenes <- moduleColors == module
head(moduleGenes)
head(moduleColors)
test_module_index <- which(moduleColors == module) # the index of genes which belong to 'blue' module
length(colnames(datExpr0)[test_module_index])
length(rownames(filtered_TPM_normalized_counts)[test_module_index])
# 注意datExpr0和filtered_normalized_counts是转置的关系，所以datExpr0的colnames和filtered_normalized_counts的rownames是一致的
# 都是基因名，相当于后面的probes
test_module_transcriptName <- rownames(filtered_TPM_normalized_counts)[test_module_index]
length(test_module_transcriptName); test_module_transcriptName


alpha_genes <- getBM(attribute=c('ensembl_gene_id', 'entrezgene','zfin_id_symbol'),
                    filters = 'ensembl_gene_id', values= test_module_transcriptName, mart = ensembl)
alpha_genes <- alpha_genes[!is.na(alpha_genes$entrezgene),]

#### The input for enrichment analysis is entrezgene id 
alpha_genes <- unique(alpha_genes$entrezgene)


#### Biological Process
##### Calculate 'BP' GO 
BP <- enrichGO(alpha_genes,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
result_BP <- simplify(BP,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_BP
##### dotplot
dotplot(result_BP,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.075)+scale_color_continuous(low = 'purple', high = 'green')  # Use scale_size to change the size of the bubble, if the bubble is too large and be cut by the edge, you can used xlim to extend the x axis
dotplot(result_BP,showCategory = 10,x = 'count')+scale_size(range = c(2,12))+xlim(NA,50) + scale_color_viridis()
##### barplot
barplot(result_BP,drop = T, showCategory = 10) + 
  scale_x_discrete(labels = function(x)str_wrap(x, width = 25)) 
barplot(result_BP,x = 'count', showCategory = 10) + 
  scale_x_discrete(labels = function(x)str_wrap(x, width = 25)) 
plotGOgraph(result_BP)
goplot(result_BP)
##### Gene-concept network (cnetplot): sometimes very complicated and not easy to read
cnetplot(result_BP)
cnetplot(result_BP, circular = T)
?cnetplot

#### Molecular Function
MF <- enrichGO(alpha_genes,'org.Dr.eg.db',pvalueCutoff = 0.05, pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'MF', readable = readable)
result_MF <- simplify(MF,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_MF
dotplot(result_MF,showCategory = 10)+scale_size(range = c(2,10))+ggplot2::xlim(NA, 0.07)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_MF,drop = T, showCategory = 10)
plotGOgraph(result_MF)
goplot(result_MF)
cnetplot(result_MF)
#### Cellular Component
CC <- enrichGO(alpha_genes,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'CC', readable = readable)
result_CC <- simplify(CC,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_CC
dotplot(result_CC,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.045)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_CC,drop = T, showCategory = 10)
plotGOgraph(result_CC)
goplot(result_CC)
cnetplot(result_CC)
#### KEGG pathway
kk <- enrichKEGG(alpha_genes,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.5)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
result_kk
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_kk,drop = T, showCategory = 10)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 5)
#### ReactomePA
pa <- enrichPathway(alpha_genes,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
result_pa <- setReadable(pa,'org.Dr.eg.db',keytype = 'ENTREZID')
result_pa
dotplot(result_pa)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.13)+scale_color_continuous(low = 'purple', high = 'green') 


############################### delta cell ##################################
# 也可以指定感兴趣的模块进行分析，每一个module都分配了一个color
# 比如对module = 'red’ 的模块进行分析
# 'red' module gene
module <- c('turquoise')
column <- match(module, modNames)
moduleGenes <- moduleColors == module
head(moduleGenes)
head(moduleColors)
test_module_index <- which(moduleColors == module) # the index of genes which belong to 'blue' module
length(colnames(datExpr0)[test_module_index])
length(rownames(filtered_TPM_normalized_counts)[test_module_index])
# 注意datExpr0和filtered_normalized_counts是转置的关系，所以datExpr0的colnames和filtered_normalized_counts的rownames是一致的
# 都是基因名，相当于后面的probes
test_module_transcriptName <- rownames(filtered_TPM_normalized_counts)[test_module_index]
length(test_module_transcriptName); test_module_transcriptName


delta_genes <- getBM(attribute=c('ensembl_gene_id', 'entrezgene','zfin_id_symbol'),
                     filters = 'ensembl_gene_id', values= test_module_transcriptName, mart = ensembl)
delta_genes <- delta_genes[!is.na(delta_genes$entrezgene),]

#### The input for enrichment analysis is entrezgene id 
delta_genes <- unique(delta_genes$entrezgene)


#### Biological Process
##### Calculate 'BP' GO 
BP <- enrichGO(delta_genes,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
result_BP <- simplify(BP,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_BP
##### dotplot
dotplot(result_BP,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.07)+scale_color_continuous(low = 'purple', high = 'green')  # Use scale_size to change the size of the bubble, if the bubble is too large and be cut by the edge, you can used xlim to extend the x axis
dotplot(result_BP,showCategory = 10,x = 'count')+scale_size(range = c(2,12))+xlim(NA,60) + scale_color_viridis()
##### barplot
barplot(result_BP,drop = T, showCategory = 10) + 
  scale_x_discrete(labels = function(x)str_wrap(x, width = 25)) 
barplot(result_BP,x = 'count', showCategory = 10) + 
  scale_x_discrete(labels = function(x)str_wrap(x, width = 25)) 
plotGOgraph(result_BP)
goplot(result_BP)
##### Gene-concept network (cnetplot): sometimes very complicated and not easy to read
cnetplot(result_BP)
cnetplot(result_BP, circular = T)
?cnetplot

#### Molecular Function
MF <- enrichGO(delta_genes,'org.Dr.eg.db',pvalueCutoff = 0.05, pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'MF', readable = readable)
result_MF <- simplify(MF,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_MF
dotplot(result_MF,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.08)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_MF,drop = T, showCategory = 10)
plotGOgraph(result_MF)
goplot(result_MF)
cnetplot(result_MF)
#### Cellular Component
CC <- enrichGO(delta_genes,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'CC', readable = readable)
result_CC <- simplify(CC,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_CC
dotplot(result_CC,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.055)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_CC,drop = T, showCategory = 10)
plotGOgraph(result_CC)
goplot(result_CC)
cnetplot(result_CC)
#### KEGG pathway
kk <- enrichKEGG(delta_genes,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.5) 
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
result_kk
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.135)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_kk,drop = T, showCategory = 10)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 5)
#### ReactomePA
pa <- enrichPathway(delta_genes,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
result_pa <- setReadable(pa,'org.Dr.eg.db',keytype = 'ENTREZID')
result_pa
dotplot(result_pa)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.1115)+scale_color_continuous(low = 'purple', high = 'green') 


############################### acinal cell ##################################
# 也可以指定感兴趣的模块进行分析，每一个module都分配了一个color
# 比如对module = 'brown’ 的模块进行分析
# 'red' module gene
module <- c('brown')
column <- match(module, modNames)
moduleGenes <- moduleColors == module
head(moduleGenes)
head(moduleColors)
test_module_index <- which(moduleColors == module) # the index of genes which belong to 'blue' module
length(colnames(datExpr0)[test_module_index])
length(rownames(filtered_TPM_normalized_counts)[test_module_index])
# 注意datExpr0和filtered_normalized_counts是转置的关系，所以datExpr0的colnames和filtered_normalized_counts的rownames是一致的
# 都是基因名，相当于后面的probes
test_module_transcriptName <- rownames(filtered_TPM_normalized_counts)[test_module_index]
length(test_module_transcriptName); test_module_transcriptName


acinal_genes <- getBM(attribute=c('ensembl_gene_id', 'entrezgene','zfin_id_symbol'),
                     filters = 'ensembl_gene_id', values= test_module_transcriptName, mart = ensembl)
acinal_genes <- acinal_genes[!is.na(acinal_genes$entrezgene),]

#### The input for enrichment analysis is entrezgene id 
acinal_genes <- unique(acinal_genes$entrezgene)


#### Biological Process
##### Calculate 'BP' GO 
BP <- enrichGO(acinal_genes,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
result_BP <- simplify(BP,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_BP
##### dotplot
dotplot(result_BP,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.11)+scale_color_continuous(low = 'purple', high = 'green')  # Use scale_size to change the size of the bubble, if the bubble is too large and be cut by the edge, you can used xlim to extend the x axis
dotplot(result_BP,showCategory = 10,x = 'count')+scale_size(range = c(2,12))+xlim(NA,110) + scale_color_viridis()
##### barplot
barplot(result_BP,drop = T, showCategory = 10) + 
  scale_x_discrete(labels = function(x)str_wrap(x, width = 25)) 
barplot(result_BP,x = 'count', showCategory = 10) + 
  scale_x_discrete(labels = function(x)str_wrap(x, width = 25)) 
plotGOgraph(result_BP)
goplot(result_BP)
##### Gene-concept network (cnetplot): sometimes very complicated and not easy to read
cnetplot(result_BP)
cnetplot(result_BP, circular = T)

#### Molecular Function
MF <- enrichGO(acinal_genes,'org.Dr.eg.db',pvalueCutoff = 0.05, pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'MF', readable = readable)
result_MF <- simplify(MF,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_MF
dotplot(result_MF,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.065)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_MF,drop = T, showCategory = 10)
plotGOgraph(result_MF)
goplot(result_MF)
cnetplot(result_MF)
#### Cellular Component
CC <- enrichGO(acinal_genes,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'CC', readable = readable)
result_CC <- simplify(CC,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_CC
dotplot(result_CC,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.065)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_CC,drop = T, showCategory = 10)
plotGOgraph(result_CC)
goplot(result_CC)
cnetplot(result_CC)
#### KEGG pathway
kk <- enrichKEGG(acinal_genes,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.5)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
result_kk
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.135)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_kk,drop = T, showCategory = 10)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 5)
#### ReactomePA
pa <- enrichPathway(acinal_genes,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
result_pa <- setReadable(pa,'org.Dr.eg.db',keytype = 'ENTREZID')
result_pa
dotplot(result_pa)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.155)+scale_color_continuous(low = 'purple', high = 'green') 



############################### grey module ##################################
# 也可以指定感兴趣的模块进行分析，每一个module都分配了一个color
# 比如对module = 'grey’ 的模块进行分析
# 'red' module gene
module <- c('grey')
column <- match(module, modNames)
moduleGenes <- moduleColors == module
head(moduleGenes)
head(moduleColors)
test_module_index <- which(moduleColors == module) # the index of genes which belong to 'blue' module
length(colnames(datExpr0)[test_module_index])
length(rownames(filtered_TPM_normalized_counts)[test_module_index])
# 注意datExpr0和filtered_normalized_counts是转置的关系，所以datExpr0的colnames和filtered_normalized_counts的rownames是一致的
# 都是基因名，相当于后面的probes
test_module_transcriptName <- rownames(filtered_TPM_normalized_counts)[test_module_index]
length(test_module_transcriptName); test_module_transcriptName


grey_genes <- getBM(attribute=c('ensembl_gene_id', 'entrezgene','zfin_id_symbol'),
                      filters = 'ensembl_gene_id', values= test_module_transcriptName, mart = ensembl)
grey_genes <- grey_genes[!is.na(grey_genes$entrezgene),]

#### The input for enrichment analysis is entrezgene id 
grey_genes <- unique(grey_genes$entrezgene)


#### Biological Process
##### Calculate 'BP' GO 
BP <- enrichGO(grey_genes,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
result_BP <- simplify(BP,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_BP
##### dotplot
dotplot(result_BP,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.22)+scale_color_continuous(low = 'purple', high = 'green')  # Use scale_size to change the size of the bubble, if the bubble is too large and be cut by the edge, you can used xlim to extend the x axis
dotplot(result_BP,showCategory = 10,x = 'count')+scale_size(range = c(2,12))+xlim(NA,110) + scale_color_viridis()
##### barplot
barplot(result_BP,drop = T, showCategory = 10) + 
  scale_x_discrete(labels = function(x)str_wrap(x, width = 25)) 
barplot(result_BP,x = 'count', showCategory = 10) + 
  scale_x_discrete(labels = function(x)str_wrap(x, width = 25)) 
plotGOgraph(result_BP)
goplot(result_BP)
##### Gene-concept network (cnetplot): sometimes very complicated and not easy to read
cnetplot(result_BP)
cnetplot(result_BP, circular = T)

#### Molecular Function
MF <- enrichGO(grey_genes,'org.Dr.eg.db',pvalueCutoff = 0.05, pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'MF', readable = readable)
result_MF <- simplify(MF,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_MF
dotplot(result_MF,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.25)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_MF,drop = T, showCategory = 10)
plotGOgraph(result_MF)
goplot(result_MF)
cnetplot(result_MF)
#### Cellular Component
CC <- enrichGO(grey_genes,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'CC', readable = readable)
result_CC <- simplify(CC,cutoff = 0.7, by = 'p.adjust', select_fun = min)
result_CC
dotplot(result_CC,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.22)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_CC,drop = T, showCategory = 10)
plotGOgraph(result_CC)
goplot(result_CC)
cnetplot(result_CC)
#### KEGG pathway
kk <- enrichKEGG(grey_genes,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.5)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
result_kk
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.28)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_kk,drop = T, showCategory = 10)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 5)
#### ReactomePA
pa <- enrichPathway(grey_genes,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
result_pa <- setReadable(pa,'org.Dr.eg.db',keytype = 'ENTREZID')
result_pa
dotplot(result_pa)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.21)+scale_color_continuous(low = 'purple', high = 'green') 


########################### Exporting to Cytoscape all one by one ##########################
## Reset workding directory to export the node and edge data
setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/lncRNA_data/WGCNA/WGCNA/WGCNA')
## Select each module
for (mod in 1:nrow(table(moduleColors)))
{
  
  modules = names(table(moduleColors))[mod]
  # Select module probes
  probes = colnames(t(datExpr0))
  inModule = (moduleColors == modules)
  modProbes = probes[inModule]
  modGenes = modProbes
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule]
  
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-TPM-edges-", modules , ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-TPM-nodes-", modules, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
}





#===================================================================================================
#
#                          Differential expressed annotated antisense lncRNA
#
#===================================================================================================
setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/lncRNA_EDA')
reads <- read.csv('counts_total.csv', header = T, row.names = 1)
reads[1:3,1:3];dim(reads)
reads <- reads[,2:15]
## In total, 23 samples with 32266 genes detected

## Read in Coldata(Phenodata)
sample <- read.csv('SampleGroup_endocrine.csv', header = T, row.names = 1, colClasses = 'factor')
sample

## Read in metadata of lincRNA and antisense (lncRNA)
setwd("/Users/mijiarui/RNA-seq/CAGE")
lincRNA <- read.table(file = "annotated.lincRNA.bed.txt", header = F, sep = '\t', quote = "")
colnames(lincRNA) <- c('chr', 'start', 'end', 'ensembl', 'geneid', 'strand')
row.names(lincRNA) <- lincRNA$ensembl
head(lincRNA); dim(lincRNA)
antisense <- read.table(file = "annotated.antisense.bed.txt", header = F, sep = '\t', quote = "")
colnames(antisense) <- c('chr', 'start', 'end', 'ensembl', 'geneid', 'strand')
row.names(antisense) <- antisense$ensembl
head(antisense); dim(antisense)

### expression matrix of lincRNA and antisense (lncRNA)
lincRNA_exprs <- reads[row.names(reads) %in% row.names(lincRNA),]
head(lincRNA_exprs); dim(lincRNA_exprs)

antisense_exprs <- reads[row.names(reads) %in% row.names(antisense),]
head(antisense_exprs); dim(antisense_exprs)



# Make DESeq object and Do the normalization (size factor)
## Create DESeqDataSet
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = antisense_exprs, colData = sample, design  = ~ celltype)

## Calculation and Normalization: get normalized_counts
dds <- DESeq(ddsFullCountTable)
normalized_counts <- counts(dds, normalized= T)
class(normalized_counts); head(normalized_counts)

## Sort according to mad
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts < normalized_counts[order(normalized_counts_mad, decreasing = T),]
## Log transformation
rld <- rlog(dds, blind = F) # This step takes some time, we do not put '+1' here
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing = T),]



## Hierarchical clustering
### Selecting colors
display.brewer.all(type = 'all') 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
### Calculating the pearson correlation matrix
#### Clustering is based on the relationship among each individual sample. Both distance and correlation
#### are accepted but correlation is preferred. The distance methods and correlation methods are chosen 
#### based on your data.
pearson_cor <- as.matrix(cor(rlogMat, method = 'pearson'))
head(pearson_cor); dim(pearson_cor)
hc <- hcluster(t(rlogMat), method="pearson")
#### Here I use 'heatmap.2' function from 'gplots' package.
heatmap.2(pearson_cor, Rowv = as.dendrogram(hc), trace = 'none',symm = T, col = hmcol, main = 'The pearson correlation of each')
pheatmap(pearson_cor)



#======================================================================================
#             Beta-cell vs alpha-cell and Beta-cell vs delta-cell
#======================================================================================

################### Perform differential gene expression ################
############### Differential Gene Expression #################
sampleA <- 'beta'
sampleB <- 'acinal'
sampleC <- 'alpha'
sampleD <- 'delta'


### Differential gene expression between sampleA (beta-cell) and sampleC (alpha-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleA, sampleC)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA


#### Calculate the mean value of gene expression in sampleC 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanA, baseMeanC, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)



#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleC,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


signifiant_genes_beta_vs_alpha <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1 ))

length(signifiant_genes_beta_vs_alpha)


### Differential gene expression between sampleA (beta-cell) and sampleD (delta-cell)
#### Comparison between sampleA and sampleD
contrastV <- c('celltype', sampleA, sampleD)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA


#### Calculate the mean value of gene expression in sampleC 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanA, baseMeanD, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)


#### Extract differential expressed gene
colnames(res)
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleD,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


signifiant_genes_beta_vs_delta <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1))

length(signifiant_genes_beta_vs_delta)


#### Pick up the lincRNAs that are differentially expressed in between beta-cell and alpha-cell and
#### between beta-cell and delta-cell
b <- signifiant_genes_beta_vs_delta[signifiant_genes_beta_vs_delta %in% signifiant_genes_beta_vs_alpha]
length(b)

coordinate_b_antisense <- antisense[antisense$ensembl %in% b,]
coordinate_b_antisense


coordinate_b_antisense$chr <- as.numeric(sub("chr","",coordinate_b_antisense$chr))


library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_b_antisense)){
  coordinate_b_antisense$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(coordinate_b_antisense$chr[i],coordinate_b_antisense$start[i]-200000,coordinate_b_antisense$end[i]+200000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_b_antisense$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)


kk <- enrichKEGG(entrez,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)




#======================================================================================
#             alpha-cell vs beta-cell and alpha-cell vs delta-cell
#======================================================================================

### Differential gene expression between sampleC (alpha-cell) and sampleA (beta-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleC, sampleA)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC


#### Calculate the mean value of gene expression in sampleC 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanC, baseMeanA, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)


#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleC,sampleA,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_alpha_vs_beta <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1))

length(signifiant_genes_alpha_vs_beta)





### Differential gene expression between sampleA (beta-cell) and sampleD (delta-cell)
#### Comparison between sampleA and sampleD
contrastV <- c('celltype', sampleC, sampleD)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC


#### Calculate the mean value of gene expression in sampleC 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanC, baseMeanD, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)


#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleC,sampleD,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_alpha_vs_delta <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1))

length(signifiant_genes_alpha_vs_delta)


#### Pick up the lincRNAs that are differentially expressed in between beta-cell and alpha-cell and
#### between beta-cell and delta-cell
a <- signifiant_genes_alpha_vs_delta[signifiant_genes_alpha_vs_delta %in% signifiant_genes_alpha_vs_beta]
length(a)

coordinate_a_antisense <- antisense[antisense$ensembl %in% a,]
coordinate_a_antisense


coordinate_a_antisense$chr <- as.numeric(sub("chr","",coordinate_a_antisense$chr))


library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_a_antisense)){
  coordinate_a_antisense$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(coordinate_a_antisense$chr[i],coordinate_a_antisense$start[i]-500000,coordinate_a_antisense$end[i]+500000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_a_antisense$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)

kk <- enrichKEGG(entrez,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)

pa <- enrichPathway(entrez,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
plot(pa)


#======================================================================================
#             delta-cell vs beta-cell and delta-cell vs alpha-cell
#======================================================================================

### Differential gene expression between sampleD (delta-cell) and sampleA (beta-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleD, sampleA)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD


#### Calculate the mean value of gene expression in sampleC 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanD, baseMeanA, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)



#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleD,sampleA,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_delta_vs_beta <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1))

length(signifiant_genes_delta_vs_beta)


### Differential gene expression between sampleD (delta-cell) and sampleC (alpha-cell)
#### Comparison between sampleA and sampleD
contrastV <- c('celltype', sampleD, sampleC)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleC 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD

#### Calculate the mean value of gene expression in sampleA 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC



#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanD, baseMeanC, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)





#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleD,sampleC,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_delta_vs_alpha <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1 ))

length(signifiant_genes_delta_vs_alpha)


#### Pick up the lincRNAs that are differentially expressed in between beta-cell and alpha-cell and
#### between beta-cell and delta-cell
d <- signifiant_genes_delta_vs_beta[signifiant_genes_delta_vs_beta %in% signifiant_genes_delta_vs_alpha]
length(d)

coordinate_d_antisense <- antisense[antisense$ensembl %in% d,]
coordinate_d_antisense


coordinate_d_antisense$chr <- as.numeric(sub("chr","",coordinate_d_antisense$chr))

library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_d_antisense)){
  coordinate_d_antisense$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                      values= list(coordinate_d_antisense$chr[i],coordinate_d_antisense$start[i]-500000,coordinate_d_antisense$end[i]+500000), 
                                      mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_d_antisense$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)

kk <- enrichKEGG(entrez,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)

pa <- enrichPathway(entrez,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
plot(pa)





abd_ensembl <- c(a,b,d)

pheatmap(rlogMat[rownames(rlogMat) %in% abd_ensembl,], color = colorRampPalette(c('navy', 'white', 'firebrick3'))(50), 
         cluster_rows = T, scale = 'row', cluster_cols = T, annotation_col = sample)




#======================================================================================
#             alpha-cell vs beta-cell and delta-cell vs beta-cell
#======================================================================================
### Differential gene expression between sampleC (alpha-cell) and sampleA (beta-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleC, sampleA)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC


#### Calculate the mean value of gene expression in sampleC 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanC, baseMeanA, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)


#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleC,sampleA,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_alpha_vs_beta <- row.names(subset(res, padj < 0.01 & res$log2FoldChange > 2))

length(signifiant_genes_alpha_vs_beta)


### Differential gene expression between sampleD (delta-cell) and sampleA (beta-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleD, sampleA)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD


#### Calculate the mean value of gene expression in sampleC 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanD, baseMeanA, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)



#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleD,sampleA,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_delta_vs_beta <- row.names(subset(res, padj < 0.01 & res$log2FoldChange > 1))

length(signifiant_genes_delta_vs_beta)


#### Pick up the lincRNAs that are differentially expressed in between beta-cell and alpha-cell and
#### between beta-cell and delta-cell
b <- signifiant_genes_alpha_vs_beta[signifiant_genes_alpha_vs_beta %in% signifiant_genes_delta_vs_beta]
length(b)

coordinate_b <- antisense[antisense$ensembl %in% b,]
coordinate_b


coordinate_b$chr <- as.numeric(sub("chr","",coordinate_b$chr))


for (i in 1:nrow(coordinate_b)){
  coordinate_b$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(coordinate_b$chr[i],coordinate_b$start[i]-100000,coordinate_b$end[i]+100000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_b$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)


#======================================================================================
#             beta-cell vs alpha-cell and delta-cell vs alpha-cell
#======================================================================================
### Differential gene expression between sampleA (beta-cell) and sampleC (alpha-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleA, sampleC)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA


#### Calculate the mean value of gene expression in sampleC 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanA, baseMeanC, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)



#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleC,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


signifiant_genes_beta_vs_alpha <- row.names(subset(res, padj < 0.01 & res$log2FoldChange > 1 ))

length(signifiant_genes_beta_vs_alpha)

### Differential gene expression between sampleD (delta-cell) and sampleC (alpha-cell)
#### Comparison between sampleA and sampleD
contrastV <- c('celltype', sampleD, sampleC)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleC 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD

#### Calculate the mean value of gene expression in sampleA 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC



#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanD, baseMeanC, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)





#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleD,sampleC,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_delta_vs_alpha <- row.names(subset(res, padj < 0.01 & res$log2FoldChange > 1 ))

length(signifiant_genes_delta_vs_alpha)


#### Pick up the lincRNAs that are differentially expressed in between beta-cell and alpha-cell and
#### between beta-cell and delta-cell
a <- signifiant_genes_beta_vs_alpha[signifiant_genes_beta_vs_alpha %in% signifiant_genes_delta_vs_alpha]
length(a)

coordinate_a <- antisense[antisense$ensembl %in% a,]
coordinate_a


coordinate_a$chr <- as.numeric(sub("chr","",coordinate_a$chr))


library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_a)){
  coordinate_a$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(coordinate_a$chr[i],coordinate_a$start[i]-100000,coordinate_a$end[i]+100000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_a$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)


#======================================================================================
#             alpha-cell vs delta-cell and beta-cell vs delta-cell
#======================================================================================
### Differential gene expression between sampleA (beta-cell) and sampleD (delta-cell)
#### Comparison between sampleA and sampleD
contrastV <- c('celltype', sampleC, sampleD)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC


#### Calculate the mean value of gene expression in sampleC 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanC, baseMeanD, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)


#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleC,sampleD,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_alpha_vs_delta <- row.names(subset(res, padj < 0.01 & res$log2FoldChange > 2 ))

length(signifiant_genes_alpha_vs_delta)

### Differential gene expression between sampleA (beta-cell) and sampleD (delta-cell)
#### Comparison between sampleA and sampleD
contrastV <- c('celltype', sampleA, sampleD)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA


#### Calculate the mean value of gene expression in sampleC 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanA, baseMeanD, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)


#### Extract differential expressed gene
colnames(res)
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleD,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


signifiant_genes_beta_vs_delta <- row.names(subset(res, padj < 0.01 & res$log2FoldChange > 1))

length(signifiant_genes_beta_vs_delta)

d <- signifiant_genes_delta_vs_beta[signifiant_genes_delta_vs_beta %in% signifiant_genes_delta_vs_alpha]
length(d)

for (i in 1:nrow(coordinate_d)){
  coordinate_b$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(coordinate_d$chr[i],coordinate_d$start[i]-100000,coordinate_d$end[i]+100000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_d$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)


























#===================================================================================================
#
#                          Differential expressed annotated lincRNA
#
#===================================================================================================
setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/lncRNA_EDA')
reads <- read.csv('counts_total.csv', header = T, row.names = 1)
reads[1:3,1:3];dim(reads)
reads <- reads[,2:15]
## In total, 23 samples with 32266 genes detected

## Read in Coldata(Phenodata)
sample <- read.csv('SampleGroup_endocrine.csv', header = T, row.names = 1, colClasses = 'factor')
sample

## Read in metadata of lincRNA and antisense (lncRNA)
setwd("/Users/mijiarui/RNA-seq/CAGE")
lincRNA <- read.table(file = "annotated.lincRNA.bed.txt", header = F, sep = '\t', quote = "")
colnames(lincRNA) <- c('chr', 'start', 'end', 'ensembl', 'geneid', 'strand')
row.names(lincRNA) <- lincRNA$ensembl
head(lincRNA); dim(lincRNA)
antisense <- read.table(file = "annotated.antisense.bed.txt", header = F, sep = '\t', quote = "")
colnames(antisense) <- c('chr', 'start', 'end', 'ensembl', 'geneid', 'strand')
row.names(antisense) <- antisense$ensembl
head(antisense); dim(antisense)

### expression matrix of lincRNA and antisense (lncRNA)
lincRNA_exprs <- reads[row.names(reads) %in% row.names(lincRNA),]
head(lincRNA_exprs); dim(lincRNA_exprs)

antisense_exprs <- reads[row.names(reads) %in% row.names(antisense),]
head(antisense_exprs); dim(antisense_exprs)



# Make DESeq object and Do the normalization (size factor)
## Create DESeqDataSet, change the count data (antisense or lincRNA)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = lincRNA_exprs, colData = sample, design  = ~ celltype)

## Calculation and Normalization: get normalized_counts
dds <- DESeq(ddsFullCountTable)
normalized_counts <- counts(dds, normalized= T)
class(normalized_counts); head(normalized_counts)

## Sort according to mad
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts < normalized_counts[order(normalized_counts_mad, decreasing = T),]
## Log transformation
rld <- rlog(dds, blind = F) # This step takes some time, we do not put '+1' here
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing = T),]



## Hierarchical clustering
### Selecting colors
display.brewer.all(type = 'all') 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
### Calculating the pearson correlation matrix
#### Clustering is based on the relationship among each individual sample. Both distance and correlation
#### are accepted but correlation is preferred. The distance methods and correlation methods are chosen 
#### based on your data.
pearson_cor <- as.matrix(cor(rlogMat, method = 'pearson'))
head(pearson_cor); dim(pearson_cor)
hc <- hcluster(t(rlogMat), method="pearson")
#### Here I use 'heatmap.2' function from 'gplots' package.
heatmap.2(pearson_cor, Rowv = as.dendrogram(hc), trace = 'none',symm = T, col = hmcol, main = 'The pearson correlation of each')
pheatmap(pearson_cor)



#======================================================================================
#             Beta-cell vs alpha-cell and Beta-cell vs delta-cell
#======================================================================================

################### Perform differential gene expression ################
############### Differential Gene Expression #################
sampleA <- 'beta'
sampleB <- 'acinal'
sampleC <- 'alpha'
sampleD <- 'delta'


### Differential gene expression between sampleA (beta-cell) and sampleC (alpha-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleA, sampleC)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA


#### Calculate the mean value of gene expression in sampleC 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanA, baseMeanC, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)



#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleC,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


signifiant_genes_beta_vs_alpha <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1 ))

length(signifiant_genes_beta_vs_alpha)


### Differential gene expression between sampleA (beta-cell) and sampleD (delta-cell)
#### Comparison between sampleA and sampleD
contrastV <- c('celltype', sampleA, sampleD)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA


#### Calculate the mean value of gene expression in sampleC 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanA, baseMeanD, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)


#### Extract differential expressed gene
colnames(res)
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleD,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


signifiant_genes_beta_vs_delta <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1))

length(signifiant_genes_beta_vs_delta)


#### Pick up the lincRNAs that are differentially expressed in between beta-cell and alpha-cell and
#### between beta-cell and delta-cell
b <- signifiant_genes_beta_vs_delta[signifiant_genes_beta_vs_delta %in% signifiant_genes_beta_vs_alpha]
length(b)

coordinate_b_lincRNA <- lincRNA[lincRNA$ensembl %in% b,]
coordinate_b_lincRNA

coordinate_b_lincRNA$chr <- as.numeric(sub("chr","",coordinate_b_lincRNA$chr))


library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_b_lincRNA)){
  coordinate_b_lincRNA$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                      values= list(coordinate_b_lincRNA$chr[i],coordinate_b_lincRNA$start[i]-200000,coordinate_b_lincRNA$end[i]+200000), 
                                      mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_b_lincRNA$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)

#### integrate lincRNA and antisense lncRNA in beta cell
coordinate_b <- rbind(coordinate_b_antisense, coordinate_b_lincRNA)
head(coordinate_b); dim(coordinate_b)
for (i in 1:nrow(coordinate_b)){
  coordinate_b$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                            values= list(coordinate_b$chr[i],coordinate_b$start[i]-1000000,coordinate_b$end[i]+1000000), 
                                            mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_b$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)

kk <- enrichKEGG(entrez,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)

pa <- enrichPathway(entrez,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
barplot(pa)



################## Pick up the promoter regions ###################
### Normally we define the promoter regions as 500 bp upstream of TSS. Because this is unstranded library, we need to check
### the both ends. It is good to use 'dplyr' package here.
library('dplyr')
test_beta_coordinate <- coordinate_b
promoter_left_strand <- mutate(.data = test_beta_coordinate, Start = start -100, End = start-90)[, c(1,8,9,5)]
promoter_left_strand
promoter_right_strand <- mutate(.data = test_beta_coordinate, Start = end + 90, End = end + 100 )[,c(1,8,9,5)]
promoter_right_strand

######
####### ++/--
####### 对正负链上5’和3’区域的序列进行标注
####### Positive labelled, 5' region
promoter_positiveLabel <- rbind(promoter_left_strand, promoter_right_strand)
promoter_positiveLabel$strand <- rep('+', times = nrow(promoter_positiveLabel))
promoter_positiveLabel$transcript_ID[1:nrow(promoter_left_strand)]<- paste0(promoter_left_strand$geneid,'-1')
promoter_positiveLabel$transcript_ID[(nrow(promoter_left_strand)+1):nrow(promoter_positiveLabel)] <- paste0(promoter_left_strand$geneid,'-2')

###### Negative labelled, 3' region
promoter_negativeLabel <- rbind(promoter_left_strand, promoter_right_strand)
promoter_negativeLabel$strand <- rep('-', times = nrow(promoter_negativeLabel))
promoter_negativeLabel$transcript_ID[1:nrow(promoter_left_strand)]<- paste0(promoter_left_strand$geneid,'-3')
promoter_negativeLabel$transcript_ID[(nrow(promoter_left_strand)+1):nrow(promoter_negativeLabel)] <- paste0(promoter_left_strand$geneid,'-4')

###### Merge 5' region and 3' region
promoter <- rbind(promoter_positiveLabel, promoter_negativeLabel)
######


promoter
for (i in 1:nrow(promoter)){## here we found that the 'chr' does not have 'chr' label, we need to add it
  promoter$Chr[i] <- paste0('chr',promoter$chr[i])
}
promoter <- data.frame(transcript_ID = promoter$transcript_ID, chr = promoter$Chr, start = promoter$Start, end = promoter$End, strand = promoter$strand) 
head(promoter); dim(promoter)
promoter; class(promoter)
write.table(x = promoter, file = '/Users/mijiarui/biosoft/HOMER/results/new/promoter_significant_beta_vs_alphaDelta_knownLncRNA.txt', 
            sep = '\t', row.names = F, col.names = F, quote = F)





#======================================================================================
#             alpha-cell vs beta-cell and alpha-cell vs delta-cell
#======================================================================================

### Differential gene expression between sampleC (alpha-cell) and sampleA (beta-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleC, sampleA)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC


#### Calculate the mean value of gene expression in sampleC 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanC, baseMeanA, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)


#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleC,sampleA,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_alpha_vs_beta <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1))

length(signifiant_genes_alpha_vs_beta)





### Differential gene expression between sampleA (beta-cell) and sampleD (delta-cell)
#### Comparison between sampleA and sampleD
contrastV <- c('celltype', sampleC, sampleD)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC


#### Calculate the mean value of gene expression in sampleC 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanC, baseMeanD, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)


#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleC,sampleD,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_alpha_vs_delta <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1 ))

length(signifiant_genes_alpha_vs_delta)


#### Pick up the lincRNAs that are differentially expressed in between beta-cell and alpha-cell and
#### between beta-cell and delta-cell
a <- signifiant_genes_alpha_vs_delta[signifiant_genes_alpha_vs_delta %in% signifiant_genes_alpha_vs_beta]
length(a)

coordinate_a_lincRNA <- lincRNA[lincRNA$ensembl %in% a,]
coordinate_a_lincRNA


coordinate_a_lincRNA$chr <- as.numeric(sub("chr","",coordinate_a_lincRNA$chr))


library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_a_lincRNA)){
  coordinate_a_lincRNA$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(coordinate_a_lincRNA$chr[i],coordinate_a_lincRNA$start[i]-1000000,coordinate_a_lincRNA$end[i]+1000000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_a_lincRNA$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)

#### Integrate lincRNA and antisense lncRNA in alpha-cell

coordinate_a <- rbind(coordinate_a_antisense, coordinate_a_lincRNA)
dim(coordinate_a)

for (i in 1:nrow(coordinate_a)){
  coordinate_a$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(coordinate_a$chr[i],coordinate_a$start[i]-1000000,coordinate_a$end[i]+1000000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_a$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)

kk <- enrichKEGG(entrez,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)

pa <- enrichPathway(entrez,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
plot(pa)


################## Pick up the promoter regions ###################
### Normally we define the promoter regions as 500 bp upstream of TSS. Because this is unstranded library, we need to check
### the both ends. It is good to use 'dplyr' package here.
library('dplyr')
test_beta_coordinate <- coordinate_a
promoter_left_strand <- mutate(.data = test_beta_coordinate, Start = start -100, End = start-90)[, c(1,8,9,5)]
promoter_left_strand
promoter_right_strand <- mutate(.data = test_beta_coordinate, Start = end + 90, End = end + 100 )[,c(1,8,9,5)]
promoter_right_strand

######
####### ++/--
####### 对正负链上5’和3’区域的序列进行标注
####### Positive labelled, 5' region
promoter_positiveLabel <- rbind(promoter_left_strand, promoter_right_strand)
promoter_positiveLabel$strand <- rep('+', times = nrow(promoter_positiveLabel))
promoter_positiveLabel$transcript_ID[1:nrow(promoter_left_strand)]<- paste0(promoter_left_strand$geneid,'-1')
promoter_positiveLabel$transcript_ID[(nrow(promoter_left_strand)+1):nrow(promoter_positiveLabel)] <- paste0(promoter_left_strand$geneid,'-2')

###### Negative labelled, 3' region
promoter_negativeLabel <- rbind(promoter_left_strand, promoter_right_strand)
promoter_negativeLabel$strand <- rep('-', times = nrow(promoter_negativeLabel))
promoter_negativeLabel$transcript_ID[1:nrow(promoter_left_strand)]<- paste0(promoter_left_strand$geneid,'-3')
promoter_negativeLabel$transcript_ID[(nrow(promoter_left_strand)+1):nrow(promoter_negativeLabel)] <- paste0(promoter_left_strand$geneid,'-4')

###### Merge 5' region and 3' region
promoter <- rbind(promoter_positiveLabel, promoter_negativeLabel)
######


promoter
for (i in 1:nrow(promoter)){## here we found that the 'chr' does not have 'chr' label, we need to add it
  promoter$Chr[i] <- paste0('chr',promoter$chr[i])
}
promoter <- data.frame(transcript_ID = promoter$transcript_ID, chr = promoter$Chr, start = promoter$Start, end = promoter$End, strand = promoter$strand) 
head(promoter); dim(promoter)
promoter; class(promoter)
write.table(x = promoter, file = '/Users/mijiarui/biosoft/HOMER/results/new/promoter_significant_alpha_vs_betaDelta_knownLncRNA.txt', 
            sep = '\t', row.names = F, col.names = F, quote = F)






#======================================================================================
#             delta-cell vs beta-cell and delta-cell vs alpha-cell
#======================================================================================

### Differential gene expression between sampleD (delta-cell) and sampleA (beta-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleD, sampleA)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD


#### Calculate the mean value of gene expression in sampleC 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanD, baseMeanA, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)



#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleD,sampleA,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_delta_vs_beta <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1))

length(signifiant_genes_delta_vs_beta)


### Differential gene expression between sampleD (delta-cell) and sampleC (alpha-cell)
#### Comparison between sampleA and sampleD
contrastV <- c('celltype', sampleD, sampleC)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleC 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD

#### Calculate the mean value of gene expression in sampleA 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC



#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanD, baseMeanC, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)





#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleD,sampleC,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_delta_vs_alpha <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1 ))

length(signifiant_genes_delta_vs_alpha)


#### Pick up the lincRNAs that are differentially expressed in between beta-cell and alpha-cell and
#### between beta-cell and delta-cell
d <- signifiant_genes_delta_vs_beta[signifiant_genes_delta_vs_beta %in% signifiant_genes_delta_vs_alpha]
length(d)

coordinate_d_lincRNA <- lincRNA[lincRNA$ensembl %in% d,]
coordinate_d_lincRNA

coordinate_d_lincRNA$chr <- as.numeric(sub("chr","",coordinate_d_lincRNA$chr))


library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_d_lincRNA)){
  coordinate_d_lincRNA$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(coordinate_d_lincRNA$chr[i],coordinate_d_lincRNA$start[i]-400000,coordinate_d_lincRNA$end[i]+400000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_d_lincRNA$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)


#### Integrate lincRNA and antisense lncRNA in delta-cell

coordinate_d <- rbind(coordinate_a_antisense, coordinate_d_lincRNA)

for (i in 1:nrow(coordinate_d)){
  coordinate_d$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(coordinate_d$chr[i],coordinate_d$start[i]-1000000,coordinate_d$end[i]+1000000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_d$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)


kk <- enrichKEGG(entrez,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)

pa <- enrichPathway(entrez,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
barplot(pa, showCategory = 20)
cnetplot(pa)



################## Pick up the promoter regions ###################
### Normally we define the promoter regions as 500 bp upstream of TSS. Because this is unstranded library, we need to check
### the both ends. It is good to use 'dplyr' package here.
library('dplyr')
test_beta_coordinate <- coordinate_d
promoter_left_strand <- mutate(.data = test_beta_coordinate, Start = start -100, End = start-90)[, c(1,8,9,5)]
promoter_left_strand
promoter_right_strand <- mutate(.data = test_beta_coordinate, Start = end + 90, End = end + 100 )[,c(1,8,9,5)]
promoter_right_strand

######
####### ++/--
####### 对正负链上5’和3’区域的序列进行标注
####### Positive labelled, 5' region
promoter_positiveLabel <- rbind(promoter_left_strand, promoter_right_strand)
promoter_positiveLabel$strand <- rep('+', times = nrow(promoter_positiveLabel))
promoter_positiveLabel$transcript_ID[1:nrow(promoter_left_strand)]<- paste0(promoter_left_strand$geneid,'-1')
promoter_positiveLabel$transcript_ID[(nrow(promoter_left_strand)+1):nrow(promoter_positiveLabel)] <- paste0(promoter_left_strand$geneid,'-2')

###### Negative labelled, 3' region
promoter_negativeLabel <- rbind(promoter_left_strand, promoter_right_strand)
promoter_negativeLabel$strand <- rep('-', times = nrow(promoter_negativeLabel))
promoter_negativeLabel$transcript_ID[1:nrow(promoter_left_strand)]<- paste0(promoter_left_strand$geneid,'-3')
promoter_negativeLabel$transcript_ID[(nrow(promoter_left_strand)+1):nrow(promoter_negativeLabel)] <- paste0(promoter_left_strand$geneid,'-4')

###### Merge 5' region and 3' region
promoter <- rbind(promoter_positiveLabel, promoter_negativeLabel)
######


promoter
for (i in 1:nrow(promoter)){## here we found that the 'chr' does not have 'chr' label, we need to add it
  promoter$Chr[i] <- paste0('chr',promoter$chr[i])
}
promoter <- data.frame(transcript_ID = promoter$transcript_ID, chr = promoter$Chr, start = promoter$Start, end = promoter$End, strand = promoter$strand) 
head(promoter); dim(promoter)
promoter; class(promoter)
write.table(x = promoter, file = '/Users/mijiarui/biosoft/HOMER/results/new/promoter_significant_delta_vs_alphaBeta_knownLncRNA.txt', 
            sep = '\t', row.names = F, col.names = F, quote = F)





#### heatmap containing differentially expressed lincRNA and antisense lncRNA
abd_ensembl <- c(a,b,d)
abd_coordinate <- rbind(coordinate_a, coordinate_b, coordinate_d)

pheatmap(rlogMat[rownames(rlogMat) %in% abd_ensembl,], color = colorRampPalette(c('navy', 'white', 'firebrick3'))(50), 
         cluster_rows = T, scale = 'row', cluster_cols = T, annotation_col = sample)

abd_coordinate$chr <- as.numeric(sub("chr","",abd_coordinate$chr))

library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(abd_coordinate)){
  abd_coordinate$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(abd_coordinate$chr[i],abd_coordinate$start[i]-1000000,abd_coordinate$end[i]+1000000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(abd_coordinate$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)


#===================================================================================================
#
#                          Differential expressed annotated lincRNA
#
#===================================================================================================

# Make DESeq object and Do the normalization (size factor)
## Create DESeqDataSet
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = lincRNA_exprs, colData = sample, design  = ~ celltype)

## Calculation and Normalization: get normalized_counts
dds <- DESeq(ddsFullCountTable)
normalized_counts <- counts(dds, normalized= T)
class(normalized_counts); head(normalized_counts)

## Sort according to mad
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts < normalized_counts[order(normalized_counts_mad, decreasing = T),]
## Log transformation
rld <- rlog(dds, blind = F) # This step takes some time, we do not put '+1' here
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing = T),]


#======================================================================================
#             alpha-cell vs beta-cell and delta-cell vs beta-cell
#======================================================================================

################### Perform differential gene expression ################
############### Differential Gene Expression #################
sampleA <- 'beta'
sampleB <- 'acinal'
sampleC <- 'alpha'
sampleD <- 'delta'

### Differential gene expression between sampleC (alpha-cell) and sampleA (beta-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleC, sampleA)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC


#### Calculate the mean value of gene expression in sampleC 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanC, baseMeanA, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)


#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleC,sampleA,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_alpha_vs_beta <- row.names(subset(res, padj < 0.01 & res$log2FoldChange > 2))

length(signifiant_genes_alpha_vs_beta)


### Differential gene expression between sampleD (delta-cell) and sampleA (beta-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleD, sampleA)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD


#### Calculate the mean value of gene expression in sampleC 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanD, baseMeanA, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)



#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleD,sampleA,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_delta_vs_beta <- row.names(subset(res, padj < 0.01 & res$log2FoldChange > 1))

length(signifiant_genes_delta_vs_beta)



b <- signifiant_genes_alpha_vs_beta[signifiant_genes_alpha_vs_beta %in% signifiant_genes_delta_vs_beta]
length(b)

coordinate_b <- lincRNA[lincRNA$ensembl %in% b,]
coordinate_b

coordinate_b$chr <- as.numeric(sub("chr","",coordinate_b$chr))


library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_b)){
  coordinate_b$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(coordinate_b$chr[i],coordinate_b$start[i]-10000,coordinate_b$end[i]+10000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_b$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)


#======================================================================================
#             beta-cell vs alpha-cell and delta-cell vs alpha-cell
#======================================================================================

### Differential gene expression between sampleA (beta-cell) and sampleC (alpha-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleA, sampleC)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA


#### Calculate the mean value of gene expression in sampleC 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanA, baseMeanC, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)



#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleC,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


signifiant_genes_beta_vs_alpha <- row.names(subset(res, padj < 0.01 & res$log2FoldChange > 1 ))

length(signifiant_genes_beta_vs_alpha)



### Differential gene expression between sampleD (delta-cell) and sampleC (alpha-cell)
#### Comparison between sampleA and sampleD
contrastV <- c('celltype', sampleD, sampleC)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleC 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD

#### Calculate the mean value of gene expression in sampleA 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC



#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanD, baseMeanC, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)





#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleD,sampleC,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_delta_vs_alpha <- row.names(subset(res, padj < 0.01 & res$log2FoldChange > 1 ))

length(signifiant_genes_delta_vs_alpha)


#### Pick up the lincRNAs that are differentially expressed in between beta-cell and alpha-cell and
#### between beta-cell and delta-cell
a <- signifiant_genes_beta_vs_alpha[signifiant_genes_beta_vs_alpha %in% signifiant_genes_delta_vs_alpha]
length(a)

coordinate_a <- lincRNA[lincRNA$ensembl %in% a,]
coordinate_a


coordinate_a$chr <- as.numeric(sub("chr","",coordinate_a$chr))


library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_a)){
  coordinate_a$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(coordinate_a$chr[i],coordinate_a$start[i]-10000,coordinate_a$end[i]+10000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_a$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)




#======================================================================================
#             alpha-cell vs delta-cell and beta-cell vs delta-cell
#======================================================================================

### Differential gene expression between sampleA (beta-cell) and sampleD (delta-cell)
#### Comparison between sampleA and sampleD
contrastV <- c('celltype', sampleC, sampleD)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC


#### Calculate the mean value of gene expression in sampleC 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanC, baseMeanD, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)


#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleC,sampleD,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)



signifiant_genes_alpha_vs_delta <- row.names(subset(res, padj < 0.01 & res$log2FoldChange > 2 ))

length(signifiant_genes_alpha_vs_delta)


### Differential gene expression between sampleA (beta-cell) and sampleD (delta-cell)
#### Comparison between sampleA and sampleD
contrastV <- c('celltype', sampleA, sampleD)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA


#### Calculate the mean value of gene expression in sampleC 
baseD <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleD]
baseMeanD <- as.data.frame(rowMeans(baseD, na.rm = T))
colnames(baseMeanD) <- sampleD

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanA, baseMeanD, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)


#### Extract differential expressed gene
colnames(res)
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleD,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


signifiant_genes_beta_vs_delta <- row.names(subset(res, padj < 0.01 & res$log2FoldChange > 1))

length(signifiant_genes_beta_vs_delta)


#### Pick up the lincRNAs that are differentially expressed in between beta-cell and alpha-cell and
#### between beta-cell and delta-cell
b <- signifiant_genes_alpha_vs_delta[signifiant_genes_alpha_vs_delta %in% signifiant_genes_beta_vs_delta]
length(b)

coordinate_b <- lincRNA[lincRNA$ensembl %in% b,]
coordinate_b

coordinate_b$chr <- as.numeric(sub("chr","",coordinate_b$chr))


library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_b)){
  coordinate_b$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                    values= list(coordinate_b$chr[i],coordinate_b$start[i]-100000,coordinate_b$end[i]+100000), 
                                    mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_b$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)










#===================================================================================================
#
#      Differential expressed annotated antisense lncRNA (between endocrine and exocrine)
#
#===================================================================================================
setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/lncRNA_EDA')
reads <- read.csv('counts_total.csv', header = T, row.names = 1)
reads[1:3,1:3];dim(reads)
## In total, 23 samples with 32266 genes detected

## Read in Coldata(Phenodata)
sample <- read.csv('SampleGroup_endocrine_vs_exocrine.csv', header = T, row.names = 1, colClasses = 'factor')
sample

## Read in metadata of lincRNA and antisense (lncRNA)
setwd("/Users/mijiarui/RNA-seq/CAGE")
lincRNA <- read.table(file = "annotated.lincRNA.bed.txt", header = F, sep = '\t', quote = "")
colnames(lincRNA) <- c('chr', 'start', 'end', 'ensembl', 'geneid', 'strand')
row.names(lincRNA) <- lincRNA$ensembl
head(lincRNA); dim(lincRNA)
antisense <- read.table(file = "annotated.antisense.bed.txt", header = F, sep = '\t', quote = "")
colnames(antisense) <- c('chr', 'start', 'end', 'ensembl', 'geneid', 'strand')
row.names(antisense) <- antisense$ensembl
head(antisense); dim(antisense)

### expression matrix of lincRNA and antisense (lncRNA)
lincRNA_exprs <- reads[row.names(reads) %in% row.names(lincRNA),]
head(lincRNA_exprs); dim(lincRNA_exprs)

antisense_exprs <- reads[row.names(reads) %in% row.names(antisense),]
head(antisense_exprs); dim(antisense_exprs)



# Make DESeq object and Do the normalization (size factor)
## Create DESeqDataSet
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = antisense_exprs, colData = sample, design  = ~ celltype)

## Calculation and Normalization: get normalized_counts
dds <- DESeq(ddsFullCountTable)
normalized_counts <- counts(dds, normalized= T)
class(normalized_counts); head(normalized_counts)

## Sort according to mad
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts < normalized_counts[order(normalized_counts_mad, decreasing = T),]
## Log transformation
rld <- rlog(dds, blind = F) # This step takes some time, we do not put '+1' here
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing = T),]




#======================================================================================
#                          endocrine higher than exocrine
#======================================================================================

################### Perform differential gene expression ################
############### Differential Gene Expression #################
sampleA <- 'endocrine'
sampleC <- 'exocrine'



### Differential gene expression between sampleA (beta-cell) and sampleC (alpha-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleA, sampleC)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA


#### Calculate the mean value of gene expression in sampleC 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanA, baseMeanC, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)



#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleC,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


signifiant_genes_endocrine_vs_exocrine <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1 ))

length(signifiant_genes_endocrine_vs_exocrine)


coordinate_antisense <- antisense[antisense$ensembl %in% signifiant_genes_endocrine_vs_exocrine,]
coordinate_antisense
dim(coordinate_antisense)

coordinate_antisense$chr <- as.numeric(sub("chr","",coordinate_antisense$chr))


library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_antisense)){
  coordinate_antisense$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                              values= list(coordinate_antisense$chr[i],coordinate_antisense$start[i]-500000,coordinate_antisense$end[i]+500000), 
                                              mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_antisense$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)

pa <- enrichPathway(entrez,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
barplot(pa)
cnetplot(pa)

kk <- enrichKEGG(entrez,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)

endocrine_highly_exprs <- rlogMat[rownames(rlogMat) %in% row.names(coordinate_antisense),]
dim(endocrine_highly_exprs)
endocrine_highly_exprs_var <- apply(endocrine_highly_exprs[,2:15], 1, var)
endocrine_highly_exprs1 <- endocrine_highly_exprs[order(endocrine_highly_exprs_var, decreasing = T),]
pheatmap(endocrine_highly_exprs1, color = colorRampPalette(c('navy', 'white', 'firebrick3'))(50), 
         cluster_rows = T, scale = 'row', cluster_cols = T, annotation_col = sample)

### add the flanking genes will alter the structure of the table, therefore, we redo the 赋值 step
coordinate_antisense <- antisense[antisense$ensembl %in% signifiant_genes_endocrine_vs_exocrine,]
write.table(x = coordinate_antisense, file = 'coordinate_antisense_endocrine_higherThan_exocrine.txt', sep = '\t', quote = F)

#======================================================================================
#                          exocrine higher than endocrine
#======================================================================================

################### Perform differential gene expression ################
############### Differential Gene Expression #################
sampleA <- 'exocrine'
sampleC <- 'endocrine'



### Differential gene expression between sampleA (beta-cell) and sampleC (alpha-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleA, sampleC)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA


#### Calculate the mean value of gene expression in sampleC 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanA, baseMeanC, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)



#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleC,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


signifiant_genes_exocrine_vs_endocrine <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1 ))

length(signifiant_genes_exocrine_vs_endocrine)


coordinate_antisense <- antisense[antisense$ensembl %in% signifiant_genes_exocrine_vs_endocrine,]
coordinate_antisense
dim(coordinate_antisense)


coordinate_antisense$chr <- as.numeric(sub("chr","",coordinate_antisense$chr))
write.table(x = coordinate_antisense, file = 'coordinate_antisense_exocrine_higherThan_endocrine.txt', quote = F)

library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_antisense)){
  coordinate_antisense$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                            values= list(coordinate_antisense$chr[i],coordinate_antisense$start[i]-500000,coordinate_antisense$end[i]+500000), 
                                            mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_antisense$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)

pa <- enrichPathway(entrez,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
barplot(pa)

kk <- enrichKEGG(entrez,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)

exocrine_highly_exprs <- rlogMat[rownames(rlogMat) %in% signifiant_genes_exocrine_vs_endocrine,]
dim(exocrine_highly_exprs)
exocrine_highly_exprs_var <- apply(exocrine_highly_exprs[,2:15], 1, var)
exocrine_highly_exprs1 <- endocrine_highly_exprs[order(exocrine_highly_exprs_var, decreasing = T),]
pheatmap(exocrine_highly_exprs1, color = colorRampPalette(c('navy', 'white', 'firebrick3'))(50), 
         cluster_rows = T, scale = 'row', cluster_cols = T, annotation_col = sample)
















#===================================================================================================
#
#      Differential expressed annotated lincRNA (between endocrine and exocrine)
#
#===================================================================================================
setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/lncRNA_EDA')
reads <- read.csv('counts_total.csv', header = T, row.names = 1)
reads[1:3,1:3];dim(reads)
## In total, 23 samples with 32266 genes detected

## Read in Coldata(Phenodata)
sample <- read.csv('SampleGroup_endocrine_vs_exocrine.csv', header = T, row.names = 1, colClasses = 'factor')
sample

## Read in metadata of lincRNA and antisense (lncRNA)
setwd("/Users/mijiarui/RNA-seq/CAGE")
lincRNA <- read.table(file = "annotated.lincRNA.bed.txt", header = F, sep = '\t', quote = "")
colnames(lincRNA) <- c('chr', 'start', 'end', 'ensembl', 'geneid', 'strand')
row.names(lincRNA) <- lincRNA$ensembl
head(lincRNA); dim(lincRNA)
antisense <- read.table(file = "annotated.antisense.bed.txt", header = F, sep = '\t', quote = "")
colnames(antisense) <- c('chr', 'start', 'end', 'ensembl', 'geneid', 'strand')
row.names(antisense) <- antisense$ensembl
head(antisense); dim(antisense)

### expression matrix of lincRNA and antisense (lncRNA)
lincRNA_exprs <- reads[row.names(reads) %in% row.names(lincRNA),]
head(lincRNA_exprs); dim(lincRNA_exprs)

antisense_exprs <- reads[row.names(reads) %in% row.names(antisense),]
head(antisense_exprs); dim(antisense_exprs)



# Make DESeq object and Do the normalization (size factor)
## Create DESeqDataSet
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = lincRNA_exprs, colData = sample, design  = ~ celltype)

## Calculation and Normalization: get normalized_counts
dds <- DESeq(ddsFullCountTable)
normalized_counts <- counts(dds, normalized= T)
class(normalized_counts); head(normalized_counts)

## Sort according to mad
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts < normalized_counts[order(normalized_counts_mad, decreasing = T),]
## Log transformation
rld <- rlog(dds, blind = F) # This step takes some time, we do not put '+1' here
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing = T),]




#======================================================================================
#                          endocrine higher than exocrine
#======================================================================================

################### Perform differential gene expression ################
############### Differential Gene Expression #################
sampleA <- 'endocrine'
sampleC <- 'exocrine'



### Differential gene expression between sampleA (beta-cell) and sampleC (alpha-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleA, sampleC)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA


#### Calculate the mean value of gene expression in sampleC 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanA, baseMeanC, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)



#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleC,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


signifiant_genes_endocrine_vs_exocrine <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1 ))

length(signifiant_genes_endocrine_vs_exocrine)


coordinate_lincRNA <- lincRNA[lincRNA$ensembl %in% signifiant_genes_endocrine_vs_exocrine,]
coordinate_lincRNA
dim(coordinate_lincRNA)

coordinate_lincRNA$chr <- as.numeric(sub("chr","",coordinate_lincRNA$chr))


library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_lincRNA)){
  coordinate_lincRNA$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                              values= list(coordinate_lincRNA$chr[i],coordinate_lincRNA$start[i]-500000,coordinate_lincRNA$end[i]+500000), 
                                              mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_lincRNA$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
barplot(GO)
cnetplot(GO)


kk <- enrichKEGG(entrez,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)


pa <- enrichPathway(entrez,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
barplot(pa)
cnetplot(pa)


endocrine_highly_exprs <- rlogMat[rownames(rlogMat) %in% row.names(coordinate_lincRNA),]
dim(endocrine_highly_exprs)
endocrine_highly_exprs_var <- apply(endocrine_highly_exprs[,2:15], 1, var)
endocrine_highly_exprs1 <- endocrine_highly_exprs[order(endocrine_highly_exprs_var, decreasing = T),]
pheatmap(endocrine_highly_exprs1, color = colorRampPalette(c('navy', 'white', 'firebrick3'))(50), 
         cluster_rows = T, scale = 'row', cluster_cols = T, annotation_col = sample)


### add the flanking genes will alter the structure of the table, therefore, we redo the 赋值 step
coordinate_lincRNA <- lincRNA[lincRNA$ensembl %in% signifiant_genes_endocrine_vs_exocrine,]
write.table(x = coordinate_lincRNA, file = 'coordinate_lincRNA_endocrine_higherThan_exocrine.txt', sep = '\t', quote = F)


#======================================================================================
#                          exocrine higher than endocrine
#======================================================================================

################### Perform differential gene expression ################
############### Differential Gene Expression #################
sampleA <- 'exocrine'
sampleC <- 'endocrine'



### Differential gene expression between sampleA (beta-cell) and sampleC (alpha-cell)
#### Comparison between sampleA and sampleC
contrastV <- c('celltype', sampleA, sampleC)
res <- results(dds, contrast = contrastV)
head(res)


#### Calculate the mean value of gene expression in sampleA 
baseA <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleA]
baseMeanA <- as.data.frame(rowMeans(baseA, na.rm = T))
colnames(baseMeanA) <- sampleA


#### Calculate the mean value of gene expression in sampleC 
baseC <- counts(dds, normalized = T)[,colData(dds)$celltype==sampleC]
baseMeanC <- as.data.frame(rowMeans(baseC, na.rm = T))
colnames(baseMeanC) <- sampleC

#### Prepare the res table of differential gene expression between sampleA and sampleC
res <- cbind(baseMeanA, baseMeanC, as.data.frame(res))
ID <- rownames(res)
res <- cbind(ID, as.data.frame(res))
res$padj[is.na(res$padj)]<- 1
res <- res[order(res$padj),]
res$significance <- (res$padj<0.05)



#### Extract differential expressed gene
res_de <- subset(res, res$padj<0.1, select=c('ID',sampleA,sampleC,'log2FoldChange','padj'))
res_de_up <- subset(res_de,res_de$log2FoldChange>=1)
res_de_down <- subset(res_de,res_de$log2FoldChange<=-1)


signifiant_genes_exocrine_vs_endocrine <- row.names(subset(res, padj < 0.05 & res$log2FoldChange > 1 ))

length(signifiant_genes_exocrine_vs_endocrine)


coordinate_lincRNA <- lincRNA[lincRNA$ensembl %in% signifiant_genes_exocrine_vs_endocrine,]
coordinate_lincRNA
dim(coordinate_lincRNA)

coordinate_lincRNA$chr <- as.numeric(sub("chr","",coordinate_lincRNA$chr))
write.table(x = coordinate_lincRNA, file = 'coordinate_lincRNA_exocrine_higherThan_endocrine.txt', quote = F)

library(biomaRt)
library(clusterProfiler)
readable = T
mart <- useMart('ensembl')
ensembl <- useDataset('drerio_gene_ensembl', mart) ## # select dataset "drerio_gene_ensembl"
for (i in 1:nrow(coordinate_lincRNA)){
  coordinate_lincRNA$flanking[i] <- getBM(attributes = 'entrezgene', filters = c('chromosome_name','start','end'),
                                            values= list(coordinate_lincRNA$chr[i],coordinate_lincRNA$start[i]-200000,coordinate_lincRNA$end[i]+200000), 
                                            mart=ensembl)
}
entrez <- as.vector(unlist(coordinate_lincRNA$flanking))
GO <- enrichGO(entrez,'org.Dr.eg.db',pvalueCutoff = 0.2,
               pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
GO
dotplot(GO)
cnetplot(GO)

pa <- enrichPathway(entrez,organism = 'zebrafish',pvalueCutoff = 0.05 ,pAdjustMethod = 'BH',qvalueCutoff = 0.2, readable = T)
summary(pa)
barplot(pa)
cnetplot(pa)

kk <- enrichKEGG(entrez,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)

exocrine_highly_exprs <- rlogMat[rownames(rlogMat) %in% signifiant_genes_exocrine_vs_endocrine,]
dim(exocrine_highly_exprs)
exocrine_highly_exprs_var <- apply(exocrine_highly_exprs[,c(1,16:23)], 1, var)
exocrine_highly_exprs1 <- exocrine_highly_exprs[order(exocrine_highly_exprs_var, decreasing = T),]
pheatmap(exocrine_highly_exprs1, color = colorRampPalette(c('navy', 'white', 'firebrick3'))(50), 
         cluster_rows = T, scale = 'row', cluster_cols = T, annotation_col = sample)


