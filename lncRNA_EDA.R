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
library("amap")
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
library("BiocParallel")
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
rld <- rlog(dds, blind = F) # This step takes some time
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing = T),]

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
signifiant_genes <- getBM(attribute=c('ensembl_gene_id', 'entrezgene','zfin_id_symbol'),filters = 'ensembl_gene_id', values= row.names(subset(res, -log10(padj) > 100)),mart = ensembl)

##### Plot VolcanoPlot
this_title <- paste('Cutoff for logFC is 1', '\nThe number of upregulated genes is', nrow(res[res$change == 'UP',]),
                     '\nThe number of downregulated gene is', nrow(res[res$change == 'DOWN',]))
ggplot(data = res, aes(x=log2FoldChange, y = -log10(padj), color = change, alpha = 0.5)) +
  geom_point()+
  scale_color_manual(values = c('blue','black', 'red')) +
  geom_hline(yintercept = -log10(0.05),lty=4, lwd=0.6,alpha=0.5)+
  geom_vline(xintercept = c(-1,1),lty=4, lwd=0.6,alpha=0.5)+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title= this_title, x= 'log2(fold change)', y = '-log10(padj)')+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text_repel(data=subset(res, -log10(padj) > 100), aes(label=signifiant_genes$zfin_id_symbol), col= 'black')



# Heatmap of beta cell and acinal cell
res_de_up_top50_id <- as.vector(head(res_de_up$ID,50))
res_de_down_top50_id <- as.vector(head(res_de_down$ID,50))
res_de_top100 <- c(res_de_up_top50_id,res_de_down_top50_id)
res_de_top100_expr <- normalized_counts[rownames(normalized_counts) %in% res_de_top100,c(2:6,1,16:20)]
head(res_de_top100_expr)
pheatmap(res_de_top100_expr,cluster_rows = T, scale = 'row',annotation_col = sample)



# Use biomaRt to concert gene ID (for enrichment analysis, entrezid should be used)
head(listMarts())
mart <- useMart('ensembl')
View(listDatasets(mart))  # select dataset "drerio_gene_ensembl"
ensembl <- useDataset('drerio_gene_ensembl', mart)
## beta higher than acinal
entrezid_beta_higherthan_acinal<- getBM(attribute=c('ensembl_gene_id', 'entrezgene'),filters = 'ensembl_gene_id', values= row.names(res_de_up),mart = ensembl)
head(entrezid_beta_higherthan_acinal)
dim(entrezid_beta_higherthan_acinal)
entrezid_beta_higherthan_acinal <- entrezid_beta_higherthan_acinal[is.na(entrezid_beta_higherthan_acinal$entrezgene) != T,]
summary(entrezid_beta_higherthan_acinal)
dim(res_de)
head(res_de)
geneList_foldchange_beta_higherthan_acinal <- res_de[res_de$ID %in% entrezid_beta_higherthan_acinal$ensembl_gene_id,c('ID','log2FoldChange')]
head(geneList_foldchange_beta_higherthan_acinal)
entrezid_beta_higherthan_acinal <- unique(entrezid_beta_higherthan_acinal$entrezgene)
### Biological Process
BP <- enrichGO(entrezid_beta_higherthan_acinal,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'BP', readable = readable)
result_BP <- simplify(BP,cutoff = 0.7, by = 'p.adjust', select_fun = min)
head(result_BP)
?dotplot
dotplot(result_BP,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.04)+scale_color_continuous(low = 'purple', high = 'green')  # Use scale_size to change the size of the bubble, if the bubble is too large and be cut by the edge, you can used xlim to extend the x axis
dotplot(result_BP,showCategory = 10,x = 'count')+scale_size(range = c(2,12))+xlim(NA,120)
barplot(result_BP,drop = T, showCategory = 10)+scale_x_discrete(labels = function(x) str_wrap(x, width = 25))
barplot(result_BP,x = 'count', showCategory = 10)
plotGOgraph(result_BP)
### Molecular Function
MF <- enrichGO(entrezid_beta_higherthan_acinal,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'MF', readable = readable)
result_MF <- simplify(MF,cutoff = 0.7, by = 'p.adjust', select_fun = min)
dotplot(result_BP,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.04)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_MF,drop = T, showCategory = 10)+scale_x_discrete(labels = function(x) str_wrap(x, width = 25))
head(result_MF)
plotGOgraph(result_MF)
### Cellular Component
CC <- enrichGO(entrezid_beta_higherthan_acinal,'org.Dr.eg.db',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2,ont = 'CC', readable = readable)
result_CC <- simplify(CC,cutoff = 0.7, by = 'p.adjust', select_fun = min)
dotplot(result_CC,showCategory = 10)+scale_size(range = c(2,15))+ggplot2::xlim(NA, 0.045)+scale_color_continuous(low = 'purple', high = 'green') 
barplot(result_CC,drop = T, showCategory = 10)+scale_x_discrete(labels = function(x) str_wrap(x, width = 25))
head(result_CC)
plotGOgraph(result_CC)
### KEGG pathway
kk <- enrichKEGG(entrezid_beta_higherthan_acinal,organism = 'dre',keyType = 'kegg',pvalueCutoff = 0.05,pAdjustMethod = 'BH',qvalueCutoff = 0.2)
result_kk <- setReadable(kk,'org.Dr.eg.db',keytype = 'ENTREZID')
head(result_kk)
dotplot(result_kk)+scale_size(range = c(2,15))+ggplot2::xlim(NA,0.09)
cnetplot(result_kk,categorySize = 'aa$geneNum', showCategory = 3)
?cnetplot
result_kk$Description
head(geneList_foldchange_beta_higherthan_acinal)
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



