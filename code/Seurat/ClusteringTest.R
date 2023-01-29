
# Step 1 Setup the Seruat Object

#  Load library 

library(dplyr) # Load dplyr package

library(Seurat) # Load seurat package


# intial file path 
scRNAexpressionfilepath = "D:/data/GSE115469/matrix_mapped.txt" # Input 
markergenepath = "D:/data/GSE115469/bioMarkers.txt"  # File that Seurat will write biomakers 




# load our matrix data 
counts<-read.table(file = scRNAexpressionfilepath, header = TRUE, sep = "\t",row.names = 1)

# Create seurat object 
pbmc <- CreateSeuratObject(counts = counts, project = "hl_test", min.cells = 3, min.features = 200)
pbmc

#############################################################################################################################################

# Step 2 Pre-processing workflow 

# Normalizing the data 

# Default normalization method is "LogNormalize". 
# Options: LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. 
#                        This is then natural-log transformed using log1p.
#
#          CLR: Applies a centered log ratio transformation
#
#          RC: Relative counts. Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. 
#              No log-transformation is applied. For counts per million (CPM) set

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize") 



# Identification of highly variable features (feature selection)

# Options for the selection.method:
#  "vst": First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). 
#        Then standardizes the feature values using the observed mean and expected variance (given by the fitted line).
#        Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter). 
#
#  "mean.var.plot": First, uses a function to calculate average expression (mean.function) and dispersion (dispersion.function) for each feature.
#                  Next, divides features into num.bin (deafult 20) bins based on their average expression, 
#                  and calculates z-scores for dispersion within each bin.
#                  The purpose of this is to identify variable features 
#                  while controlling for the strong relationship between variability and average expression.
#
#  "dispersion": selects the genes with the highest dispersion values
#
# nfeatures: Number of features to select as top variable features; 
#           only used when selection.method is set to 'dispersion' or 'vst'
#
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",nfeatures = 2000)
pbmc

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2




# Scaling the data

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)






# Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Dimensional reduction results visualization 
DimPlot(pbmc, reduction = "pca")




# Elbowplot

ElbowPlot(pbmc)

# Cluster the cells

pbmc <- FindNeighbors(pbmc, dims = 1:15)

pbmc <- FindClusters(pbmc, resolution = 0.5)


DimPlot(pbmc, reduction = "pca")






# Visualization of clustering
pbmc <- RunTSNE(pbmc, dims = 1:15)

DimPlot(pbmc, reduction = "tsne")

#pbmc <- RunUMAP(pbmc, dims = 1:15)

#DimPlot(pbmc, reduction = "umap")


markergenepath = "D:/data/GSE115469/bioMarkers.txt"  # File that Seurat will write biomakers 

# Finding differentially expressed features (cluster biomarkers)
pbmc.markers <- FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 2)

# Write the expressed features into a file
write.table(pbmc.markers,markergenepath, row.names = FALSE,sep = "\t")


#violin plot for the biomarker

VlnPlot(pbmc, features = c("CD3D","IL32"))

#Feature plot

FeaturePlot(pbmc, features = c("CD3D"))

# write clusters labels for each cell in to a file

clustersfilepath = "D:/clusterlabel.txt"
clusters<-data.frame(Idents(pbmc))
write.table(clusters,clustersfilepath,sep="\t",col.names = FALSE)
