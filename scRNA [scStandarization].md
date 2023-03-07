# scStandarization [scRNA]

Programs: Cellranger, R, Seurat
Project: https://www.notion.so/Office-9c560aa95c594ffd9e83ea794d6e8d62
Topic: Pills

## Pipeline:

### 0. Sample Sheets

- SampleSheet.csv

### 1. CellRanger

- `mkfastq` w/ SampleSheets.csv
    
    ```r
    cellranger mkfastq --id=FolderName\
                         --run=/path/to/tiny_bcl \
                         --csv=SampleSheet.csv
    ```
    
- `count`
    
    ```r
    cellranger count --id=FolderName \
                       --transcriptome=/opt/refdata-gex-GRCh38-2020-A \
                       --fastqs=/home/path/to/outs/fastq_path \
                       --sample=SampleName 
    ```
    

### 2. SoupX

Remove Ambiental Contamination

```r
# Remove Ambiental Contamination -----------------------------------------
library(SoupX)
sc = load10X(path_count_outs)
sc = autoEstCont(sc)
soupx_out = adjustCounts(sc)
```

### 3. scDblFinder

Detection and handling of doublets/multiplets

```r
# Remove doublets -----------------------------------------------
library(scDblFinder)
counts <- Read10X(soupx_out)
sce <- SingleCellExperiment(list(counts=counts))
sce <- scDblFinder(sce)
```

### 4. Seurat

- A. Build Seurat Object
    
    ```r
    # Build Seurat Object -----------------------------------------
    SeuObj <- as.seurat(sce)
    ```
    
- B. Quality Control Metrics
    
    ```r
    # Generate Quality Metrics --------------------------------------------------- 
      # Mitochondrial percent
    SeuObj$percent.mt <- PercentageFeatureSet(SeuObj, pattern = "^MT-")
      # Ribosomal percent
    SeuObj$percent.rb <- PercentageFeatureSet(SeuObj, pattern = "^RP[SL]")
    
    # Make VlnPlot() for c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb")
    ```
    
- C. Filtering
    
    
    ```r
    # Filtering ---------------------------------------------------
        # filter out Mins 
    SeuObj<- subset(x = SeuObj, 
             subset= (nCount_RNA >= 500) & 
                     (nFeature_RNA >= 250) &
                     (percent.mt < 15))
    ```
    
    ```r
    # Filtering ---------------------------------------------------
        # filter out Mins and Max
    nFeature_RNA > 1000 & 
    nFeature_RNA < 4000 & 
    nCount_RNA > 1500 & 
    nCount_RNA < 12.5e03 & 
    percent.mt < 8 & 
    percent.rb < 45)
    ```
    
    ```r
    # Filtering ---------------------------------------------------
        # filter out quantile Mins and Max 
    gup<-quantile(list_all[[n]][[total_values]],probs=(.9))
    gdown<-quantile(list_all[[n]][[total_values]],probs=(.1))
    filtered <- (list_all[[n]][[total_values]] > gdown & 
    							list_all[[n]][[total_values]] < gup)
    ```
    
- E. Gene Expression Processing
    
    ```r
    # Gene Expression processing -------------------------------------------------
    		# IF FOR RNA ASSAY ---------------
    				# Normalize the counts
    SeuObj<- NormalizeData(SeuObj)
    
    			  # Discovers the most variable features (genes)
    SeuObj<- FindVariableFeatures(SeuObj, selection.method= "vst", nfeatures = 2000) 
    
    			   # Cell Cycle Scores
    s.genes <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    SeuObj <- CellCycleScoring(SeuObj, s.features= s.genes, g2m.features= g2m.genes)
    
    		     # Scale Data
    SeuObj <- ScaleData(SeuObj, features= rownames(SeuObj),
                   vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
    
    		     # Dimensionality Reduction 
    SeuObj <- RunPCA(SeuObj, features = VariableFeatures(object = SeuObj))
    
    		     # Clustering 
    # Scclusteval Tool
    SeuObj <- FindNeighbors(SeuObj, dims = 1:30)
    SeuObj <- FindClusters(SeuObj, resolution = 0.5)
    SeuObj <- RunUMAP(SeuObj, dims = 1:30, verbose = F)
    
    		# If FOR SCT ASSAY ---------------
    				 # SCTransform
    SeuObj <- SCTransform(SeuObj, method = "glmGamPoi", verbose = F,
    							vars.to.regress = c("nFeature_RNA", "percent.mt"))
    SeuObj <- RunPCA(SeuObj, verbose = F)
    SeuObj <- RunUMAP(SeuObj, dims = 1:30, verbose = F)
    SeuObj <- FindNeighbors(SeuObj, dims = 1:30, verbose = F)
    SeuObj <- FindClusters(SeuObj, verbose = F, resolution = 0.5)
    ```
    

### 5. Annotation

- Azimuth                → (w/ ref sc, cell level)
- SingleR                 → (w/ ref sc or bulk, cell level)
- Label transfer       → (w/ own ref)
- sc-type                 → (w/ marker genes, cluster level)

### 6. Integration with other samples

- Seurat
    
    ```r
    # Seurat Integration samples ---------------------------------------------------	
    			# RNA assay integration 
    DefaultAssay(SeuObj_1) <- "RNA"
    DefaultAssay(SeuObj_2) <- "RNA"
    SeuObj_list <- list()
    SeuObj_list[[samples_names[1]]] <- SeuObj_1
    SeuObj_list[[samples_names[2]]] <- SeuObj_2
    
    SeuObj_anchors <- FindIntegrationAnchors(object.list = SeuObj_list, dims = 1:30)
    integrated_seurat <- IntegrateData(anchorset = SeuObj_anchors, dims = 1:30)
    
    DefaultAssay(integrated_seurat) <- "integrated"
    integrated_seurat <- ScaleData(integrated_seurat, verbose = F)
    integrated_seurat <- RunPCA(integrated_seurat, npcs = 30, verbose = F)
    integrated_seurat <- RunUMAP(integrated_seurat, reduction= "pca", dims= 1:30)
    integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:30, k.param = 20)
    integrated_seurat <- FindClusters(integrated_seurat, resolution= 0.5))
    
    			# SCT assay integration 
    DefaultAssay(SeuObj_1) <- "SCT"
    DefaultAssay(SeuObj_2) <- "SCT"
    SCT_merge <- merge(x = SeuObj_1, y = SeuObj_2, add.cell.id= samples_names)
    
    SCT_list <- SplitObject(SCT_merge, split.by = "orig.ident")
    SCT_list <- lapply(X = SCT_list, FUN = SCTransform)
    
    features <- SelectIntegrationFeatures(object.list = SCT_list, nfeatures = 3000)
    SCT_list <- PrepSCTIntegration(object.list= SCT_list, anchor.features= features)
    sct_anchors <- FindIntegrationAnchors(object.list = SCT_list, normalization.method = "SCT", anchor.features = features)
    
    integ_sct <- IntegrateData(anchorset= sct_anchors, normalization.method= "SCT")
    integ_sct <- RunPCA(integ_sct, verbose = FALSE)
    integ_sct <- RunUMAP(integ_sct, reduction = "pca", dims = 1:30)
    integ_sct <- FindNeighbors(integ_sct, dims = 1:30, verbose = F)
    integ_sct <- FindClusters(integ_sct, resolution = 0.5)
    
    ```
    
- Harmony
    
    ```r
    # Load Sample Setting ----------------------------------------------------------
        # Set sample names and paths
    samples_names <- c(sample1,sample2)
    seu_objs_path <- c(seu_objs_path1,seu_objs_path2)
    
        # Load Clustered Seurat Objects
    SeuObj_1 <- readRDS(paste0(seu_objs_path[1],"Clustered_SeuObj.rds"))
    SeuObj_2 <- readRDS(paste0(seu_objs_path[2],"Clustered_SeuObj.rds"))
    DefaultAssay(SeuObj_1) <- "RNA"
    DefaultAssay(SeuObj_2) <- "RNA"
    
        # Merge 
    harmony <- merge(SeuObj_1,SeuObj_2)
    ```
    
    ```r
    # Harmony Integration samples --------------------------------------------------	
    			# Normalize, Scale, Cluster 
    harmony <- NormalizeData(harmony, verbose= F)
    harmony <- FindVariableFeatures(harmony, selection.method="vst", nfeatures=2000)
    harmony <- ScaleData(harmony, verbose = F)
    harmony <- RunPCA(harmony, npcs = 30, verbose = F)
    harmony <- RunUMAP(harmony, reduction = "pca", dims = 1:30, verbose = F)
    
    				# Integration
    harmony <- harmony %>% RunHarmony("orig.ident", plot_convergence= F)
    harmony <- harmony %>% 
      RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
      FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
      FindClusters(resolution = 0.5) %>% 
      identity()
    ```
    
- Rliger
    
    ```r
    # Load Sample Setting ----------------------------------------------------------
        # Set sample names and paths
    samples_names <- c(sample1,sample2)
    seu_objs_path <- c(seu_objs_path1,seu_objs_path2)
    
        # Load Clustered Seurat Objects
    SeuObj_1 <- readRDS(paste0(seu_objs_path[1],"Clustered_SeuObj.rds"))
    SeuObj_2 <- readRDS(paste0(seu_objs_path[2],"Clustered_SeuObj.rds"))
    DefaultAssay(SeuObj_1) <- "RNA"
    DefaultAssay(SeuObj_2) <- "RNA"
    
        # Merge 
    liger <- merge(SeuObj_1,SeuObj_2)
    ```
    
    ```r
    # Rliger Integration -----------------------------------------------------------
    liger <- NormalizeData(liger)
    liger <- FindVariableFeatures(liger)
    liger <- ScaleData(liger, split.by = "orig.ident", do.center = F)
    
    liger <- RunOptimizeALS(liger, k = 30, lambda = 5, split.by = "orig.ident") 
    liger <- RunQuantileNorm(liger, split.by = "orig.ident")
    liger <- FindNeighbors(liger,reduction = "iNMF",k.param = 10,dims = 1:30)
    liger <- FindClusters(liger, resolution = 0.5)
    liger <- RunUMAP(liger, dims = 1:ncol(liger[["iNMF"]]), reduction = "iNMF")
    ```