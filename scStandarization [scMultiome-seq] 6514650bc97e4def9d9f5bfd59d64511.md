# scStandarization [scMultiome-seq]

Programs: Cellranger, R, Seurat
Project: https://www.notion.so/Office-9c560aa95c594ffd9e83ea794d6e8d62
Topic: Pills

## Pipeline Guide:

### 0. Sample Sheets

- RNA_SampleSheet.csv
- ATAC_SampleSheet.csv
- Count_SampleSheet.csv

### 1. CellRanger

- `mkfastq` w/ RNA_SampleSheets.csv
    
    ```r
    cellranger  mkfastq --run=$Path_BCL_Input$RNAtag/  
    									  --output-dir=$F_out/CellRanger/Fastq/RNA/ 
                        --samplesheet=$F_out/config/RNA_SampleSheet_mk.csv
    ```
    
- `mkfastq` w/ ATAC_SampleSheets.csv
    
    ```r
    cellranger  mkfastq --run=$Path_BCL_Input$ATACtag/  
                        --output-dir=$F_out/CellRanger/Fastq/ATAC/ 
                        --samplesheet=$F_out/config/ATAC_SampleSheet_mk.csv
    ```
    
- `count` w/ Count_SampleSheet.csv
    
    ```r
    cellranger  count --id=$SampleName 
                      --reference=$Path_Ref_Genome_CR 
                      --libraries=$F_out/config/SampleSheet_cnt.csv
    ```
    

### 2. SoupX

Remove Ambiental Cont.

```r
# Remove Ambiental Contamination -----------------------------------------
		# If you have some 10X data which has been mapped with cellranger
library(SoupX)
sc = load10X(path_count_outs)
sc = autoEstCont(sc)
soupx_out = adjustCounts(sc)
```

### 3. scDblFinder

Detection and handling of doublets/multiplets

RNA assay

```r
# Remove doublets in RNA -----------------------------------------------
library(scDblFinder)
soupx_out <- as.SingleCellExperiment(soupx_out)
sce_RNA <- scDblFinder()
```

ATAC assay

```r
# Remove doublets in ATAC -----------------------------------------------
library(scDblFinder)
seurat_data <- Read10X_h5("path/outs/filtered_feature_bc_matrix.h5") 
atac_counts <- as.SingleCellExperiment(seurat_data$Peaks)
sce_ATAC <- scDblFinder(atac_counts)
```

### 4. Seurat

- A. Build Seurat Object
    
    
    RNA assay
    
    ```r
    # Build RNA Seurat Object ----------------------------------------------------------
    sce_RNA <- addPerCellQCMetrics(sce_RNA)
        # Create Seurat object with RNA no-contaminated
    SeuObj <- as.seurat(sce_RNA)                
    
    ```
    
    ATAC assay
    
    ```r
    # Build ATAC Seurat Object ----------------------------------------------------------
    sce_ATAC <- addPerCellQCMetrics(sce_ATAC)
    	 # Create Seurat object with ATAC
    grange.counts <- StringToGRanges(rownames(sce_ATAC), sep = c(":", "-"))
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use), ]
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevelsStyle(annotations) <- 'UCSC'                                                                  
    genome(annotations) <- "hg38"
    
    frag.file <- "path/outs/atac_fragments.tsv.gz"         
    chrom_assay <- CreateChromatinAssay(counts = sce_ATAC , sep = c(":", "-"),
       genome= 'hg38', fragments= frag.file, min.cells= 10, annotation= annotations)
    
    SeuObj[["ATAC"]] <- chrom_assay
    ```
    
- B. Quality Control Metrics
    
    
    RNA assay
    
    ```r
    # Quality Control RNA -----------------------------------------------
    DefaultAssay(SeuObj) <- "RNA"
    SeuObj[["percent.mt"]] <- PercentageFeatureSet(SeuObj, pattern = "^MT-")
    SeuObj[["percent.rb"]] <- PercentageFeatureSet(SeuObj, pattern = "^RP[SL]")
    
    # Make VlnPlot() plots for c(nCount_RNA, nFeature_RNA, nCount_ATAC, nFeature_ATAC, percent.mt, percent.rb, TSS.enrichment, nucleosome_signal...)
    ```
    
    ATAC assay
    
    ```r
    # Quality Control ATAC -----------------------------------------------
    DefaultAssay(SeuObj) <- "ATAC"
    SeuObj <- NucleosomeSignal(SeuObj)
    SeuObj <- TSSEnrichment(SeuObj)
    
    # Make VlnPlot() plots for c(nCount_RNA, nFeature_RNA, nCount_ATAC, nFeature_ATAC, percent.mt, percent.rb, TSS.enrichment, nucleosome_signal...)
    ```
    
- C. Filtering
    
    ```r
    # Filtering Version --------------------------------------------------------
    SeuObj <- subset(x = SeuObj,
     subset = pct_reads_in_peaks > 15 & 
    	        blacklist_ratio < 0.05 & 
    	        nucleosome_signal < 2 &
    	        TSS.enrichment > 2 &
    	        FRIP > 0.7 & 
    	        percent.mt < 20 & 
    	        scDblFinder.class == "singlet" & 
    	        nCount_ATAC > 1000 & 
    	        nCount_RNA < 40000 & 
    	        nFeature_RNA > 250 &
    	        nCount_RNA > 500)
    ```
    
- D. Mono-Omic Processing
    
    
    Gene Expression Precessing
    
    ```r
    # Gene Expression processing -------------------------------------------------
    DefaultAssay(SeuObj) <- "RNA"
    
        # Cell Cycle Scores
    s.genes <- cc.genes.updated.2019$s.genes                                                         
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    SeuObj <- NormalizeData(SeuObj)
    SeuObj <- CellCycleScoring(SeuObj, s.features= s.genes, g2m.features= g2m.genes)
    
        # SCTransform
    SeuObj <- SCTransform(SeuObj, method = "glmGamPoi", 
    											vars.to.regress = c("nFeature_RNA", "percent.mt"))
    SeuObj <- RunPCA(SeuObj)
    SeuObj <- RunUMAP(SeuObj, dims = 1:30)
    ```
    
    ATAC assay DNA Accessibility Processing
    
    ```r
    # Peak Calling ---------------------------------------------------------------
        # call peaks using MACS2
    peaks <- CallPeaks(SeuObj, macs2.path= Path_Macs2, outdir= some/path)
    
        # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    peaks <- subsetByOverlaps(x= peaks, ranges= blacklist_hg38_unified, invert= T)
    
        # quantify counts in each peak
    macs2_counts <- FeatureMatrix(fragments = Fragments(SeuObj), features = peaks,
                                  cells = colnames(SeuObj))
    
        # create new assay using the MACS2 peak set and add it to the Seurat object
    SeuObj[["peaks"]] <- CreateChromatinAssay(counts = macs2_counts, fragments = frag.file, annotation = annotations)
    
    #SCENIC!
    # DNA accessibility processing -----------------------------------------------
    DefaultAssay(SeuObj) <- "peaks"
    SeuObj <- FindTopFeatures(SeuObj, min.cutoff = q0) 
    SeuObj <- RunTFIDF(SeuObj)
    SeuObj <- RunSVD(SeuObj)
    ```
    
- E. Mono-Omic Integration
    
    
    RNA 
    
    ```r
    # RNA -------------------------------------------------------------------------
    list_Obj_RNA <- list(seuObj1,seuObj1,seuObj1)
    list_Obj_RNA <- SelectIntegrationFeatures(list_Obj_RNA)
    
    # Instration for 
    SeuObj_anchors <- FindIntegrationAnchors(object.list = SeuObj_list, dims = 1:30)
    integrated_seurat <- IntegrateData(anchorset = SeuObj_anchors, dims = 1:30)
    
    # For quality Check
    DefaultAssay(integrated_seurat) <- "integrated"
    integrated_seurat <- ScaleData(integrated_seurat, verbose = F)
    integrated_seurat <- RunPCA(integrated_seurat, npcs = 30, verbose = F)
    integrated_seurat <- RunUMAP(integrated_seurat, reduction= "pca", dims= 1:30)
    integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:30, k.param = 20)
    integrated_seurat <- FindClusters(integrated_seurat, resolution= 0.5))
    
    ```
    
    ATAC
    
    ```r
    
    # ATAC -------------------------------------------------------------------------
    combined.peaks <- reduce(x = c(peaks_sample1, peaks_sample2, peaks_sampleN)) # Calculate peak consensus
    
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevelsStyle(annotations) <- 'UCSC'                                                                  
    genome(annotations) <- "hg38"
    
    # make for statement for each sample
    	atac_basal <- FeatureMatrix(fragments = frags.basal, features = combined.peaks, cells = rownames(metadata_features_basal)) # Recuantify matrisx for peaks consensus
    
        # Add in the ATAC-seq data
    grange.counts <- StringToGRanges(rownames(atac_basal), sep = c(":", "-"))
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use), ]
    
    frag.file <- "path/outs/atac_fragments.tsv.gz"         
    chrom_assay <- CreateChromatinAssay(counts = atac_counts, sep = c(":", "-"),
       genome= 'hg38', fragments= frag.file, min.cells= 10, annotation= annotations)
    
    ```
    

### 5. Annotation

- Azimuth
- SingleR
- sc-type

### 6. Integration

```r
# Build a joint neighbor graph using both assays -------------------------------
SeuObj <- SelectIntegrationFeatures(object = SeuObj)

SeuObj <- FindMultiModalNeighbors(object = SeuObj, 
					reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30), 
					modality.weight.name = "RNA.weight", verbose = T)

  # build a joint UMAP visualization
SeuObj <- RunUMAP(object = SeuObj, nn.name = "weighted.nn",
								  assay = "RNA", verbose = TRUE)

# Clustering
SeuObj <- FindNeighbors(SeuObj, dims = 1:30, verbose = FALSE)
SeuObj_AllRes <- SeuObj
SeuObj <- FindClusters(SeuObj, graph.name= "wsnn", resolution= 0.5, verbose= FALSE)
```

### 7. Linking peaks to genes

```r
# Peak 2 Gene ------------------------------------------------------------------
    # first compute the GC content for each peak
DefaultAssay(SeuObj) <- "peaks"
SeuObj <- RegionStats(SeuObj, genome = BSgenome.Hsapiens.UCSC.hg38)

    # link peaks to genes
SeuObj <- LinkPeaks(object = SeuObj,peak.assay = "peaks",
  expression.assay = "SCT")
```

---

### B**ibliography**:

[6. Quality Control — Single-cell best practices](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html)

[Common Considerations for Quality Control Filters for Single Cell RNA-seq Data - 10x Genomics](https://www.10xgenomics.com/resources/analysis-guides/common-considerations-for-quality-control-filters-for-single-cell-rna-seq-data)

[Joint RNA and ATAC analysis: 10x multiomic](https://stuartlab.org/signac/articles/pbmc_multiomic.html)

[What is Cell Ranger ARC? -Software -Single Cell Multiome ATAC + Gene Exp. -Official 10x Genomics Support](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc)

[https://github.com/constantAmateur/SoupX](https://github.com/constantAmateur/SoupX)

https://github.com/plger/scDblFinder

[https://github.com/satijalab/seurat/issues/5916](https://github.com/satijalab/seurat/issues/5916)