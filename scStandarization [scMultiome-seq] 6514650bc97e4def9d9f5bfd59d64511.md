# scStandarization [scMultiome-seq]

Programs: Cellranger, R, Seurat
Project: https://www.notion.so/Office-9c560aa95c594ffd9e83ea794d6e8d62
Topic: Pills

<aside>
💡 [WORD TEMPLATES HERE](https://drive.google.com/drive/folders/18ojbHb-yUYKsjkiwaYkuE1qfxOi0q9As?usp=share_link)

</aside>

## Pipeline:

### 0. Sample Sheets

RNA_SampleSheet.csv

![Untitled](scStandarization%20%5BscMultiome-seq%5D%206514650bc97e4def9d9f5bfd59d64511/Untitled.png)

ATAC_SampleSheet.csv

![Untitled](scStandarization%20%5BscMultiome-seq%5D%206514650bc97e4def9d9f5bfd59d64511/Untitled%201.png)

Count_SampleSheet.csv

![Untitled](scStandarization%20%5BscMultiome-seq%5D%206514650bc97e4def9d9f5bfd59d64511/Untitled%202.png)

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
sc = load10X(path_count_outs)
sc = autoEstCont(sc)
soupx_out = adjustCounts(sc)
```

### 3. scDblFinder

Detection and handling of doublets/multiplets

```r
# Remove doublets -----------------------------------------------
library(scDblFinder)
sce <- scDblFinder(sce)
```

### 4. Seurat

- A. Build Seurat Object
    
    ```r
    # Build Seurat Object ----------------------------------------------------------
        # Extract Cell Ranger Data 
    seurat_data <- Read10X_h5("path/outs/filtered_feature_bc_matrix.h5") 
    atac_counts <- seurat_data$Peaks                   
    #seurat_data$`Gene Expression`  # We switch gene expression with soupx matrix.
    
        # Create Seurat object with RNA no-contaminated
    SeuObj <- CreateSeuratObject(counts= soupx_out, min.cell= 10)                
    
        # Add in the ATAC-seq data
    grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use), ]
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevelsStyle(annotations) <- 'UCSC'                                                                  
    genome(annotations) <- "hg38"
    
    frag.file <- "path/outs/atac_fragments.tsv.gz"         
    chrom_assay <- CreateChromatinAssay(counts = atac_counts, sep = c(":", "-"),
       genome= 'hg38', fragments= frag.file, min.cells= 10, annotation= annotations)
    
    # Remove cells from Chrom_assay that are not in SeuO due to min.features=200
       
    SeuObj[["ATAC"]] <- chrom_assay
    ```
    
- B. Quality Control Metrics
    
    
    ```r
    # Quality Control -----------------------------------------------
    DefaultAssay(SeuObj) <- "RNA"
    SeuObj[["percent.mt"]] <- PercentageFeatureSet(SeuObj, pattern = "^MT-")
    SeuObj[["percent.rb"]] <- PercentageFeatureSet(SeuObj, pattern = "^RP[SL]")
    
    DefaultAssay(SeuObj) <- "ATAC"
    SeuObj <- NucleosomeSignal(SeuObj)
    SeuObj <- TSSEnrichment(SeuObj)
    
    # Make VlnPlot() plots for c(nCount_RNA, nFeature_RNA, nCount_ATAC, nFeature_ATAC, percent.mt, percent.rb, TSS.enrichment, nucleosome_signal)
    ```
    
    ```r
    # Itziar Version -------------------------------------------------------------
    
    #  Summary Multiome Metrics Definitions 
    summary_metrics<-list()
    for (i in 1:length(seurat_list)){
      pathmeta<-paste0("/SCMultiom_analysis/",seurat_list[i],"/DATA/summary.csv") #From cellranger 
      summary_metrics[[i]]<-read.csv(file=pathmeta, sep=",")
    }
    summary_metrics<- do.call(rbind, summary_metrics)
    
    setwd("/SCMultiom_analysis/AL_samples/")
    ## Estimated Number of Cells
    #Let us start by plotting the estimated number of cells per library. Note this this number will differ a lot from the final number of cells after applying future QC filters.
    pdf("Estimated Number of Cells.pdf",width = 14,height = 10)
    barplot_data(reshape2::melt(summary_metrics[,c("Sample.ID","Estimated.number.of.cells")])) + coord_flip()
    print(paste("Estimated number of cells for all the samples:", sum(summary_metrics$Estimated.number.of.cells)))
    dev.off()
    ## Joint Metrics
    #Here, we can visualize the number of the features (peaks and genes) and feature linkages detected.
    #Note that the feature linkage is calculated taking in account genes with its proximal peaks and between pairs of proximal peaks across the genome (based on correlation). 
    pdf("Linked_genes_features.pdf",width = 14,height = 10)
    barplot_data(reshape2::melt(summary_metrics[,c("Sample.ID",
                                    "Linked.genes",
                                    "Linked.peaks")]))
    dev.off()
    
    linkage_melt <- reshape2::melt(summary_metrics[,c("Sample.ID",
                                       "Feature.linkages.detected")])
    pdf("Feature_linkages_detected.pdf",width = 14,height = 10)
    barplot_data(multiome_melted = linkage_melt)
    dev.off()
    
    # Chromatin Accessibility Metrics Definitions
    
    ## Sequencing 
    #Our first objective is to evaluated the quality and the quantity of our sequenced libraries prior to mapping.
    
    ### Sequenced read pairs
    #Here, you can see the total number of sequenced read pairs assigned to the Multiome ATAC library. The optimal number of read pairs are 25,000 per cell taking in account the technical note.
    
    summary_metrics$ATAC.Mean.raw.read.pairs.per.cell <- round(summary_metrics$ATAC.Mean.raw.read.pairs.per.cell,2)
    
    pdf("ATAC.Mean.raw.read.pairs.per.cell.pdf",width = 14,height = 10)
    ggbarplot(.....)
    dev.off()
    
    #Here, you can see the quality (leveraging the "Q30" variables)  of sequenced read pairs assigned to the Multiome ATAC library.
    
    multiome_melted <- reshape2::melt(summary_metrics[,c("Sample.ID", "ATAC.Q30.bases.in.barcode",
                                          "ATAC.Q30.bases.in.read.1",        
                                          "ATAC.Q30.bases.in.read.2")])
    pdf("ATAC_Q30quality.pdf",width = 14,height = 10)
    ggplot(......)
    dev.off()
    
    ### Valid barcodes
    #Fraction of read pairs with barcodes that match the whitelist after error correction. A value higher than 85% represent a high quality library.
    
    pdf("ATAC_validbarcodes.pdf",width = 14,height = 10)
    ggbarplot(.....)
    dev.off()
    
    ### Percent of duplicates
    #The percentatge of the duplicated data is correlated with the library size and the sequencing depth.
    
    pdf("correlation_ATAC_dupli_reads.pdf",width = 14,height = 10)
    correlation_plot(summary_metrics,"ATAC.Percent.duplicates","ATAC.Sequenced.read.pairs")
    dev.off()
    
    pdf("correlation_ATAC_dupli_peaks.pdf",width = 14,height = 10)
    correlation_plot(summary_metrics,"ATAC.Percent.duplicates","ATAC.Number.of.peaks")
    dev.off()
    
    ## Cells
    
    ###  Median number of fragments per cell
    #High-Quality Fragment: A read-pair with mapping quality > 30, that is not chimerically mapped, has a valid 10x barcode, and maps to any nuclear contig (not mitochondria) that contains at least one gene.
    
    qc_atac_read_pairs <- reshape2::melt(summary_metrics[,c("Sample.ID","ATAC.Mean.raw.read.pairs.per.cell",                      "ATAC.Median.high.quality.fragments.per.cell")])
    pdf("qc_atac_read_pairs.pdf",width = 14,height = 10)
    barplot_data(qc_atac_read_pairs)
    dev.off()
    qc_atac_fraction <- reshape2::melt(summary_metrics[,c("Sample.ID","ATAC.Fraction.of.high.quality.fragments.in.cells","ATAC.Fraction.of.transposition.events.in.peaks.in.cells")])
    
    pdf("qc_atac_fraction.pdf",width = 14,height = 10)
    barplot_data(qc_atac_fraction)
    dev.off()
    
    ### Targeting Metrics
    #The targetting metrics can be summarized by these 4 main score:
      
      #### Total number of peaks on primary contigs. 
     # This number presents a high correlation with the sequencing depth specifically with the high quality fragments.
    pdf("qc_atac_npeaks.pdf",width = 14,height = 10)
    ggbarplot(.......)
    dev.off()
    
    pdf("correlation_atac_peaks_hfragments.pdf",width = 14,height = 10)
    correlation_plot(summary_metrics,"ATAC.Number.of.peaks",
                     "ATAC.Median.high.quality.fragments.per.cell")
    dev.off()
    
    #### Fraction of high quality fragments that overlapp in peaks.
    #The fraction of high quality fragments in cell are expected to be higher 40%. The fraction of transposition events that fall within peaks > 25%
      
    multiome_melted <- reshape2::melt(summary_metrics[,c("Sample.ID", "ATAC.Fraction.of.high.quality.fragments.in.cells",
                                          "ATAC.Fraction.of.high.quality.fragments.overlapping.peaks",        "ATAC.Fraction.of.high.quality.fragments.overlapping.TSS")])
    multiome_melted$value = round(multiome_melted$value,2)    
    
    pdf("Hfragments_in_peaks.pdf",width = 14,height = 10)
    ggbarplot(......)
    dev.off()
    
    #### Fraction of genome in peaks
    
    #This fraction is quite low in all samples. We do not have a reference value to be able to compare them. However, the examples we reviewed had a value of around 2%. 
    
    pdf("Fraction of genome in peaks.pdf",width = 14,height = 10)
    ggplot(......)
    dev.off()
    
    #### Transcription Start Site (TSS)
    #It is expected to obtain a large enrichment around TSS (> 5% EV, red dashed line), as these regions are known to show a high degree of chromatin accessibility compared to the flanking regions.
    
    data_tss <- summary_metrics[, c("Sample.ID", "ATAC.TSS.enrichment.score")]
    data_tss_melted <- reshape2::melt(data_tss)
    pdf("ATAC.TSS.enrichments.pdf",width = 14,height = 10)
    barplot_data(data_tss_melted) + geom_hline(yintercept = 5, linetype = 2, color = "red")
    dev.off()
    
    # Gene Expression Metrics Definitions
    
    ## Median genes per cell
    median_genes_gg <- reshape2::melt(summary_metrics[,c("Sample.ID","GEX.Median.genes.per.cell")])
    pdf("RNA_mediangenes_cells.pdf",width = 14,height = 10)
    barplot_data(median_genes_gg)
    dev.off()
    
    ### Sequenced read pairs
    #Here, you can see the total number of sequenced read pairs assigned to the Multiome Gene Expression library. The optimal number of read pairs are ~20,000 per cell taking in account the technical note.
    summary_metrics$GEX.Mean.raw.reads.per.cell <- round(summary_metrics$GEX.Mean.raw.reads.per.cell,2)
    pdf("RNA_meanrawreads_cell.pdf",width = 14,height = 10)
    ggbarplot(........)
    dev.off()
    
    #Here, you can see the quality (leveraging the "Q30" variables)  of sequenced read pairs assigned to the Multiome Gene Expresion library.
    
    multiome_melted <- reshape2::melt(summary_metrics[,c("Sample.ID", "GEX.Q30.bases.in.barcode",
                                          "GEX.Q30.bases.in.UMI",        
                                          "GEX.Q30.bases.in.read.2")])
    pdf("RNA_Q30metrics.pdf",width = 14,height = 10)
    ggplot(........)
    dev.off()
    
    ### Valid barcodes
    #Fraction of read pairs with barcodes that match the whitelist after error correction. A value higher than 80% represent a high quality library,
    pdf("RNA_valid_barcodes.pdf",width = 14,height = 10)
    ggbarplot(........)
    dev.off()
    
    ### Percent of duplicates
    #The percentatge of the duplicated data is correlated with the library size and the sequencing depth.
    pdf("RNA_correla_dupli_readpairs.pdf",width = 14,height = 10)
    correlation_plot(summary_metrics,"GEX.Percent.duplicates","GEX.Sequenced.read.pairs")
    dev.off()
    pdf("RNA_correla_dupli_totalreads.pdf",width = 14,height = 10)
    correlation_plot(summary_metrics,"GEX.Percent.duplicates","GEX.Total.genes.detected")
    dev.off()
    
    ## Mapping
    #Secondly, we will assess the quality of cellranger's mapping by comparing the percentage of reads mapping to the genome, intergenic regions, intronic and exonic regions across libraries. 
    #Reads mapping to intergenic regions suggest contamination of ambient DNA, while reads mapping to intronic regions may come from pre-mRNAs or mature splice isoforms that retain the introns.
    
    #The fraction of sequenced reads that map to a unique gene in the transcriptome is expected to be higher than the 50%.
    summary_metrics$GEX.Fraction.of.transcriptomic.reads.in.cells <- round(summary_metrics$GEX.Fraction.of.transcriptomic.reads.in.cells,2)
    pdf("RNA_Fraction.of.transcriptomic.read.pdf",width = 14,height = 10)
    ggbarplot(.....)
    dev.off()
    
    multiome_melted <- reshape2::melt(summary_metrics[,c("Sample.ID", "GEX.Reads.mapped.to.genome",  "GEX.Reads.mapped.confidently.to.genome")])
                    
    multiome_melted$value <- round(multiome_melted$value,2)                    
    pdf("RNA_reads_mapped_genome.pdf",width = 14,height = 10)
    ggbarplot(......)
    dev.off()
    
    #The fraction of sequenced reads that map to a unique gene in the transcriptome is expected to be higher than 50% (blue dashed line). 
    #However, the percentatge of reads mapped confidently to intergenic regions and reads mapped antisense to gene are expected to be lower than 30% (red dashed line).
    
    mapping_qc_vars <- c(
      "Sample.ID",
      "GEX.Reads.mapped.confidently.to.intergenic.regions",
      "GEX.Reads.mapped.confidently.to.intronic.regions",
      "GEX.Reads.mapped.confidently.to.exonic.regions",
      "GEX.Reads.mapped.confidently.to.transcriptome",
      "GEX.Reads.mapped.antisense.to.gene")
    mapping_qc_melted <- reshape2::melt(summary_metrics[,mapping_qc_vars])
    pdf("RNA_reads_mapped_info.pdf",width = 14,height = 10)
    ggplot(........)
    dev.off()
    
    # Sequencing Saturation
    #Thirdly, we will plot the number of detected genes per library as a function of the total reads sequenced. We know that this function reaches a plateau, whereby more sequenced reads does not result in more detected genes. In those scenarios, we say that we have sequenced until saturation:
    pdf("RNA_sequencing_saturation.pdf",width = 14,height = 10)
    ggplot(.......)
    dev.off()
    
    ### Some other parameters QC
    
    ig_genes<-getBM(attributes = c('external_gene_name','gene_biotype'), 
                    mart = human)
    ig_genes<-ig_genes[grep(c("^IG"),ig_genes$gene_biotype),]
    ig_genes_names<-ig_genes[ig_genes$external_gene_name%in%rownames(seurat.obj),1]
      percent.ig<-Matrix::colSums(GetAssayData(seurat.obj, "counts")[ig_genes_names,])/Matrix::colSums(seurat.obj)
      percent.ig<-percent.ig*100
      seurat.obj$percent.ig<-percent.ig
    
      seurat.obj$FRIP<-seurat.obj$atac_peak_region_fragments/seurat.obj$atac_fragments
     seurat.obj$FRIP_group[ seurat.obj$FRIP>0.7]<-"high_FRIP"
       seurat.obj$FRIP_group[metad seurat.obj$FRIP<=0.7]<-"low_FRIP"
    
    #file_list [[i]] <- Seurat object 
     file_list [[i]]$log10GenesPerUMI <- log10(file_list [[i]]$nFeature_RNA) / log10(file_list [[i]]$nCount_RNA)
      file_list [[i]]$complexity<-ifelse(file_list [[i]]$log10GenesPerUMI<0.8,"low","ok")
    file_list [[i]]$blacklist_fraction <- FractionCountsInRegion( object = file_list [[i]], regions = blacklist_hg19)
      file_list [[i]]$pct_reads_in_peaks <- file_list [[i]]$atac_peak_region_fragments / file_list [[i]]$atac_fragments * 100
      file_list [[i]]$blacklist_ratio <- file_list [[i]]$blacklist_fraction / file_list [[i]]$atac_peak_region_fragments
      file_list [[i]]$high.tss <- ifelse(file_list [[i]]$TSS.enrichment > 2, 'High', 'Low')
      file_list [[i]]$nucleosome_group <- ifelse(file_list [[i]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
    ```
    
- C. Filtering
    
    
    ```r
    # Filtering -----------------------------------------------------
        # filter out low quality cells
    SeuObj <- subset(x = SeuObj,
      subset = nCount_ATAC < 100000 &
               nCount_ATAC > 1000 &
               nCount_RNA < 25000 &
               nCount_RNA > 1000 &
               nFeature_RNA > 250 &
               percent.mt < 15 &
               nucleosome_signal < 2 &
               TSS.enrichment > 1)
    ```
    
    ```r
    # Doublet Finder OTHER !!!!!
    
    # Itziar Version ----------------------------------------------------------------
    #file_list [[i]] <- Seurat object 
    list_all<-list()
    for (i in 1:length(file_list)){
      list_all [[i]]<-as.SingleCellExperiment(file_list[[i]])
      list_all [[i]] <- scDblFinder(list_all [[i]])
    }
    
    for (i in 1:length(file_list)){
      list_all [[i]] <- addPerCellQCMetrics(list_all [[i]],flatten=FALSE)
    }
    list_2<-list()
    for (i in 1:length(list_all)){
      list_2 [[i]] <- as.Seurat(list_all [[i]], counts = "counts")
      list_2 [[i]]@project.name<-as.character(seurat_list[i])
      list_2 [[i]]@active.ident<-list_2 [[i]]$orig.ident
    }
    
    for (i in 1:length(file_list)){
      file_list [[i]]@meta.data<-list_2 [[i]]@meta.data
    }
    rm(list_all)
    rm(list_2)
    ```
    
    ```r
    # Itziar Version ------------------------------------------------------------
        # filter out low quality cells
    # Filtering Version --------------------------------------------------------
    for (i in 1:length(file_list)){
      new_fil [[i]]<-subset(file_list [[i]], 
                            subset= pct_reads_in_peaks>15 & 
                                    blacklist_ratio<0.05 & 
                                    nucleosome_signal<2 &
                                    TSS.enrichment>2 &
                                    FRIP>0.7 & 
                                    percent.mt<20 & 
                                    scDblFinder.class == "singlet" & 
                                    nCount_ATAC>1000 & 
                                    nCount_RNA<40000 & 
                                    nFeature_RNA>250 &
                                    nCount_RNA>500)
      }
    ```
    
- D. Peak Calling
    
    ```r
    # Signac: "The set of peaks identified using Cellranger often merges distinct peaks that are close together. This can create a problem for certain analyses, particularly motif enrichment analysis and peak-to-gene linkage. To identify a more accurate set of peaks, we can call peaks using MACS2 with the CallPeaks() function"
    
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
    
    # Peak Calling by celltypes -> Scenic Link: https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html
    ```
    
- E. Gene Expression Processing
    
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
    
- F.  DNA Accessibility Processing
    
    ```r
    # DNA accessibility processing -----------------------------------------------
    DefaultAssay(SeuObj) <- "peaks"
    SeuObj <- FindTopFeatures(SeuObj, min.cutoff = q0) 
        # More Peaks More Information -> "q0"
    		# min.cutoff='q5'-> top 95% most common features as the VariableFeatures.
    		# min.cutoff=5 will include features in >5 cells of VariableFeatures.
    SeuObj <- RunTFIDF(SeuObj)
    SeuObj <- RunSVD(SeuObj)
    ```
    

### 5. mono-omic integration

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

### 6. Annotation

- Azimuth
- SingleR
- sc-type

### 7. Integration

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

### 8. Linking peaks to genes

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

[Joint RNA and ATAC analysis: 10x multiomic](https://stuartlab.org/signac/articles/pbmc_multiomic.html)

[What is Cell Ranger ARC? -Software -Single Cell Multiome ATAC + Gene Exp. -Official 10x Genomics Support](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc)

[https://github.com/constantAmateur/SoupX](https://github.com/constantAmateur/SoupX)

https://github.com/plger/scDblFinder

[https://github.com/satijalab/seurat/issues/5916](https://github.com/satijalab/seurat/issues/5916)