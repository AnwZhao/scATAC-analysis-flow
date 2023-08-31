library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(magick)
library(showtext)
library(shadowtext)
showtext_auto()

#source("updir.R")


#111111111111111111111111111
#######read#######
read_1 <- function(path,upp,min.cells=10,min.features=200){
  
  counts <- Read10X_h5(filename = paste0(path,'/filtered_feature_bc_matrix.h5'))
  metadata <- read.csv(
  file =  paste0(path,'/singlecell.csv'),
  header = TRUE,
  row.names = 1
  
)


  chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments =  paste0(path,'/fragments.tsv.gz'),
  min.cells = min.cells,
  min.features = min.features
)
  

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
  setwd(upp)
  dir.create('figure')

#gtf <- rtracklayer::import(paste0(path,'/ITAG4.0_gene_models_add_biotype.gtf'))
gtf_name <- list.files(path=path, pattern='gtf$')
gtf <- rtracklayer::import(paste0(path,'/',gtf_name)) 
  list.files(path='C:/Users/Titanic/Desktop/GUI_data/ATAC/data', pattern='gtf$')
gene.coords <- gtf[gtf$type == 'gene']
#seqlevelsStyle(gene.coords) <- 'UCSC'
#gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
Annotation(pbmc) <- gene.coords

saveRDS(pbmc, paste0(path,'/process.rds'))

}


#222222222222222222222222222222
#######TSS#######
TSS_2 <- function(path,upp,TSS.enrichment=1.2){
  
pbmc <- readRDS(paste0(path,'/process.rds'))
DefaultAssay(pbmc) <- "peaks"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc, fast = FALSE)
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > TSS.enrichment, 'High', 'Low')

#setwd(paste0(upp,'/figure'))
#png( "2_TSSPlot.png")          
p1 = TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
ggsave(paste0(upp,'/figure/2_TSSPlot.png'),p1 ,dpi = 300)
#dev.off()

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

saveRDS(pbmc, paste0(path,'/process.rds'))

}


#3333333333333333333333333333333333
#######QCVlnPlot#######
QCVlnPlot_3 <- function(path,upp){
  
pbmc <- readRDS(paste0(path,'/process.rds'))
Idents(pbmc) <- "all"  # group all cells together, rather than by replicate

#setwd(paste0(upp,'/figure'))
#png( "/3_QCVlnPlot.png") 
p2 = VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5)
ggsave(paste0(upp,'/figure/3_QCVlnPlot.png'),p2 ,dpi = 300)
#dev.off()
saveRDS(pbmc, paste0(path,'/process.rds'))

}


#444444444444444444444444444444444
Dimpro_4 <- function(path,upp,k.param=40,resolution=1,dims=30){
  
pbmc <- readRDS(paste0(path,'/process.rds'))
########cluster DimPlot#########
pbmc <-  subset(
  x = pbmc,
  subset = peak_region_fragments > 500 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 10 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1.2
)

######################
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:dims, reduction.name = 'umap.atac')
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', k.param=k.param,dims = 2:dims)
pbmc <- FindClusters(object = pbmc, resolution=resolution,verbose = FALSE, algorithm = 3)

p3 = DimPlot(pbmc, reduction = "lsi")
ggsave(paste0(upp,'/figure/4_DimPlot.png'),p3 ,dpi = 300)

p4 = VizDimLoadings(pbmc, dims = 1:2, reduction = "lsi")
ggsave(paste0(upp,'/figure/4_DimLoad.png'),p4 ,dpi = 300)

saveRDS(pbmc, paste0(path,'/process.rds'))
}


#5555555555555555555555555555555
select_table_5 <- function(path){
pbmc <- readRDS(paste0(path,'/process.rds'))
pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

saveRDS(pbmc, paste0(path,'/process.rds'))
return(c(top10))
}



#66666666666666666666666666666666
show_select_6 <- function(path,upp,gene='SL4.0ch00-3282049-3282954'){

pbmc <- readRDS(paste0(path,'/process.rds'))

p6 = VlnPlot(pbmc, 
        features = c(gene), 
        slot = "counts", log = TRUE)
ggsave(paste0(upp,'/figure/6_VlnPlot.png'),p6 )

p7 = FeaturePlot(pbmc, features = c(gene))
ggsave(paste0(upp,'/figure/6_FeaturePlot.png'),p7 )

}



#7777777777777777777777777777777777
show_cluster_7 <- function(path,upp){
pbmc <- readRDS(paste0(path,'/process.rds'))

p8 = DimPlot(pbmc, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP")
ggsave(paste0(upp,'/figure/7_cluster.png'),p8 )


}