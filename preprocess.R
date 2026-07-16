#!/usr/bin/env Rscript

# Xenium preprocessing: 11 JIA samples and 1 STAT5b LOF patient.
# Input : one Xenium output directory per sample (JIA1 ... JIA11, STAT5B)
# Output: merged, Harmony-integrated Seurat object
#
# Cell type annotation and TLS delineation are performed manually on this object,
# as described in Methods. The annotated object is the input of
# Figure3_STAT5b_unified_pipeline.R.

source("utils.R")

RAW_ROOT <- "xenium_outputs"
OUT_RDS  <- "Xenium_JIA11_STAT5B_merged.rds"

file_list  <- sort(list.files(RAW_ROOT, pattern = "^(JIA[0-9]+|STAT5B)$", full.names = TRUE))
sample_ids <- str_extract(file_list, "JIA[0-9]+|STAT5B")

# Per-sample ReadXenium -> QC -> SCTransform
obj_list <- list()

for (i in seq_along(file_list)) {
  data <- ReadXenium(file_list[i], outs = c("matrix","microns"), type = c("centroids","segmentations"))

  segmentations.data <- list(
    centroids    = CreateCentroids(data$centroids),
    segmentation = CreateSegmentation(data$segmentations)
  )
  coords <- CreateFOV(
    coords    = segmentations.data,
    type      = c("segmentation","centroids"),
    molecules = data$microns,
    assay     = "Spatial"
  )

  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = "Spatial")
  xenium.obj[["BlankCodeword"]]   <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]]    <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  xenium.obj[["fov"]] <- coords

  xenium.obj <- subset(
    xenium.obj,
    subset = nCount_Spatial > 0 &
             nFeature_Spatial > 20 & nFeature_Spatial < 5000 &
             nCount_Spatial < 25000
  )

  xenium.obj <- SCTransform(
    xenium.obj,
    assay           = "Spatial",
    conserve.memory = TRUE,
    vst.flavor      = "v2",
    verbose         = FALSE
  )

  obj_list[[sample_ids[i]]] <- xenium.obj
}

var_genes <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 2000, verbose = FALSE)

# Merge and assign sample-level metadata
ste <- merge(x = obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
names(ste@images) <- names(obj_list)

ste$sample_id <- sub("^([^_]+)_.*$", "\\1", rownames(ste@meta.data))

slide_map <- c(
  JIA1="slide1", JIA3="slide1",
  JIA2="slide2", JIA4="slide2",
  JIA5="slide3", JIA6="slide3",
  JIA7="slide4", JIA8="slide4", JIA9="slide4",
  JIA10="slide5", STAT5B="slide5",
  JIA11="slide6"
)
ste$slide   <- unname(slide_map[ste$sample_id])
ste$disease <- ifelse(ste$sample_id == "STAT5B", "STAT5b LOF", "JIA")

ste[["Spatial"]] <- JoinLayers(ste[["Spatial"]])

ste <- SCTransform(
  ste,
  assay           = "Spatial",
  conserve.memory = TRUE,
  vst.flavor      = "v2",
  verbose         = FALSE
)

# Integration features from the per-sample objects are kept as the variable features of the merged object
VariableFeatures(ste[["SCT"]]) <- var_genes

set.seed(1234)
ste <- RunPCA(ste, assay = "SCT", features = var_genes, verbose = FALSE)

for (nm in names(ste@images)) {
  ste@images[[nm]]@key <- paste0(nm, "_")
}

# Harmony integration over slide and sample
ste <- harmony::RunHarmony(
  object        = ste,
  group.by.vars = c("slide","sample_id"),
  assay.use     = "SCT"
)

ste <- FindNeighbors(ste, reduction = "harmony", dims = 1:30, k.param = 30)
ste <- RunUMAP(ste, reduction = "harmony", dims = 1:30, n.neighbors = 30L, min.dist = 0.30)

# Louvain clustering over a range of resolutions
snn_pcs <- BuildSNNSeurat(ste[["harmony"]]@cell.embeddings[, 1:30], nn.eps = 0)

resolution_list <- c(0.20, 0.40, 0.60, 0.80, 1.00)

ids_mat <- Reduce(
  cbind,
  parallel::mclapply(
    resolution_list,
    function(res_use) {
      Seurat:::RunModularityClustering(
        SNN                = snn_pcs,
        modularity         = 1,
        resolution         = res_use,
        algorithm          = 3,
        n.start            = 10,
        n.iter             = 10,
        random.seed        = 0,
        print.output       = FALSE,
        temp.file.location = NULL,
        edge.file.name     = NULL
      )
    },
    mc.cores = length(resolution_list)
  )
)

ids_df <- data.frame(ids_mat)
colnames(ids_df) <- sprintf("res_%.2f", resolution_list)
rownames(ids_df) <- rownames(ste@meta.data)
ids_df <- ids_df %>%
  dplyr::mutate(across(everything(), ~ factor(.x, levels = sort(unique(.x)))))

ste <- AddMetaData(ste, ids_df)

ste$seurat_clusters <- ste@meta.data[["res_0.20"]]
Idents(ste) <- ste$seurat_clusters

data("cc.genes", package = "Seurat")
ste <- CellCycleScoring(
  ste,
  s.features   = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident    = FALSE
)
Idents(ste) <- ste$seurat_clusters

# Moran's I per gene
ste <- spatial_autocorr(
  ste,
  neighbors.k      = 30,
  connectivity_key = "nn",
  genes            = NULL,
  mode             = "moran",
  transformation   = TRUE,
  n_perms          = 50,
  corr_method      = "BH",
  assay            = "SCT",
  attr             = "data",
  seed             = 1938493,
  copy             = FALSE
)

saveRDS(ste, file = OUT_RDS)
