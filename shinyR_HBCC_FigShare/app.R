###################################################################
library("iSEE")
#library("SingleCellExperiment") # dont need them as "iSEE" have these all
#library("shiny") # dont need them as "iSEE" have these all
###########################################
### Fetch the data from FigShare/ MendeleyData

#To retrieve an option
#getOption('timeout')
#To set an option
options(timeout=600)
# FigShare:
# DietSeura Object
dat <- ("https://figshare.com/ndownloader/files/39246662/sce_dlpfc_sgacc_final.RDS")
# regular latest Seurat object
#dat <- ("https://figshare.com/ndownloader/files/39305303/sce_dlpfc_sgacc_final.RDS")
download.file(dat, destfile = "sce_dlpfc_sgacc_final.RDS")
sce_small <- readRDS("sce_dlpfc_sgacc_final.RDS")

# FigShare:
#dat <- ("https://figshare.com/ndownloader/files/39261101/sce_dlpfc_sgacc_final_DietSuerat.RDS")
#download.file(dat, destfile = "sce_dlpfc_sgacc_final_DietSuerat.RDS")
#sce_small <- readRDS("sce_dlpfc_sgacc_final_DietSuerat.RDS")

# ################################################
# Specify number of colurs for each cell type
library(RColorBrewer)
n <- 47
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- sample(col_vector, n)
names(col_vector) <- as.vector(unique(sce_small$final_celltype))
# ################################################
initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "TSNE", XAxis = 1L, YAxis = 2L, 
                                          FacetRowByColData = "sample_ID", FacetColumnByColData = "sample_ID", 
                                          ColorByColumnData = "broad.class", ColorByFeatureNameAssay = "logcounts", 
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "sample_ID", 
                                          SizeByColumnData = "frac.mito", FacetRowBy = "None", FacetColumnBy = "None", 
                                          ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                          ColorByFeatureName = "SNAP25", ColorByFeatureSource = "---", 
                                          ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "2543_sgACC_2_AAACCTGAGATAGGAG", 
                                          ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                          ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                          ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                          VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                          PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                          CustomLabels = FALSE, CustomLabelsText = "2543_sgACC_2_AAACCTGAGATAGGAG", 
                                          FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
                                          HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "sample_ID", 
                                          LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                            c(2L, 6L, 0L)), class = c("package_version", "numeric_version"
                                            ))), PanelId = c(ReducedDimensionPlot = 1L), PanelHeight = 600L, 
                                          PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---", 
                                          ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, 
                                          ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE, 
                                          ColumnSelectionRestrict = TRUE, SelectionHistory = list())


######################################

sce_small <- registerAppOptions(sce_small, color.maxlevels = 47)

iSEE(
  sce_small,
  appTitle = "HBCC sgACC-DLPFC snRNA-seq study 2023",
  initial = initial,
  colormap = ExperimentColorMap(colData = list(
    celltype = function(n) {
      col_vector[!grepl("drop", names(col_vector))]
    }
  ))
)
############################################################################
