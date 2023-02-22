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
#dat <- ("https://figshare.com/ndownloader/files/39246662/sce_dlpfc_sgacc_final.RDS")
# regular latest Seurat object
dat <- ("https://figshare.com/ndownloader/files/39305303/sce_dlpfc_sgacc_final.RDS")
download.file(dat, destfile = "sce_dlpfc_sgacc_final.RDS")
sce_small <- readRDS("sce_dlpfc_sgacc_final.RDS")
tour <- read.delim("tour.txt", sep=";", stringsAsFactors = FALSE, row.names = NULL)
# FigShare:
#dat <- ("https://figshare.com/ndownloader/files/39307748/sce_dlpfc_sgacc_final_DietSuerat_v1.RDS")
#download.file(dat, destfile = "sce_dlpfc_sgacc_final_DietSuerat_v1.RDS")
#sce_small <- readRDS("sce_dlpfc_sgacc_final_DietSuerat_v1.RDS")

# ################################################
# Specify number of colurs for each cell type
library(RColorBrewer)
n <- 47
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- sample(col_vector, n)
names(col_vector) <- as.vector(unique(sce_small$celltype))
# ################################################
initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "TSNE", XAxis = 1L, 
                                          YAxis = 2L, FacetRowByColData = "Barcode", FacetColumnByColData = "Barcode", 
                                          ColorByColumnData = "ID", ColorByFeatureNameAssay = "logcounts", 
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "sample_ID", 
                                          SizeByColumnData = "sum", FacetRowBy = "None", FacetColumnBy = "None", 
                                          ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                          ColorByFeatureName = "SNAP25", ColorByFeatureSource = "---", 
                                          ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "donor4_AAACCCAAGAGTCTTC.1", 
                                          ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                          ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                          ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE, 
                                          VisualChoices = c("Color", "Shape"), ContourAdd = FALSE, 
                                          ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1, 
                                          Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE, 
                                          CustomLabelsText = "donor4_AAACCCAAGAGTCTTC.1", FontSize = 1, 
                                          LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE, 
                                          LabelCenters = FALSE, LabelCentersBy = "Barcode", LabelCentersColor = "#000000", 
                                          VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version", 
                                                                                                             "numeric_version"))), PanelId = c(ReducedDimensionPlot = 1L), 
                                          PanelHeight = 600L, PanelWidth = 6L, SelectionBoxOpen = FALSE, 
                                          RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                          DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                          RowSelectionRestrict = FALSE, ColumnSelectionRestrict = TRUE, 
                                          SelectionHistory = list())



################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data",
                                      XAxisColumnData = "broad.class", XAxisFeatureName = "SNAP25",
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
                                      YAxisFeatureName = "SNAP25", YAxisFeatureSource = "RowDataTable1",
                                      YAxisFeatureDynamicSource = TRUE, FacetRowByColData = "Barcode",
                                      FacetColumnByColData = "Barcode", ColorByColumnData = "broad.class",
                                      ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
                                      ShapeByColumnData = "sample_ID", SizeByColumnData = "sum", FacetRowBy = "None",
                                      FacetColumnBy = "None", ColorBy = "Column data", ColorByDefaultColor = "#000000",
                                      ColorByFeatureName = "SNAP25", ColorByFeatureSource = "---",
                                      ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "{{cellone}}",
                                      ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
                                      ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
                                      ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
                                      VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
                                      PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
                                      CustomLabels = FALSE, CustomLabelsText = "{{cellone}}",
                                      FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom",
                                      HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "Barcode",
                                      LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
                                        c(2L, 4L, 0L)), class = c("package_version", "numeric_version"
                                        ))), PanelId = c(FeatureAssayPlot = 1L), PanelHeight = 600L,
                                      PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
                                      ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
                                      ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
                                      ColumnSelectionRestrict = TRUE, SelectionHistory = list())
################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE, 
    CustomRowsText = "SATB2\nGAD2\nAQP4\nMOG\nMEGF11\nPTPRC\nFLT1", 
    ClusterRows = TRUE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2", 
    DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = c("neuron", 
    "celltype"), RowData = character(0), CustomBounds = FALSE, 
    LowerBound = NA_real_, UpperBound = NA_real_, AssayCenterRows = FALSE, 
    AssayScaleRows = FALSE, DivergentColormap = "purple < black < yellow", 
    ShowDimNames = "Rows", LegendPosition = "Right", LegendDirection = "Horizontal", 
    VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10, 
    ShowColumnSelection = FALSE, OrderColumnSelection = TRUE, 
    VersionInfo = list(iSEE = structure(list(c(2L, 6L, 0L)), class = c("package_version", 
    "numeric_version"))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 12L, 
    SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
    RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
    SelectionHistory = list())

################################################################################
# Settings for Column data plot 1
################################################################################

initial[["ColumnDataPlot1"]] <- new("ColumnDataPlot", XAxis = "Column data", YAxis = "nFeature_RNA", 
    XAxisColumnData = "celltype", FacetRowByColData = "sample_ID", 
    FacetColumnByColData = "sample_ID", ColorByColumnData = "celltype", 
    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
    ShapeByColumnData = "sample_ID", SizeByColumnData = "nCount_RNA", 
    FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data", 
    ColorByDefaultColor = "#000000", ColorByFeatureName = "RP11-34P13.3", 
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
    ColorBySampleName = "2543_sgACC_2_AAACCTGAGATAGGAG", ColorBySampleSource = "---", 
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
    VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE, 
    ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1, 
    Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE, 
    CustomLabelsText = "2543_sgACC_2_AAACCTGAGATAGGAG", FontSize = 1, 
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE, 
    LabelCenters = FALSE, LabelCentersBy = "sample_ID", LabelCentersColor = "#000000", 
    VersionInfo = list(iSEE = structure(list(c(2L, 6L, 0L)), class = c("package_version", 
    "numeric_version"))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 6L, 
    SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
    SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "SNAP25", Search = "", SearchColumns = c("",
                                                                                                      "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
                                                                                                      "", "", "", "", "", "", "", "", "", "", ""), HiddenColumns = character(0),
                                  VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
                                                                                                     "numeric_version"))), PanelId = c(RowDataTable = 1L), PanelHeight = 600L,
                                  PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
                                  ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
                                  ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
                                  ColumnSelectionRestrict = FALSE, SelectionHistory = list())


######################################

sce_small <- registerAppOptions(sce_small, color.maxlevels = 47)

iSEE(
  sce_small,
  tour = tour,
  appTitle = "HBCC sgACC-DLPFC snRNA-seq study 2023",
  initial = initial,
  colormap = ExperimentColorMap(colData = list(
    celltype = function(n) {
      col_vector[!grepl("drop", names(col_vector))]
    }
  ))
)

############################################################################
