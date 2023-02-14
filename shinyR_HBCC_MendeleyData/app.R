###################################################################
library("iSEE")
#library("SingleCellExperiment") # dont need them as "iSEE" have these all
#library("shiny") # dont need them as "iSEE" have these all
###########################################
# Fetch the data from MendeleyData
dat <- ("https://data.mendeley.com/api/datasets/4pmcfgy9ss/draft/files/eae29684-04fc-43af-b41d-2a7631a5cac0/SCE_DLPFC-n3_tran-etal.rda")
download.file(dat, destfile = "SCE_DLPFC-n3_tran-etal.rda")

sce_small <- load("SCE_DLPFC-n3_tran-etal.rda")
sce_small <- sce.dlpfc.tran


###########################################
# read colors
cell_colors <- readRDS("cell_colors_dlpfc.rds")
stopifnot(packageVersion("iSEE") >= "2.4.0")

## Related to https://github.com/iSEE/iSEE/issues/568
colData(sce_small) <- cbind(
  colData(sce_small)[, !colnames(colData(sce_small)) %in% c("donor", "cellType")],
  colData(sce_small)[, c("cellType", "donor")]
)

initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "PCA_corrected", XAxis = 1L,
                                          YAxis = 2L, FacetRowByColData = "Barcode", FacetColumnByColData = "Barcode",
                                          ColorByColumnData = "cellType", ColorByFeatureNameAssay = "logcounts",
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "donor",
                                          SizeByColumnData = "sum", FacetRowBy = "None", FacetColumnBy = "None",
                                          ColorBy = "Column data", ColorByDefaultColor = "#000000",
                                          ColorByFeatureName = "SNAP25", ColorByFeatureSource = "---",
                                          ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "{{cellone}}",
                                          ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
                                          ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
                                          ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
                                          VisualChoices = c("Color", "Shape"), ContourAdd = FALSE,
                                          ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
                                          Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
                                          CustomLabelsText = "{{cellone}}", FontSize = 1,
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
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "sum", CustomRows = TRUE,
                                        CustomRowsText = "SNAP25
SLC17A6
SLC17A7
SLC17A8", ClusterRows = FALSE,
                                        ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2",
                                        DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = c("cellType",
                                                                                                           "donor"), RowData = character(0), CustomBounds = FALSE,
                                        LowerBound = NA_real_, UpperBound = NA_real_, AssayCenterRows = FALSE,
                                        AssayScaleRows = FALSE, DivergentColormap = "purple < black < yellow",
                                        ShowDimNames = "Rows", LegendPosition = "Bottom", LegendDirection = "Horizontal",
                                        VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10,
                                        ShowColumnSelection = TRUE, OrderColumnSelection = TRUE,
                                        VersionInfo = list(iSEE = structure(list(c(2L, 4L, 0L)), class = c("package_version",
                                                                                                           "numeric_version"))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 6L,
                                        SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---",
                                        RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
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

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data",
                                      XAxisColumnData = "cellType", XAxisFeatureName = "SNAP25",
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
                                      YAxisFeatureName = "SNAP25", YAxisFeatureSource = "RowDataTable1",
                                      YAxisFeatureDynamicSource = TRUE, FacetRowByColData = "Barcode",
                                      FacetColumnByColData = "Barcode", ColorByColumnData = "cellType",
                                      ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
                                      ShapeByColumnData = "donor", SizeByColumnData = "sum", FacetRowBy = "None",
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




iSEE(
  sce_small,
  appTitle = "M.N. Tran et al 2021, {{region}} region https://bit.ly/LIBD10xHuman",
  initial = initial,
  colormap = ExperimentColorMap(colData = list(
    donor = function(n) {
      cols <- RColorBrewer::brewer.pal(8, "Dark2")
      names(cols) <- paste0("donor", seq_len(8))
      return(cols)
    },
    cellType = function(n) {
      cell_colors[!grepl("drop", names(cell_colors))]
    }
  ))
)

