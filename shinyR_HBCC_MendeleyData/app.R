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

# FigShare :
#dat <- ("https://figshare.com/ndownloader/files/39246662/sce_dlpfc_sgacc_final.RDS")
#download.file(dat, destfile = "sce_dlpfc_sgacc_final.RDS")

# MendeleyData:

#dat <- ("https://data.mendeley.com/public-files/datasets/4pmcfgy9ss/files/83d489c7-2773-4cc4-afcf-62e2d9f8f07d/file_downloaded")
#download.file(dat, destfile = "sce_dlpfc_sgacc_final.RDS")
#sce_small <- readRDS("sce_dlpfc_sgacc_final.RDS")

dat <- ("https://data.mendeley.com/api/datasets/4pmcfgy9ss/draft/files/645e92de-3314-4745-9514-b922d0ece2d2/file_downloaded")
download.file(dat, destfile = "sce_dlpfc_sgacc_final_DietSuerat.RDS")
sce_small <- readRDS("sce_dlpfc_sgacc_final_DietSuerat.RDS")

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

