###################################################################
library("iSEE")
#library("SingleCellExperiment") # dont need them as "iSEE" have these all
#library("shiny") # dont need them as "iSEE" have these all
###########################################
# Fetch the data from FigShare
dat <- ("https://figshare.com/ndownloader/files/39246662/sce_dlpfc_sgacc_final.RDS")
download.file(dat, destfile = "sce_dlpfc_sgacc_final.RDS")

sce_small <- load("sce_dlpfc_sgacc_final.RDS")
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
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE, 
                                        CustomRowsText = "SNAP25\nSLC17A6\nSLC17A7\nSLC17A8\nGAD1\nGAD2\nDRD1\nDRD2\nAQP4\nGFAP\nCLDN5\nFLT1\nCD163\nSIGLEC1\nC3\nCD74\nCOL1A2\nPDGFRB\nMBP\nPDGFRA\nVCAN\nSKAP1\nCD247", 
                                        ClusterRows = FALSE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2", 
                                        DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = c("broad.class", "neuron"), RowData = character(0), CustomBounds = FALSE, LowerBound = NA_real_, 
                                        UpperBound = NA_real_, AssayCenterRows = FALSE, AssayScaleRows = FALSE, 
                                        DivergentColormap = "purple < black < yellow", ShowDimNames = "Rows", 
                                        LegendPosition = "Bottom", LegendDirection = "Horizontal", 
                                        VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10, 
                                        ShowColumnSelection = FALSE, OrderColumnSelection = FALSE, 
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



######################################

iSEE(sce_small, initial = initial)
