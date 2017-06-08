
#################################################
## FAST Genomics data normalization using SCnorm
#################################################

# libraries
library(yaml)
library(reshape2)
# library(Matrix)
library(SCnorm)

# variables
appTag <- "SCnorm"

# functions
readFG <- function(countFile) {
  # require(Matrix)
  countData <- read.table(countFile, header = TRUE)
  geneData <- unique(countData[, 2])
  countData$index <- match(countData[, 2], geneData)
  cellIDs <- sort(unique(countData[, 1]))
  # countMatrix <- Matrix(matrix(data = 0, 
  #                              ncol = length(cellIDs),
  #                              nrow = length(geneData),
  #                              dimnames = list(geneData, cellIDs)))
  countMatrix <- matrix(data = 0, 
                        ncol = length(cellIDs),
                        nrow = length(geneData),
                        dimnames = list(geneData, cellIDs))
  for (i in cellIDs){
    tmp <- countData[countData[, 1]==i,]
    countMatrix[tmp$index, as.character(i)] <- tmp[, 3]
  }
  colnames(countMatrix) <- paste("cellID", colnames(countMatrix), sep  = ".")
  return(countMatrix)
}

writeFG <- function(countMatrix, countFile) {
  # require(Matrix)
  require(reshape2)
  # expression_data <- data.frame(gene = rownames(countMatrix),
  #                               as.matrix(countMatrix),
  #                               stringsAsFactors = F,
  #                               row.names = NULL)
  expression_data <- data.frame(gene = rownames(countMatrix),
                                countMatrix,
                                stringsAsFactors = F,
                                row.names = NULL)
  expression_data <- melt(expression_data, id.vars="gene", measure.vars = colnames(countMatrix))
  expression_data <- as.matrix(expression_data[, c(2, 1, 3)])
  colnames(expression_data) <- c("cellId*Ganzzahl", "entrezId*Ganzzahl", "expressionValue*Zahl")
  expression_data <- expression_data[as.numeric(expression_data[,3])!=0,]
  expression_data <- gsub(" ", "", expression_data)
  expression_data[, 1] <- sub("cellID.", "", expression_data[, 1])
  write.table(expression_data, countFile, row.names = F, col.names = T, quote = F, sep = "\t")
}


#################################################
## Data Input
#################################################

# load the parameters required for analysis

setwd("/opt/config/")
setwd("C:/Users/scholzcl/Documents/GitHub/fg_scnorm/config")
SCnormParams <- yaml.load_file("SCnorm_config.yml")


setwd("/fastgenomics/input")
setwd("C:/Users/scholzcl/Downloads/Geissmann_with_published_clusters_and_custom_clustering_Schema_Version_2_1")
# load:
# - expression matrix
countMatrix <- readFG(countFile = SCnormParams$EXPRESSION_FILE_NAME)

# - cell metadata (incl. batch definition)
cell_metadata <- read.table(SCnormParams$CELL_META_FILE_NAME, 
                            header = TRUE,
                            stringsAsFactors = FALSE,
                            sep = "\t")
tmp <- strsplit(colnames(cell_metadata), 
                split = ".",
                fixed = TRUE)
colnames(cell_metadata) <- sapply(tmp, "[[", 1)

# - dataset manifest
manifest <- yaml.load_file(SCnormParams$MANIFEST_FILE_NAME)

cell_conditions <- if("batch_column" %in% names(manifest$data$cell_metadata)) {
  cell_metadata[, manifest$data$cell_metadata$batch_column]
} else {
  NULL
}

#################################################
## Calculations
#################################################

# normalize the data
normCountMatrix <- SCnorm(Data = countMatrix, 
                          Conditions = cell_conditions, 
                          OutputName = NULL,
                          SavePDF = TRUE, 
                          PropToUse = SCnormParams$PROPTOUSE, 
                          Tau = SCnormParams$TAU, 
                          reportSF = SCnormParams$REPORTSF,
                          FilterCellNum = SCnormParams$FILTERCELLNUM, 
                          K = SCnormParams$FILTERCELLNUM, 
                          NCores = NULL, 
                          FilterExpression = SCnormParams$FILTEREXPRESSION,
                          Thresh = SCnormParams$THRESH, 
                          ditherCounts = SCnormParams$DITHERCOUNTS, 
                          withinSample = SCnormParams$WITHINSAMPLE,
                          useSpikes = SCnormParams$USESPIKES)


#################################################
## Data Output
#################################################

setwd("/fastgenomics/output")
# write:
# - normalized expression matrix
