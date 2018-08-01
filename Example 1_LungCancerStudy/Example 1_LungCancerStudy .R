
##################################################################################################################################
#####################################               Author : Ramyar Molania                #######################################
##################################################################################################################################

### The R scripts below are provided to reproduce all results and figures for example 1: lung cancer study reported in our paper:
### 'A new normalization for the Nanostring nCounter gene expression assay'
### https://www.biorxiv.org/content/early/2018/07/23/374173


### Example 1: lung cancer study
### For full detail of this study, please refer to our paper.


### R version 3.3.3 (2017-03-06)
##################################################################################################################################
####################################################################### Required R packages #######################################

Libraries <- c('matrixStats', 'ruv', 'NanoStringNorm', 'data.table','NormqPCR',
               'ComplexHeatmap','dplyr', 'plotrix', 'tuple', 'scales', 'ggplot2')
lapply(Libraries, require, character.only = TRUE)
rm(Libraries)

##################################################################################################################################
##################################################################  Reading sample annotation and raw counts #####################

#***** Reading Nanostring gene expression raw data
Nano_ExpressionMatrix <- read.delim('LungCancer_Nanostring_RawCounts.txt',
                                    stringsAsFactors = FALSE, header = TRUE, as.is = TRUE)
dim(Nano_ExpressionMatrix) # 614 165
# all Endogenous genes (600) and Nanostring negative and positive spike-in controls (14)


#***** Reading sample and clinical information
Nano_SampleInfo <- read.delim('LungCancer_Nanostring_SampleInformation.txt',
                              stringsAsFactors = FALSE, header = TRUE, as.is = TRUE)
dim(Nano_SampleInfo) # 162  17
table(Nano_SampleInfo$Tissues)
# lung cancer    normal_lung tissue
#  146          16


#***** Expression matrix
Nano_RawCounts <- as.matrix(Nano_ExpressionMatrix[ , 4:ncol(Nano_ExpressionMatrix)])
dim(Nano_RawCounts) # 614 162
sum(Nano_RawCounts == 0) # 0
row.names(Nano_RawCounts) <- Nano_ExpressionMatrix$Name
Nano_SampleInfo$SampleNames <- colnames(Nano_RawCounts)


#***** Frequency of cancer and normal tissues
length(unique(Nano_SampleInfo$Patient.unique.barcodes [ Nano_SampleInfo$Tissues == 'normal-lung'])) # 12
length(unique(Nano_SampleInfo$Patient.unique.barcodes [ Nano_SampleInfo$Tissues == 'Cancer']))  # 96



##################################################################################################################################
##################################################################  Exploratory  data analysis – raw counts #####################
#***** Colors for each cartridges
Color_Batches <- c('purple','orange','darkred','blue','chartreuse',
                   'darkgoldenrod4','tan2','darkgreen','red3','darkmagenta',
                   'deeppink','violet','navy','red','dodgerblue')

#***** box plot of raw data - Endogenous genes only
RawCounts_log <- log2(Nano_RawCounts[1:587 , ]) # Excluding 13 housekeeping genes and Nanostring spike-ins
par(mar = c(6.5,6.5,2.3,0), mgp = c(3.7 , 1 , 0))
boxplot(RawCounts_log, las = 1, cex.axis = 2, ylab = '' , xlab = '', cex.lab = 4,
        xaxt = 'n', yaxt = 'n', main = 'Unnormalized counts', cex.main = 3.5,
        outline = FALSE, names = FALSE, frame = FALSE,
        whisklty = 3, whisklwd = 1.5, staplelty = 1, notch = TRUE, boxlwd = 2,
        staplelwd = 0 , boxcol = Color_Batches[factor(Nano_SampleInfo$Cartridges)],
        border = Color_Batches[factor(Nano_SampleInfo$Cartridges)] , col='gray87')
box(lwd = 7, bty = 'l')
axis(1, cex.axis=1, at = c(1, seq(20,166,20)), cex.axis = 2.5, lwd.ticks = 4, mgp = c(3.5,1.6,0))
axis(2, at = c(0, seq(3,15,3)), mgp = c(3.5,.9,0), lwd.ticks = 4, las = 1, cex.axis=3)
mtext(expression(paste('Samples', '(', 'n'[samples], '=', '162', ')')), 1, line = 4.5, cex = 2.5)
mtext(expression(paste(Log[2],' (raw counts)')), 2, line = 3.5, cex = 2.7)


#***** RLE plot - Figure 1 A - Unormalized
par(mar = c(6.5,6.5,2.3,0))
boxplot(RawCounts_log - rowMedians(RawCounts_log),
        main = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', ylim = c(-4.1,4.1),
        outline = FALSE, names = FALSE, frame = FALSE, whisklty = 3, whisklwd = 1.5, staplelty = 1, notch = TRUE, boxlwd = 2,
        staplelwd = 0 , boxcol = Color_Batches[factor(Nano_SampleInfo$Cartridges)],
        border = Color_Batches[factor(Nano_SampleInfo$Cartridges)], col = 'gray87')
box(lwd = 7, bty = 'l')
title('Unnormalized counts', line = -2, cex.main = 3.5)
Median_RawData <- apply(RawCounts_log - rowMedians(RawCounts_log), 2, median)
points(c(1:ncol(RawCounts_log)), Median_RawData, col = Color_Batches[factor(Nano_SampleInfo$Cartridges)], pch = 19, cex = 1.2)
axis(2, mgp = c(3.5, .9 ,0), lwd.ticks=6, las=1, cex.axis=3)
mtext('RLE', 2, line = 3.5, cex = 3.5)
abline(h = 0, col = 'black', lwd = 5, lty = 2)
par(lwd = 3)
axis.break(2, -4.2, style = 'zigzag', brw = .02)
legend(160, 4.1, legend = c(1,2,3,'.','.','.', 13,14,15),
       col = c(Color_Batches[1:3], rep('white', 3), Color_Batches[13:15]),
       pch = 19, bty = 'n', cex = 1.4)
text(x = 162, y = 4.2 ,labels  = 'Cartridges', cex = 1.5)
rm(Median_RawData)


#***** Average plot
### Average of Nanostring positive spike-ins controls
Mean_NegativeControlProbes <- apply(Nano_RawCounts[601:608, ], 2, mean)
### Average of Nanostring positive spike-ins controls
Mean_PositiveControlProbes <- apply(Nano_RawCounts[609:614, ], 2, mean)
### Average of housekeeping genes
Mean_HousekeepingGenes <- apply(Nano_RawCounts[588:600, ], 2, mean)
### Library size
LibrarySize <- colSums(Nano_RawCounts [ 1:600, ])

#### Average plots -  Supplementary Figure 2
par(mar = c(6,7,0,0))
plot(log2(LibrarySize), ylim = c(0,22), bty = 'l', typ = 'n', ylab = '', xlab = '', xaxt = 'n', yaxt = 'n')
axis(1, cex.axis = 1, at = c(1, seq(20,166,20)), cex.axis = 2, lwd.ticks = 4, mgp = c(3.5,1.4,0))
axis(2, cex.axis = 1, cex.axis = 2, lwd.ticks = 4, mgp = c(3.5,1,0), las = 1)
mtext(expression(paste(Log [2], ' (raw counts)')), 2, line = 3.5, cex = 3)
mtext(expression(paste('Samples', ' (', 'n'[samples], '=', '162', ')')), 1, line = 4.5, cex = 2.5)
X <- c(0,6,18,30,42,53,65,77,89,101,113,125,137,143,155,167)
GradiantColors <- paste0('gray', seq(90,20, by = -5))
for(i in 1:15) rect(X[i],-5, X[i+1], 22, col = GradiantColors[i], lty = 0 )
points(log2(LibrarySize), col = alpha('darkgoldenrod1', .8), pch = 15, cex = 1.8, lwd = 1.5)
lines(smooth.spline(c(1:162), log2(LibrarySize), df = 40), col = 'darkgoldenrod1', lwd = 4)
points(log2(Mean_HousekeepingGenes) + 3, cex = 2, col = alpha('red', .6), pch = 19)
lines(smooth.spline(c(1:162), log2(Mean_HousekeepingGenes) + 3, df = 40), col = 'red', lwd = 4)
points(log2(Mean_PositiveControlProbes) - 1, col = alpha('cyan2', .8), pch = 18, cex = 2.5)
lines(smooth.spline(c(1:162), log2(Mean_PositiveControlProbes) - 1, df = 40), col = 'cyan1', lwd = 4)
points(log2(Mean_NegativeControlProbes), col = alpha('green2', .6), pch = 17, cex = 1.8)
lines(smooth.spline(c(1:162), log2(Mean_NegativeControlProbes), df = 40), col = 'green2', lwd = 4)
box(lwd = 5, bty = 'l')
rm(X, GradiantColors, Mean_NegativeControlProbes,
   Mean_PositiveControlProbes, Mean_HousekeepingGenes,
   LibrarySize)


#***** Log ratio between all pairs of duplicated samples
DuplicatedSamples <- tuplicate(Nano_SampleInfo$Patient.barcodes, 2)
length(DuplicatedSamples) # 17

LogRatio_TechRep_RawCounts <- vector()
for(i in 1:17){
  index <- Nano_SampleInfo$Patient.barcodes == DuplicatedSamples[i]
  RepData <- RawCounts_log [ , index]
  LogRatio <- as.data.frame(RepData[ , 1] - RepData[ , 2] )
  colnames(LogRatio) <- 'LogRatio'
  LogRatio$Replicates <- rep( i , 587)
  LogRatio_TechRep_RawCounts <- rbind(LogRatio_TechRep_RawCounts , LogRatio)
  rm(i, index, RepData, LogRatio)
}
nrow(LogRatio_TechRep_RawCounts)/587 # 17
LogRatio_TechRep_RawCounts$Datasets <- 'Unnormalized'


##################################################################################################################################
#################################################  Nanostring normalization - uisng NanoStringNorm R package #####################

#***** Housekeeping genes - 13 genes
Selected_HKgenes <- Nano_ExpressionMatrix$Name[ 588:600]
Nano_ExpressionMatrix$Code.Class[Nano_ExpressionMatrix$Name %in% Selected_HKgenes] <- "Housekeeping"
table(Nano_ExpressionMatrix$Code.Class)
# Endogenous  Housekeeping  Negative     Positive
# 587           13            8            6
rm(Selected_HKgenes)

#***** Nanostring normalization
### These options recomended by Nanostring and are commonly used:
# CodeCount = geo.mean
# Background = mean.2sd
# SampleContent = housekeeping.geo.mean
dim(Nano_ExpressionMatrix) ##  614 165
Nanostring_normalized <- NanoStringNorm(x = Nano_ExpressionMatrix,
                                        CodeCount = 'geo.mean',
                                        Background = "mean.2sd",
                                        SampleContent = 'housekeeping.geo.mean',
                                        round.values = FALSE,
                                        take.log = FALSE ,
                                        return.matrix.of.endogenous.probes = TRUE)

sum(is.na(Nanostring_normalized)) ## 0
NanostringNormalized <- as.matrix(log2(Nanostring_normalized + 1))
all(colnames(NanostringNormalized) == Nano_SampleInfo$SampleNames) # TRUE


#***** RLE plots - Figure 1 A - nCounter normalization
par(mar = c(6.5,6.5,2.3,0))
boxplot(NanostringNormalized - rowMedians(NanostringNormalized), main = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', ylim = c(-4.1,4.1),
        outline = FALSE, names = FALSE, frame = FALSE, whisklty = 3, whisklwd = 1.5, staplelty = 1, notch = TRUE, boxlwd = 2,
        staplelwd = 0 , boxcol = Color_Batches[factor(Nano_SampleInfo$Cartridges)],
        border = Color_Batches[factor(Nano_SampleInfo$Cartridges)], col = 'gray87')
box(lwd = 7, bty = 'l')
title('nCounter normalized', line = -2, cex.main = 3.5)
Median_Nano <- apply(NanostringNormalized - rowMedians(NanostringNormalized), 2, median)
points(c(1:ncol(NanostringNormalized)), Median_Nano, col = Color_Batches[factor(Nano_SampleInfo$Cartridges)], pch = 19, cex = 1.2)
axis(2, mgp = c(3.5, .9 ,0), lwd.ticks=6, las=1, cex.axis=3)
mtext('RLE', 2, line = 3.5, cex = 3.5)
abline(h = 0, col = 'black', lwd = 5, lty = 2)
par(lwd = 3)
axis.break(2, -4.2, style = 'zigzag', brw = .02)
rm(Median_Nano)


#***** Correlation between ERCC1 and RRM1 in cancer samples - Figure 1 C, nCounter normalization
all(colnames(NanostringNormalized) == Nano_SampleInfo$SampleNames)
NanoNor_CancerData <- NanostringNormalized[ , Nano_SampleInfo$Tissues=="Cancer"]
par(mar = c(5,5,2,1))
BioCo_Nano <- cor.test(unlist(t(NanoNor_CancerData['RRM1' , ])) , unlist(t(NanoNor_CancerData[ 'ERCC1' , ])) , method = "spearman")[[4]]
plot (NanoNor_CancerData[ "RRM1", ], NanoNor_CancerData [ "ERCC1" , ],
      las = 1, pch = 21, col = "cyan", bg = 'blue', lwd = 1.2,
      main = "nCounter normalized", bty = "l", cex.main = 2, mgp = c(3.2,.8,0),
      cex.axis = 1.7, cex.lab = 1.5, cex = 2.5, lwd.ticks = 3, ylab = "" ,
      xlab = expression(paste(italic(RRM1) , " (", log[2] , " normalized counts" , ")")))
mtext(expression(paste(italic(ERCC1) , " (", log[2] , " normalized counts" , ")")), 2, line = 2.5, cex = 1.5)
title(paste('(r = ', round(BioCo_Nano, digits = 2), ')'), line = -.5, cex = 1.4)
box(lwd = 6, bty = "l")


#*****  Log ratio between technical replicates sample - duplicated samples
DuplicatedSamples <- tuplicate(Nano_SampleInfo$Patient.barcodes, 2)
length(DuplicatedSamples) # 17

LogRatio_TechRep_NanoNormalized <- vector()
for(i in 1:length(DuplicatedSamples)){
  index <- Nano_SampleInfo$Patient.barcodes == DuplicatedSamples[i]
  RepData <- NanostringNormalized [ , index]
  LogRatio <- as.data.frame(RepData[ , 1] - RepData[ , 2])
  colnames(LogRatio) <- 'LogRatio'
  LogRatio$Replicates <- rep(i , 587)
  LogRatio_TechRep_NanoNormalized <- rbind(LogRatio_TechRep_NanoNormalized , LogRatio)
  rm(i, index, RepData, LogRatio)
}
nrow(LogRatio_TechRep_NanoNormalized)/587
LogRatio_TechRep_NanoNormalized$Datasets <- 'nCounter normalization'




##################################################################################################################################
############################################### RUV-III normalization - using technical replicates and control genes #############

#***** Creating replicate matrix
length(Nano_SampleInfo$Patient.barcodes) # 162
length(unique(Nano_SampleInfo$Patient.barcodes)) # 135

ReplicateMatrix <- ruv::replicate.matrix(Nano_SampleInfo$Patient.barcodes)
dim(ReplicateMatrix) # 162 133

### Making sure that every row has got only one number
par(mfrow = c(2,1))
barplot(colSums(ReplicateMatrix))
barplot(rowSums(ReplicateMatrix))
par(mfrow = c(1,1))


#***** Finding control genes
### Step 1: Using all genes as a set of negative control genes
# Performing RUVIII
dataRUV <- t(log2(Nano_RawCounts[1:587, ]))
RUVcorrected <- ruv::RUVIII(Y = dataRUV, M = ReplicateMatrix, ctl = c(1:587))
RUVcorrected <- t(RUVcorrected)

### Step 2: Selecting the most stable genes from step 1
LowVarGenes <- apply(RUVcorrected, 1, var)
ControlGenes <- which(LowVarGenes < .5)
length(ControlGenes) # 492 genes


### Step3: Performing RUV-III using the ControlGenes set
RUVcorrected <- RUVIII(Y = dataRUV, M = ReplicateMatrix, ctl = ControlGenes, k = 6)
RUVcorrected <- t(RUVcorrected)
dim(RUVcorrected)
# We tried different numbers of k based on our positive controls
all(colnames(RUVcorrected) == Nano_SampleInfo$SampleNames) ## TRUE


#***** RLE plots - Figure 1 A, RUV-III normalization
par(mar = c(6.5,6.5,2.3,0))
boxplot(RUVcorrected - rowMedians(RUVcorrected), main = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', ylim = c(-4.1,4.1),
        outline = FALSE, names = FALSE, frame = FALSE, whisklty = 3, whisklwd = 1.5, staplelty = 1, notch = TRUE, boxlwd = 2,
        staplelwd = 0 , boxcol = Color_Batches[factor(Nano_SampleInfo$Cartridges)],
        border = Color_Batches[factor(Nano_SampleInfo$Cartridges)], col = 'gray87')
box(lwd = 7, bty = 'l')
title('RUV-III normalized', line = -2, cex.main = 3.5)
Median_RUV <- apply(RUVcorrected - rowMedians(RUVcorrected), 2, median)
points(c(1:ncol(RUVcorrected)), Median_RUV, col = Color_Batches[factor(Nano_SampleInfo$Cartridges)], pch = 19, cex = 1.2)
axis(1, cex.axis = 1, at = c(1, seq(20,166,20)), cex.axis = 2.5, lwd.ticks = 6,  mgp = c(3.5,1.6,0))
axis(2, mgp = c(3.5,.9,0), lwd.ticks = 6, las = 1, cex.axis = 3)
mtext(expression(paste('Samples','(', 'n'[samples], '=', '162', ')')),1 ,line = 5, cex = 3)
mtext('RLE', 2, line = 3.5, cex = 3.5)
abline(h = 0, col = 'black', lwd = 5, lty = 2)


#***** Correlation between ERCC1 and RRM1 in cancer samples
RUVIII_CancerData <- RUVcorrected[ , Nano_SampleInfo$Tissues=="Cancer"]
dim(RUVIII_CancerData) ## 587 150
BioCo_RUV <- cor.test(unlist(t(RUVIII_CancerData['RRM1' , ])) , unlist(t(RUVIII_CancerData[ 'ERCC1' , ])), method = "spearman")[[4]]
par(mar = c(5,5,2,1))
plot (RUVIII_CancerData[ "RRM1", ], RUVIII_CancerData [ "ERCC1" , ],
      las = 1, pch = 21, col = "cyan", bg = 'blue', lwd = 1.2,
      main = "RUV-III normalized", bty = "l", cex.main = 2, mgp = c(3.2,.8,0),
      cex.axis = 1.7, cex.lab = 1.5, cex = 2.5, lwd.ticks = 3, ylab = "" ,
      xlab = expression(paste(italic(RRM1) , " (", log[2] , " normalized counts" , ")")))
mtext(expression(paste(italic(ERCC1) , " (", log[2] , " normalized counts" , ")")), 2, line = 3, cex = 1.5)
title(paste('(r = ', round(BioCo_RUV, digits = 2), ')'), line = -.5, cex = 1.4)
box(lwd = 6, bty = "l")


#***** Log ratio between technical duplicated
### Assessing the perfromance of RUV-III using 'leave out one duplicate'
DuplicatedSamples <- tuplicate(Nano_SampleInfo$Patient.barcodes, 2)
length(DuplicatedSamples) # 17
LogRatio_TechRep_RUVIII <- vector()
dataRUV <- t(log2(Nano_RawCounts[1:587, ]))

for(l in 1:length(DuplicatedSamples)){
  index <- Nano_SampleInfo$Patient.barcodes == DuplicatedSamples[l]
  info <- Nano_SampleInfo
  info$Patient.barcodes[index] <- c('Rep1', 'Rep2')
  X <- length(info$Patient.barcodes)
  Y <- length(unique(info$Patient.barcodes))

  ReplicateMatrix <- matrix(0, nrow = X, ncol = Y)
  row.names(ReplicateMatrix) <- info$Patient.barcodes
  colnames(ReplicateMatrix) <- unique(info$Patient.barcodes)
  for(i in 1:Y){
    n <- colnames(ReplicateMatrix)[i]
    sa <- which(row.names(ReplicateMatrix)==n)
    for(j in 1:length(sa)){
      ReplicateMatrix [sa[j] , i] <- 1
    }
  }
  RUVcorrected_Data <- ruv::RUVIII(dataRUV , ReplicateMatrix , ctl = ControlGenes, k = 6)
  RUVcorrected_Data <- t(RUVcorrected_Data)
  index <- which(index)
  LogRatio <- as.data.frame( (RUVcorrected_Data[ , index[1]] ) - (RUVcorrected_Data[ , index[2]]) )
  colnames(LogRatio) <- 'LogRatio'
  LogRatio$Replicates <- rep(l , 587)
  LogRatio_TechRep_RUVIII <- rbind(LogRatio_TechRep_RUVIII , LogRatio)
  rm(i, index, info, X, Y, LogRatio)

}
nrow(LogRatio_TechRep_RUVIII)/587
LogRatio_TechRep_RUVIII$Datasets <- 'RUV-III normalized'


#***** Technical Replicate Agreement plot for unnormalized, Nanostring Normalization and RUV-III (leave out one technical replicate)
### Figure 1 B
TRAallNormalizations <- rbind(LogRatio_TechRep_RawCounts , LogRatio_TechRep_NanoNormalized, LogRatio_TechRep_RUVIII)
TRAallNormalizations$Replicates <- as.factor(TRAallNormalizations$Replicates)
TRAallNormalizations$Datasets <- as.factor(TRAallNormalizations$Datasets)
TRAallNormalizations$Datasets  <- factor(TRAallNormalizations$Datasets , levels = c('Unnormalized','nCounter normalization','RUV-III normalized'))
ggplot(data = TRAallNormalizations) +
  geom_boxplot(aes(y = LogRatio, x = Replicates, fill = Datasets, color=Datasets), outlier.size = .5,
               position = position_dodge(width = .7),
               width = 1, alpha = 0.8, lwd = .7) +
  ylab('Log ratio')+ylim(c(-10,10)) +
  scale_color_manual(values = c('lightblue', 'royalblue', 'navy')) +
  scale_fill_manual(values = c('lightblue', 'royalblue', 'navy')) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2, lwd = .6, col = 'red')+
  theme(axis.title = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.line = element_line(size = 1.5),
        legend.title = element_text(size = 20, face = c('italic', 'bold')),
        legend.text = element_text(size = 20),
        legend.position = "top", legend.direction = "horizontal")



##################################################################################################################################
###################################################### All 84 different normalization - using NanoStringNorm package #############

#***** All different Nanostring normalization options using NanoStringNorm.
### Step 1: Creating a matrix of normalization options
CodeCount_Option <- c('none', 'sum', 'geo.mean')
BackGround_Option <- c('none', 'mean', 'mean.2sd', 'max')
HK_Option <- c('none', 'housekeeping.sum', 'housekeeping.geo.mean', 'total.sum', 'low.cv.geo.mean', 'top.mean', 'top.geo.mean')
AllOptions <- expand.grid(CodeCount_Option, BackGround_Option, HK_Option)
dim(AllOptions) # 84  3


### Step 2: Extracting medians of rle data, log ration of duplicates and correlation between ERCC1 and RRM1
Medians_RLE <- vector()
All_TechReo_LogRatio <- vector()
BiologyCor <- vector()
for(i in 1:84){
  NanoString_mRNA_norm <- NanoStringNorm(x = Nano_ExpressionMatrix,
                                         CodeCount = as.character(AllOptions[i , 1]),
                                         Background = as.character(AllOptions[i , 2]),
                                         SampleContent = as.character(AllOptions[i , 3]),
                                         round.values = FALSE,
                                         take.log = FALSE ,
                                         return.matrix.of.endogenous.probes = TRUE)

  NanoNormalized <- as.matrix(log2(NanoString_mRNA_norm + 1))
  ### RLE transformation
  MediansOfRLE <- apply( NanoNormalized - rowMedians(NanoNormalized), 2, median)
  Medians_RLE <- cbind(Medians_RLE, MediansOfRLE)
  ### Correlation between ERRC1 and RRM1
  CancerSamples <- Nano_SampleInfo$Tissues == "Cancer"
  NanoNormalized_Cancer <- NanoNormalized[ , CancerSamples]
  BiologyCor[i] <- cor.test(unlist(t(NanoNormalized_Cancer['RRM1' , ])) , unlist(t(NanoNormalized_Cancer[ 'ERCC1' , ])), method = "spearman")[[4]]
  ### Log ratio of duplicates samples
  data <- NanoString_mRNA_norm + 1
  LogRa <- vector()
  for(j in 1:length(DuplicatedSamples)){
    index <- Nano_SampleInfo$Patient.barcodes == DuplicatedSamples[j]
    RepData <- data [ , index]
    LogRatio <- log2(RepData[ , 1]/RepData[ , 2])
    LogRa <- c(LogRatio, LogRa)
  }
  All_TechReo_LogRatio <- cbind(All_TechReo_LogRatio , LogRa)
}


### Step 3: Abbreviations for different normalizations
AllOptions$Var1 <- ifelse(AllOptions$Var1=='none', 'N', ifelse(AllOptions$Var1 == 'sum', 'S', 'GM'))
AllOptions$Var2 <- ifelse(AllOptions$Var2=='none', 'N', ifelse(AllOptions$Var2 == 'mean', 'M', ifelse(AllOptions$Var2 == 'mean.2sd', 'M2SD', 'Ma') ))
HK_OP <- c("none", "housekeeping.sum", "housekeeping.geo.mean", "total.sum", "low.cv.geo.mean", "top.mean", "top.geo.mean")
H <- as.character(c("N", "S", "GM", "TS", "LCGM", "TM", "TGM"))

AllOptions$HKoptions <- 'NA'
for(i in 1:length(HK_OP)){
  index <- grep(HK_OP[i], AllOptions$Var3)
  AllOptions$HKoptions [index] <- H[i]
}
AllOptions <- AllOptions[, c(1,2,4)]
FinalName <- paste(AllOptions$Var1, AllOptions$Var2, AllOptions$HKoptions, sep = ',')

### Step 4: Generating all elemnets for fiagure
# RLE
rle <- as.data.frame(Medians_RLE)
colnames(rle) <- FinalName
rle$RUV_III <- apply( RUVcorrected - rowMedians(RUVcorrected) , 2 , median)

# Correlation between ERCC1 amd RRM1
bioCor <- unlist(c(BiologyCor, BioCo_RUV))

# Replicates agreements
Rep <- as.data.frame(All_TechReo_LogRatio)
Rep$RUV_III <- LogRatio_TechRep_RUVIII$LogRatio


### Step 5: Generating the Supplementary 3
Correlation <- HeatmapAnnotation(points = anno_points(bioCor, ylim = c(-.7, .7), pch = 16, size = unit(.5, 'cm'), border = FALSE,
                                                      cex.axis =  4 , axis = TRUE, axis_side = 'left' ,axis_gp = gpar(fontsize = 16, lwd = 6),
                                                      gp = gpar(col = 'darkgreen') ))

Replicates <- HeatmapAnnotation(boxplot = anno_boxplot(Rep, outline = T, ylim = c(-10.5 , 10.5),
                                                       axis_gp = gpar(fontsize = 16, lwd = 6),pch=21, size = unit(.3, "mm"),
                                                       axis = TRUE, axis_side = 'left', border = FALSE,
                                                       gp = gpar(col= 'navy', lwd = 3, lty = 1)))

Heatmap(rle, name = '', cluster_columns = FALSE, cluster_rows = FALSE, show_row_names = FALSE,
        top_annotation = Correlation, bottom_annotation = Replicates,
        top_annotation_height = unit(7, "cm"), bottom_annotation_height = unit(12, "cm"),
        col = c('firebrick1', 'green', 'black', 'green'),
        column_title = expression(paste('Spearman correlation between ', italic(ERCC1), ' and ', italic(RRM1))),
        row_title = 'Medians of RLE plots',
        column_title_gp = gpar(fontsize = 23), row_title_gp = gpar(fontsize = 22),
        row_title_side = 'left', rect_gp = gpar(col= "white", lwd = .1) ,
        heatmap_legend_param = list(title = '', color_bar = "discrete", labels_gp = gpar(fontsize = 24),
                                    grid_width = unit(1, 'cm'), grid_height = unit(1, 'cm')) ,
        width = unit(42, "cm") )

decorate_annotation("points", {
  grid.text("Correlation coefficient",  y = unit(35, "mm"), x = unit(-.05, "npc"), just = "center", rot = 90, gp=gpar(fontsize = 20, col="black"))
  grid.lines(c(0, 1), unit(c(0, 0), "native"), gp = gpar(lty = 2, col = "black"))
})
decorate_annotation("boxplot", {
  grid.text("Log ratio",  y = unit(42, "mm"), x = unit(-.03, "npc"), unit( -1, "mm"), just = 'center', rot = 90, gp=gpar(fontsize = 25, col="black"), hjust = 0)
  grid.lines(c(0, 1), unit(c(0, 0), "native"), gp = gpar(lty = 1, col = "red"))

})




