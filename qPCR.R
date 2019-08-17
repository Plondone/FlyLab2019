#TJ Wiegman, SLSTP Research Associate
#Ames Research Center
#Bhattacharya Lab, Summer 2019
#R version 3.6.1

#Data Entry ----
rawData <- list() #store each csv as a data.frame in a list
directoryName <- "/Users/TJ Wiegman/Google Drive/School/SU19/qPCR Results/CSV"
fileVec <- list.files(path = directoryName)
fileVec <- fileVec[fileVec != "desktop.ini"]

for (i in 1:length(fileVec)) {
  #Read each csv in, set columns to proper names and data types
  rawData[[i]] <- read.csv(file = paste0(directoryName, "/", fileVec[i]),
                           stringsAsFactors = FALSE)
  colnames(rawData[[i]]) <- c("Sample", "Gene", "CT") #unify column names
  rawData[[i]]$CT[rawData[[i]]$CT == "Undetermined"] <- NA #kill undefined vals
  rawData[[i]]$CT <- as.numeric(rawData[[i]]$CT) #convert to numbers
  
  #I ran a plate with some samples from Nicole; need to exclude her data
  rawData[[i]]$Sample[grep(pattern = "Nicole",
                           fixed = TRUE,
                           x = rawData[[i]]$Sample)] <- ""
  #Water controls looked good, so let's get rid of them
  rawData[[i]]$Sample[grep(pattern = "Water",
                           fixed = TRUE,
                           x = rawData[[i]]$Sample)] <- ""
}

#Get each gene separately; assume 1:48 are gene1 and 49:96 are gene2
geneData <- list()
for (i in 1:length(rawData)) {
  geneData[[i]] <- rawData[[i]][c(1:48),c(1,3)] #col 1 = sample, col 3 = CT
  geneData[[i + length(rawData)]] <- rawData[[i]][c(49:96),c(1,3)]
  names(geneData)[i] <- rawData[[i]][1,2]
  names(geneData)[i + length(rawData)] <- rawData[[i]][49,2]
}

#Label each sample with unique prep number, and maintain master list
for (j in 1:length(geneData)) {
  sampleFrame <- data.frame("Sample" = NULL, "No" = NULL)
  masterSample <- c()
  for (i in 1:nrow(geneData[[j]])) {
    currentSample <- geneData[[j]]$Sample[i]
    if (currentSample == "") {
      #don't do anything to blank entries
    } else if (currentSample %in% sampleFrame$Sample) { #if seen it before
      n <- which(currentSample == sampleFrame$Sample)
      sampleFrame$No[n] <- sampleFrame$No[n] + 1
      geneData[[j]]$Sample[i] <- paste(currentSample,
                                       sampleFrame$No[n])
    } else { #if we haven't seen it before for this gene
      sampleFrame <- rbind(sampleFrame, data.frame("Sample" = currentSample,
                                                   "No" = 1))
      geneData[[j]]$Sample[i] <- paste(currentSample, 1)
    }
    
    #if we haven't seen it before for any genes
    if (!(geneData[[j]]$Sample[i] %in% masterSample)) {
      masterSample <- c(masterSample, geneData[[j]]$Sample[i])
    }
  }
}

#Data Format ----
#What are the samples we're looking at?
masterSample <- sort(masterSample[masterSample != ""]) #ignore blank entries

#Format all data into one data.frame
ctData <- data.frame(Sample = masterSample, stringsAsFactors = FALSE)
for (j in 1:length(geneData)) {
  ctData[names(geneData)[j]] <- rep(NA, times = length(masterSample))
  for (i in 1:length(masterSample)) {
    n <- geneData[[j]]$CT[geneData[[j]]$Sample == masterSample[i]]
    if (length(n) == 1) ctData[i,j+1] <- n #add CT value if it exists
  }
}

#Parse Sample Names ----
#Experimental vs Control Genotype
ctData["Genotype"] <- ifelse(
  test = grepl(pattern = "TH/+",
               x = ctData$Sample,
               fixed = TRUE),
  yes = "TH/+",
  no = ifelse(test = grepl(pattern = "TH/SOD2",
                           x = ctData$Sample,
                           fixed = TRUE),
              yes = "TH/SOD2",
              no = ctData$Sample)
)

#Gravity Level
ctData["Gravity"] <- ifelse(
  test = grepl(pattern = "1g",
               x = ctData$Sample,
               fixed = TRUE),
  yes = 1,
  no = ifelse(test = grepl(pattern = "1.2g",
                           x = ctData$Sample,
                           fixed = TRUE),
              yes = 1.2,
              no = ifelse(test = grepl(pattern = "3g",
                                       x = ctData$Sample,
                                       fixed = TRUE),
                          yes = 3,
                          no = NA))
)

#Sex
ctData["Sex"] <- ifelse(
  test = grepl(pattern = "M",
               x = ctData$Sample,
               fixed = TRUE),
  yes = "M",
  no = ifelse(test = grepl(pattern = "F",
                           x = ctData$Sample,
                           fixed = TRUE),
              yes = "F",
              no = "")
)
#Subset - just grab the last character of sample name
ctData["Subset"] <- sapply(X = ctData$Sample, FUN = function(x) {
  return(as.numeric(substr(x = x, start = nchar(x), stop = nchar(x))))
})

#Error correction ----
#Looking at the data, it's obvious that TH/SOD2 M 1.2g 1 is unreasonable
ctData <- ctData[-c(13),]

#delta-CT: 'GAPDH-new' ----
dct.GAPDHnew <- data.frame("Genotype" = ctData$Genotype,
                           "Gravity" = ctData$Gravity,
                           "Sex" = ctData$Sex,
                           "Subset" = ctData$Subset)

#which columns need compared to housekeeping gene?
columnCompare <- which(!(colnames(ctData) %in% c("Sample", "Genotype",
                                                 "Gravity", "Sex",
                                                 "Subset", "Gapdh")))
for (j in columnCompare) {
  dct.GAPDHnew[colnames(ctData)[j]] <- ctData[[j]] - ctData[["GAPDH-new"]]
}

#delta-delta CT ----
ddct.GAPDHnew <- dct.GAPDHnew #start with a copy of delta-CT

for (i in 1:nrow(ddct.GAPDHnew)) {
  for (j in 5:ncol(dct.GAPDHnew)) { #skip the "genotype", "gravity",
                                    #"sex", and "subset" columns
    controlValue <- mean(subset(dct.GAPDHnew,
                                dct.GAPDHnew$Genotype == "TH/+"
                                #& dct.GAPDHnew$Sex == "F"
                                & dct.GAPDHnew$Gravity == 1)[,j],
                         na.rm = TRUE)
    ddct.GAPDHnew[i,j] <- dct.GAPDHnew[i,j] - controlValue
    ddct.GAPDHnew[i,j] <- 2^(-ddct.GAPDHnew[i,j]) #convert to fold-change scale
  }
}

#quick overall plots
pdf(file = "qPCR.pdf", paper = "USr", width = 10, height = 7)
omar <- par("mar")
par(mar = c(7.6, omar[-1]))
for (j in 5:ncol(ddct.GAPDHnew)) {
  plotData <- ddct.GAPDHnew[order(ddct.GAPDHnew$Sex, #get everything in order
                                  ddct.GAPDHnew$Genotype,
                                  ddct.GAPDHnew$Gravity,
                                  ddct.GAPDHnew$Subset),]
  barplot(height = plotData[,j], names.arg = paste0(plotData$Genotype, " ",
                                                    plotData$Sex, " ",
                                                    plotData$Gravity, "g",
                                                    " ", plotData$Subset),
          main = colnames(plotData)[j], las = 2)
}
par(mar = omar)
dev.off()

#Careful plots ----
#SOD2
newPlotData <- plotData[-c(33,42),]
par(mar = c(7.2, omar[-1]))
barplot(newPlotData$SOD2, names.arg = paste0(newPlotData$Genotype, " ",
                                          newPlotData$Sex, " ",
                                          newPlotData$Gravity, "g",
                                          " ", newPlotData$Subset),
        main = "SOD2", las = 2, ylab = "Relative mRNA Expression")

newPlotData$Group <- paste(newPlotData$Genotype,
                           newPlotData$Gravity
                           #newPlotData$Sex
                           )
avgData <- data.frame("Sample" = NULL, "Average" = NULL, "SEM" = NULL)
for (group in unique(newPlotData$Group)) {
  currentSample <- newPlotData$SOD2[newPlotData$Group == group]
  avgData <- rbind(avgData,
                   data.frame("Sample" = group,
                              "Average" = mean(currentSample,
                                               na.rm = TRUE),
                              "SEM" = sd(currentSample, na.rm = TRUE)/
                                sqrt(length(currentSample))))
}
n <- barplot(avgData$Average, names.arg = avgData$Sample,
             col = rainbow(3, s=0.5), main = "SOD2 Relative mRNA Expression",
             las = 2, ylim = c(0,5))
arrows(x0 = n, y0 = avgData$Average, y1 = (avgData$Average + 2*avgData$SEM),
       angle = 90, length = 0.1)
arrows(x0 = n, y0 = avgData$Average, y1 = (avgData$Average - 2*avgData$SEM),
       angle = 90, length = 0.1)
legend("topright", legend = c("1g", "1.2g", "3g"), fill = rainbow(3, s=0.5))
par(mar = omar)

#Heatmaps ----
#Make a nice color palette
greenramp <- colorRampPalette(colors = c("#E6FFF1", "#63D297", "#084D28"))
directoryName <- paste(substr(directoryName,
                              start = 1,
                              stop = nchar(directoryName)-3))
oxidative <- read.csv(file = paste0(directoryName, "Oxidative.csv"),
                      row.names = 1)
oxidative <- t(scale(as.matrix(oxidative)))
rownames(oxidative) <- c(rownames(oxidative)[1:2], "DmDRP",
                         rownames(oxidative)[4:7])
require("lattice")
print(levelplot(oxidative, xlab = "Oxidative Stress Genes", ylab = "",
                col.regions = greenramp))

downstream <- read.csv(file = paste0(directoryName, "Downstream.csv"),
                       row.names = 1)
downstream <- t(scale(as.matrix(downstream)))
print(levelplot(downstream, xlab = "Downstream Genes", ylab = "",
                col.regions = greenramp))
      
