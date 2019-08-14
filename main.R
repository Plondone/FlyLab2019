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
}

#Data Format ----
#How many different samples are we working with?
sampleVec <- c()
for (i in 1:length(fileVec)) {
  sampleVec[i] <- length(unique(rawData[[i]]$Sample))
}

#What are the sample names?
sampleNames <- unique(rawData[[which.max(sampleVec)]]$Sample)
sampleNames <- sampleNames[sampleNames != ""] #ignore blank entries
sampleNames <- sampleNames[sampleNames != "Water"] #ignore water control

#Format all data into one data.frame
ctData <- data.frame(Sample = sort(sampleNames), stringsAsFactors = FALSE)
for (i in 1:length(fileVec)) {
  geneNames <- unique(rawData[[i]]$Gene)
  for (j in 1:length(geneNames)) {
    ctData[geneNames[j]] <- sapply(
      X = sampleNames,
      FUN = function(x) {
        output <- mean(
          rawData[[i]]$CT[rawData[[i]]$Sample == x
                          & rawData[[i]]$Gene == geneNames[j]],
          na.rm = TRUE)
        if (!is.nan(output)) return(output) else return(NA)
      }
    )
  }
}

#Parse Sample Names ----
#Experimental vs Control Genotype
ctData["Genotype"] <- ifelse(
  test = grepl(pattern = "+",
               x = ctData$Sample,
               fixed = TRUE),
  yes = "TH/+",
  no = ifelse(test = grepl(pattern = "SOD2",
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

#delta-CT: 'GAPDH-new' ----
dct.GAPDHnew <- data.frame("Genotype" = ctData$Genotype,
                           "Gravity" = ctData$Gravity,
                           "Sex" = ctData$Sex)

#which columns need compared to housekeeping gene?
columnCompare <- which(!(colnames(ctData) %in% c("Sample", "Genotype",
                                                 "Gravity", "Sex")))
for (j in columnCompare) {
  dct.GAPDHnew[colnames(ctData)[j]] <-
    ctData[[j]] - ctData[["GAPDH-new"]]
}

#delta-delta CT
ddct.GAPDHnew <- dct.GAPDHnew #start with a copy of delta-CT

for (i in 1:nrow(ddct.GAPDHnew)) {
  for (j in 4:ncol(dct.GAPDHnew)) { #skip the "genotype", "gravity", "sex" col
    controlValue <- subset(dct.GAPDHnew,
                           dct.GAPDHnew$Genotype == dct.GAPDHnew$Genotype[i]
                           & dct.GAPDHnew$Sex == dct.GAPDHnew$Sex[i]
                           & dct.GAPDHnew$Gravity == 1)[,j]
    ddct.GAPDHnew[i,j] <- dct.GAPDHnew[i,j] - controlValue
    ddct.GAPDHnew[i,j] <- 2^(-ddct.GAPDHnew[i,j]) #convert to fold-change scale
  }
}

#make some plots!
for (j in 4:ncol(ddct.GAPDHnew)) {
  
}





