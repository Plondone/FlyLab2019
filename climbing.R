#TJ Wiegman, SLSTP Research Associate
#Ames Research Center
#Bhattacharya Lab, Summer 2019
#R version 3.6.1

#Data Entry ----
directoryName <- "/Users/TJ Wiegman/Google Drive/School/SU19/Climbing Results"
data1 <- read.csv(file = paste0(directoryName, "/1g.csv"),
                  stringsAsFactors = FALSE, skip = 2)
data3 <- read.csv(file = paste0(directoryName, "/3g.csv"),
                  stringsAsFactors = FALSE, skip = 2)

#Data Formatting ----
columnNames <- c(sort(paste0("Day", 1:9, "n", rep(1:3, each = 9))),
              "Day10n1", "Day10n2", "Day10n3") #assuming 30 col of data
climb1 <- data.frame("Genotype" = data1$Tape, "Gravity" = 1,
                     "Sex" = data1$Vial, "Group" = NA,
                     stringsAsFactors = FALSE)
climb3 <- data.frame("Genotype" = data3$Tape, "Gravity" = 3,
                     "Sex" = data3$Vial, "Group" = NA,
                     stringsAsFactors = FALSE)

#assuming data starts on column 4 and goes to edge
for (j in 4:ncol(data1)) climb1[columnNames[j-3]] <- data1[,j]
for (j in 4:ncol(data3)) climb3[columnNames[j-3]] <- data3[,j]
rm(list = c("data1", "data3"))

#Get everything into one data.frame
climbData <- rbind(climb1, climb3)
climbData$Sex <- sapply(X = climbData$Genotype,
                        FUN = function(x) {
                          switch(x, "Blue" = "M", "Green" = "M",
                                 "Orange" = "F", "Pink" = "F")
                        })
climbData$Genotype <- sapply(X = climbData$Genotype,
                             FUN = function(x) {
                               switch(x, "Blue" = "TH/SOD2", "Green" = "TH/+",
                                      "Orange" = "TH/+", "Pink" = "TH/SOD2")
                             })
climbData$Group <- paste(climbData$Genotype, climbData$Sex, climbData$Gravity)

#Data Analysis ----
#*Daily Trends ----
dayAverages <- climbData[,c(1:4)]
dayVar <- climbData[,c(1:4)]
columnNames <- paste0("Day", 1:10)
for (j in 1:10) {
  n1 <- climbData[,3*j +2]
  n2 <- climbData[,3*j +3]
  n3 <- climbData[,3*j +4]
  data1 <- cbind(n1, n2, n3) #combine day into one dataframe, then apply by row
  dayAverages[columnNames[j]] <- apply(data1, MARGIN = 1, FUN = mean)
  dayVar[columnNames[j]] <- apply(data1, MARGIN = 1, FUN = var)
}

#get the avg pass-rate and associated variance for each group
plotAverages <- as.data.frame(apply(X = dayAverages[5:ncol(dayAverages)],
                                    MARGIN = 2,
                                    FUN = function(x) {
                                      by(x, INDICES = dayAverages$Group,
                                         FUN = mean, na.rm = TRUE)
                                    }))
plotVar <- as.data.frame(apply(X = dayVar[5:ncol(dayVar)],
                               MARGIN = 2,
                               FUN = function(x) {
                                 by(x, INDICES = dayVar$Group,
                                    FUN = mean, na.rm = TRUE)
                               }))

#Plot daily trends by group (incl. handy wrapper function)
errorPlot <- function(x, y, stdev, col, rows = 1:nrow(y),
                      xlim = c(min(x), max(x)), ylim = c(min(y), max(y)),
                      type = "p", xlab = "", ylab = "", main = "") {
  #x is a vector of x-values to be used by all sets of y-values
  #y is a matrix, where each row is to be plotted as a set of y-values
  #stdev is a matrix of same dimension as y, assumed to be symmetric
  #col is a vector of color values
  arrowAng <- 90
  arrowLen <- 0.1
  i <- min(rows)
  
  #Plot first set of data; lower, and upper error bars
  plot(x = x, y = y[i,], xlim = xlim, ylim = ylim, col = col[i],
       main = main, xlab = xlab, ylab = ylab, type = type)
  arrows(x0 = x, y0 = as.numeric(y[i,]), y1 = as.numeric(y[i,] - stdev[i,]),
         angle = arrowAng, length = arrowLen, col = col[i])
  arrows(x0 = x, y0 = as.numeric(y[i,]), y1 = as.numeric(y[i,] + stdev[i,]),
         angle = arrowAng, length = arrowLen, col = col[i])
  
  #Add rest of data and error bars to plot
  for (i in rows[-1]) {
    points(x = x, y = y[i,], col = col[i], type = type)
    arrows(x0 = x, y0 = as.numeric(y[i,]), y1 = as.numeric(y[i,] - stdev[i,]),
           angle = arrowAng, length = arrowLen, col = col[i])
    arrows(x0 = x, y0 = as.numeric(y[i,]), y1 = as.numeric(y[i,] + stdev[i,]),
           angle = arrowAng, length = arrowLen, col = col[i])
  }
  
}

colorVec <- rainbow(nrow(plotAverages))
fem <- c(1,2,5,6)
mal <- c(3,4,7,8)

errorPlot(x = 1:10, y = plotAverages, stdev = sqrt(plotVar),
          col = colorVec, type = "l", rows = fem,
          ylim = c(0, 1.1),
          main = "Climbing Assay, Females",
          xlab = "Day", ylab = "Proportion Passed")

legend(x = "bottomleft", legend = paste0(rownames(plotAverages)[fem], "g"),
       fill = colorVec[fem])

errorPlot(x = 1:10, y = plotAverages, stdev = sqrt(plotVar),
          col = colorVec, type = "l", rows = mal,
          ylim = c(0, 1.1),
          main = "Climbing Assay, Males",
          xlab = "Day", ylab = "Proportion Passed")
legend(x = "bottomleft", legend = paste0(rownames(plotAverages)[mal], "g"),
       fill = colorVec[mal])

#*Weekly Trends ----
weekAverages <- dayAverages[,c(1:4)] #start with same 4 label columns
weekAverages["Days1to5"] <- apply(X = dayAverages[,c(5:9)],
                                  MARGIN = 1, FUN = mean)
weekAverages["Days6to10"] <- apply(X = dayAverages[,c(10:14)],
                                   MARGIN = 1, FUN = mean)
weekVar <- dayVar[,c(1:4)]
weekVar["Days1to5"] <- apply(X = dayVar[,c(5:9)],
                                  MARGIN = 1, FUN = mean)
weekVar["Days6to10"] <- apply(X = dayVar[,c(10:14)],
                                   MARGIN = 1, FUN = mean)

plotAverages <- as.data.frame(apply(X = weekAverages[5:ncol(weekAverages)],
                                    MARGIN = 2,
                                    FUN = function(x) {
                                      by(x, INDICES = weekAverages$Group,
                                         FUN = mean, na.rm = TRUE)
                                    }))
plotVar <- as.data.frame(apply(X = weekVar[5:ncol(weekVar)],
                               MARGIN = 2,
                               FUN = function(x) {
                                 by(x, INDICES = weekVar$Group,
                                    FUN = mean, na.rm = TRUE)
                               }))
errorPlot(x = c(5,10), y = plotAverages, stdev = sqrt(plotVar), col = colorVec,
          xlim = c(1,10), ylim = c(0, 1.1), rows = fem,
          main = "Climbing Assay, Females, Pooled d1-5 and d6-10",
          ylab = "Proportion Passed")
legend(x = "bottomleft", legend = paste0(rownames(plotAverages)[fem], "g"),
       fill = colorVec[fem])



