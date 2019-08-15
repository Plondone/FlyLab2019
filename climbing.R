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
climbData$Group <- paste(climbData$Genotype,
                         #climbData$Sex,
                         climbData$Gravity)

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
plotAverages < rbind(plotAverages[1,],
                     plotAverages[3,],
                     plotAverages[2,],
                     plotAverages[4,])
plotErr <- rbind(sqrt(plotVar[1,])/sqrt(138), #TH/+ 1g total n
                 sqrt(plotVar[3,])/sqrt(132), #TH/SOD2 1g total n
                 sqrt(plotVar[2,])/sqrt(94), #TH/+ 3g total n
                 sqrt(plotVar[4,])/sqrt(129)) #TH/SOD2 3g total n

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

colorVec <- rainbow(nrow(plotAverages), s = 0.5)

errorPlot(x = 1:10, y = plotAverages, stdev = 2*plotErr,
          col = colorVec, type = "l",
          ylim = c(0.4, 1.1),
          main = "Climbing Assay",
          xlab = "Day", ylab = "Proportion Passed")
legend(x = "bottomleft", legend = paste0(rownames(plotAverages), "g"),
       fill = colorVec)

#*Barplots
omar <- par("mar")
#par(mar = c(6.9, omar[-1]))

for (j in 1:ncol(plotAverages)) {
  n1 <- barplot(height = plotAverages[,j],
                names.arg = paste0(rownames(plotAverages), "g"),
                #las = 2,
                col = rainbow(4, s = 0.5), ylab = "Proportion Passed",
                main = paste("Climbing Assay, Day", j), ylim = c(0,1))
  arrows(x0 = n1, y0 = plotAverages[,j], angle = 90, length = 0.1,
         y1 = plotAverages[,j] + 2*plotErr[,j])
}
par(mar = omar)



