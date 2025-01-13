#feed in data and allocate to KLB according to zones
data.pri<-read.csv("data/PRI_Cretaceous_site_by_taxa_Oct30.csv")

data.klb7<-data.pri[data.pri$KLB==7,]
data.klb8 <- data.pri[data.pri$KLB==8,]
data.klb9 <- data.pri[data.pri$KLB==9,]

library(dplyr)
library(plyr)
df7.8.9 <- rbind.fill(data.klb7, data.klb8, data.klb9)
rownames(df7.8.9) <- df7.8.9$station
df7.8.9 <- df7.8.9[,-c(1:12, 59, 67, 68)]
df7.8.9[is.na(df7.8.9)] <- 0 
df7.8.9 <- remove.zeroes(df7.8.9)

source("code/palaeoFunctions.R")
df <- remove.zeroes(df7.8.9)

library(metacom)
df2 <- OrderMatrix(df, scores =1, outputScores = FALSE)
metacom::Imagine(df2, col = c(0, "grey"),fill=FALSE,
                 xlab = "", ylab = "", speciesnames = FALSE,
                 sitenames = FALSE, order = TRUE)
#add taxa names
axis(3, at = c(1:ncol(df2)), labels = FALSE)
axis(side = 3, las = 2, mgp = c(3, 0.75, 0), labels = FALSE)
text(x = 1:ncol(df2), y = par("usr")[2] + 160, 
     labels = colnames(df2), 
     srt = -35, cex = 0.8, xpd = NA, adj = 1)
axis(side = 2, at = c(1:nrow(df2)), labels= FALSE)

#colour part of the labels 
all.sites <- rownames(df2) #ordinated station names

#klb9 
klb9.st <- data.klb9$station
ind.9 <- match(klb9.st, all.sites)
axis(side = 2, at = ind.9[1:length(ind.9)], labels = FALSE, 
     col.ticks = "orange", lwd = 2)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0), labels = FALSE, 
     cex.axis = 0.5)
text(y = ind.9[1:length(ind.9)], x = par("usr")[1] - 1, 
     labels = klb9.st, xpd = NA, cex = 0.5, col = "orange")


#klb7
klb7.st <- data.klb7$station
ind.7 <- match(klb7.st, all.sites)
axis(side = 2, at = ind.7[1:length(ind.7)], labels = FALSE, 
     col.ticks = "forestgreen", lwd = 2)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0), labels = FALSE, 
     cex.axis = 0.5)
text(y = ind.7[1:length(ind.7)], x = par("usr")[1] - 1, 
     labels = klb7.st, xpd = NA, cex = 0.5, col = "forestgreen")

#klb8
klb8.st <- data.klb8$station
ind.8 <- match(klb8.st, all.sites)
axis(side = 2, at = ind.8[1:length(ind.8)], labels = FALSE, 
     col.ticks = "blue", lwd = 2)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0), labels = FALSE, 
     cex.axis = 0.5)
text(y = ind.8[1:length(ind.8)], x = par("usr")[1] - 1, 
     labels = klb8.st, xpd = NA, cex = 0.5, col = "blue")


legend("bottomright", legend = c("ORDINATED", "klb7", 
                                 "klb8", "klb9"),col = c("black", "forestgreen",
                                                         "blue", "orange"),
       bg = "white",lwd = 2, lty = c(0, 1, 1, 1))