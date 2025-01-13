library(spatstat)
library(geosphere)
source("code/palaeoFunctions.R")
#feed in data and allocate to KLB according to zones
data.pri<-read.csv("data/PRI_Cretaceous_site_by_taxa_Oct30.csv")

data.klb9<-data.pri[data.pri$KLB==9,]

#read in Huw's Shapes
line4<-read.delim("data/Line4.txt")
line5<-read.delim("data/Line5.txt")

#line1<-line1[,c(4,5)]
line4<-line4[,c(4,5)]
line5<-line5[,c(5,4)]

#make KLB 9 window --------
rev.line5 <-cbind(rev(line5[,1]),rev(line5[,2]))
colnames(rev.line5)<-colnames(line5[,1:2])

coords.klb9<-rbind(rev.line5,(line4[,1:2]))
coords.klb9.closed <- rbind(coords.klb9, line5[nrow(line5),])

plot(coords.klb9, type = "l")
# get stations using spatstat

win.klb9<-owin(poly=list(x=rev(coords.klb9$Long),y=rev(coords.klb9$Lat)))
ppp.klb9<-ppp(data.pri$Longitude, data.pri$Latitude,
              marks= data.pri$station, window = win.klb9)

klb9.st.ppp <- ppp.klb9$marks


#treatment 1 - widen KLB ----

#move bottom of klb 
res<-move.line(line4,0.007) #have dist shows 579m
left.shift.bottom <- get.new.line(ori.line = line4, moved.line = res[[2]])

coords.klb9.treat1 <-rbind(rev.line5,(left.shift.bottom[,1:2]))
coords.klb9.treat1.closed <- rbind(coords.klb9.treat1, line5[nrow(line5),])

coords.klb9.treat1 <- as.data.frame(coords.klb9.treat1)
colnames(coords.klb9.treat1) <- c("Long", "Lat")

win.klb9.treat1 <- owin(poly=list(x=rev(coords.klb9.treat1$Long),
                                  y=rev(coords.klb9.treat1$Lat)))
ppp.klb9.treat1 <- ppp(data.pri$Longitude, data.pri$Latitude,
                       marks= data.pri$station, window = win.klb9.treat1)

klb9.st.treat1 <- ppp.klb9.treat1$marks

#treatment 2 - contract KLB ----

#move bottom of klb 
res<-move.line(line4,0.007) #have dist shows 579m
right.shift.bottom <- get.new.line(ori.line = line4, moved.line = res[[1]])

coords.klb9.treat2 <-rbind(rev.line5,(right.shift.bottom[,1:2]))
coords.klb9.treat2.closed <- rbind(coords.klb9.treat2, line5[nrow(line5),])

coords.klb9.treat2 <- as.data.frame(coords.klb9.treat2)
colnames(coords.klb9.treat2) <- c("Long", "Lat")

win.klb9.treat2 <- owin(poly=list(x=rev(coords.klb9.treat2$Long),
                                  y=rev(coords.klb9.treat2$Lat)))
ppp.klb9.treat2 <- ppp(data.pri$Longitude, data.pri$Latitude,
                       marks= data.pri$station, window = win.klb9.treat2)

klb9.st.treat2 <- ppp.klb9.treat2$marks


# create sensitivity data frames ------------------------------------------

#raw
rownames(data.klb9) <- data.klb9$station
data.klb9 <- data.klb9[12:ncol(data.klb9)]
data.klb9_1 <- data.klb9[sample(nrow(data.klb9), size = 40, replace = FALSE), ]
data.klb9_1  <- remove.zeroes(data.klb9_1 )
turn.res1 <- Turnover(data.klb9_1)


data.klb9_1 <- data.klb9[sample(nrow(data.klb9), size = 40, replace = FALSE), ]
data.klb9_1  <- remove.zeroes(data.klb9_1 )
turn.res1 <- Turnover(data.klb9_1)


data.klb9_1 <- data.klb9[sample(nrow(data.klb9), size = 40, replace = FALSE), ]
data.klb9_1  <- remove.zeroes(data.klb9_1 )
turn.res1 <- Turnover(data.klb9_1)

#exclude nektonic
data.klb9 <- remove.zeroes(data.klb9)

library(metacom)
turn.res <- Turnover(data.klb9)


remove.zeroes(data.klb9)
#treatment 1
data.klb9.treat1 <- data.pri[data.pri$station %in% klb9.st.treat1,]
rownames(data.klb9.treat1) <- klb9.st.treat1
#exclude nektonic 
data.klb9.treat1 <- data.klb9.treat1[,12:ncol(data.klb9.treat1)]

#treatment 2 
data.klb9.treat2 <- data.pri[data.pri$station %in% klb9.st.treat2,]
rownames(data.klb9.treat2) <- klb9.st.treat2
#exclude nektonic 
data.klb9.treat2 <- data.klb9.treat2[,12:ncol(data.klb9.treat2)]


# Save stations -----------

save(data.pri, klb9.col, data.klb9, klb9.st.ppp, klb9.st.treat1, 
     klb9.st.treat2, data.klb9.treat1, data.klb9.treat2,
     file = "results/sensitivity_klb9.RData")

# Open here analyses ------------------------------------------------------

load(file = "results/sensitivity_klb9.RData") 

library(metacom)
library(cooccur)
source("code/palaeoFunctions.R")

data.klb9 <- data.klb9[,-c(48, 56, 57)]
data.klb9 <- remove.zeroes(data.klb9)

data.klb9.treat1 <- data.klb9.treat1[,-c(48, 56, 57)]
data.klb9.treat1 <- remove.zeroes(data.klb9.treat1)

data.klb9.treat2 <- data.klb9.treat2[,-c(48, 56, 57)]
data.klb9.treat2 <- remove.zeroes(data.klb9.treat2)

#make list of KLB 9 res 

klb9.metacom.res <- vector(mode = "list", length = 3)
names(klb9.metacom.res) <- c("raw", "treat1", "treat2")

klb9.metacom.res$raw <- get.metacom.res(data.klb9)
klb9.metacom.res$treat1 <- get.metacom.res(data.klb9.treat1)
klb9.metacom.res$treat2 <- get.metacom.res(data.klb9.treat2)

klb9.metacom.res.df <- do.call(cbind, klb9.metacom.res)
klb9.metacom.res.df2 <- do.call(rbind, klb9.metacom.res)

write.csv(klb9.metacom.res.df, 
          file = "results/klb9_metacom_res.csv",
          row.names = FALSE)

klb9.metacom.outputs <- metacom.output(klb9.metacom.res)

save(klb9.metacom.res, klb9.metacom.outputs, klb9.metacom.res.df,
     file = "results/klb9.metacom.results.RData")


# do co-occur -------------------------------------------------------------

library(cooccur)
klb9.cooc.res <- vector("list", length = 3)
names(klb9.cooc.res) <- c("raw", "expand", "contract")

data.klb9 <- remove.sing(data.klb9)

data.klb9 <- as.data.frame(t(data.klb9)) #transpose the matrix

klb9.mat <- cooccur(data.klb9, type = "spp_site", spp_names = TRUE, 
                    thresh = TRUE)
summary(klb9.mat)
plot(klb9.mat)

klb9.cooc.res[[1]] <- klb9.mat

data.klb9.treat1 <- remove.sing(data.klb9.treat1)
data.klb9.treat1 <- as.data.frame(t(data.klb9.treat1)) #transpose the matrix
klb9.treat1.mat <- cooccur(data.klb9.treat1, type = "spp_site", spp_names = TRUE, 
                           thresh = TRUE)
plot(klb9.treat1.mat)
summary(klb9.treat1.mat)
klb9.cooc.res[[2]] <- klb9.treat1.mat

data.klb9.treat2 <- remove.sing(data.klb9.treat2)
data.klb9.treat2 <- as.data.frame(t(data.klb9.treat2)) #transpose the matrix
klb9.treat2.mat <- cooccur(data.klb9.treat2, type = "spp_site", spp_names = TRUE, 
                           thresh = TRUE)
plot(klb9.treat2.mat)
summary(klb9.treat2.mat)

klb9.cooc.res[[3]] <- klb9.treat2.mat

save(klb9.cooc.res, file = "results/klb9.cooccur.res.RData")



