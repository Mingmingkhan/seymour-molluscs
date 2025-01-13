source("code/palaeoFunctions.R")

modified.remove.zeroes <- function(df){
  #remove taxa with 0 counts
  c <- as.numeric(colSums(df))
  c.ind <- which(c == 0)
  if (length(c.ind > 0)){
    temp <- df[,-c.ind ]
  }
  r <- as.numeric(rowSums(df))
  r.ind <- which(r==0)
  if (length(r.ind > 0)){
    temp2 <- temp[-r.ind,]
  }
  else temp2 <- temp
  #temp2[temp2 >= 1] <- 1 #change to presence absence
  #df <- temp2
  return(temp2)
}

# Supplementary Figure 1: Rarefaction curves ----


data.pri<-read.csv("data/PRI_Cretaceous_site_by_taxa_Oct30.csv")

data.7 <- data.pri[data.pri$KLB == 7,]
#keep benthic only 
data.7 <- data.7[,12:66]
data.7 <- data.7[,-c(48, 52)]

data.8 <- data.pri[data.pri$KLB == 8,]
data.8 <- data.8[,12:66]
data.8 <- data.8[,-c(48, 52)]

data.9 <- data.pri[data.pri$KLB == 9,]
rownames(data.9) <- data.9$station
data.9 <- data.9[,12:66]
data.9 <- data.9[,-c(48, 52)]

data.789 <- rbind(data.7, data.8, data.9)
data.789 <- modified.remove.zeroes(data.789)

data.7 <- modified.remove.zeroes(data.7)
data.8 <- modified.remove.zeroes(data.8)
data.9 <- modified.remove.zeroes(data.9)

#do the rarefaction (Fig s1)
library(vegan)

sp.789 <- specaccum(data.789, method = "rarefaction")
sp.7 <- specaccum(data.7, method = "rarefaction")
sp.8 <- specaccum(data.8, method = "rarefaction")
sp.9 <- specaccum(data.9, method = "rarefaction")

#plot 
windows()
svg(file = "figures/specaccum.svg", h = 10, w = 10)
plot(sp.789, main = "Accumulation curves", xlab = "number of sites", 
     ylab = "observed richness (number of families)")
grid(col = "lightgray")
plot(sp.9, add = TRUE, col = "#E66100", lwd = 1.5)
plot(sp.7, add = TRUE, col = "#0C7BDC", lwd = 1.5)
plot(sp.8, add = TRUE, col = "#009E73", lwd = 1.5)
dev.off()


# spatial thinning -------------------------------------------------------

#first KLB 7 ----

data.7 <- data.pri[data.pri$KLB == 7,]
rownames(data.7) <- data.7$station
klb7.st <- data.7[1:3]
data.7 <- data.7[,12:66]
data.7 <- data.7[,-c(48, 52)]
data.7 <- modified.remove.zeroes(data.7)

#turn coords into an sf object 
library(sf)
library(tidysdm)

stations <- rownames(data.7)
stations <- as.numeric(stations)
klb7.st <- klb7.st[klb7.st$station %in% stations,]

data_coords <- st_as_sf(klb7.st, coords = c("Longitude", "Latitude"))
st_crs(data_coords) <- 4326

#klb7 thinning 30 ----

klb7_thinning30 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb7_thinning30)[i] <- paste("trial_",i)
  klb7_thinning30[[i]]$klb7_thinned_st <- thin_by_dist(data_coords, 30)
  klb7_thinned_dat <- data.pri[data.pri$station %in% klb7_thinning30[[i]]$klb7_thinned_st$station,]
  rownames(klb7_thinned_dat) <- klb7_thinned_dat$station
  #keep benthic only 
  klb7_thinned_dat <- klb7_thinned_dat[,12:66]
  klb7_thinned_dat <- klb7_thinned_dat[,-c(48, 52)]
  klb7_thinned_dat2 <- modified.remove.zeroes(klb7_thinned_dat)
  klb7_thinned_dat <- remove.sing.zeroes(klb7_thinned_dat)
  klb7_thinning30[[i]]$klb7_thinned_dat <- klb7_thinned_dat
  klb7_thinning30[[i]]$klb7_thinned_sites <- nrow(klb7_thinned_dat)
  #do metacom
  klb7_thinning30[[i]]$coh.value <- metacom::Coherence(klb7_thinned_dat)
  klb7_thinning30[[i]]$coh.z <- klb7_thinning30[[i]]$coh.value[2, 2]
  klb7_thinning30[[i]]$turn.value <- metacom::Turnover(klb7_thinned_dat)
  klb7_thinning30[[i]]$turn.z <- klb7_thinning30[[i]]$turn.value[2,2]
  klb7_thinning30[[i]]$bc.value <- metacom::BoundaryClump(klb7_thinned_dat)
  klb7_thinning30[[i]]$bc.i <- klb7_thinning30[[i]]$bc.value[1,2]
  klb7_thinning30[[i]]$bc.p <- klb7_thinning30[[i]]$bc.value[2,2]
  
  #do cooccur
  temp2 <- as.data.frame(t(klb7_thinned_dat))
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb7_thinning30[[i]]$cooc.mat <- cooc.mat
  klb7_thinning30[[i]]$cooc.percent <- temp.res
  klb7_thinning30[[i]]$sp_dat <- klb7_thinned_dat2
}


klb7_thinned.site.number <- lapply(klb7_thinning30, `[[`,3)
klb7_thinned.site.number <- do.call(rbind, klb7_thinned.site.number)
colnames(klb7_thinned.site.number) <- "sampled_sites"

coh.res.all <- lapply(klb7_thinning30, `[[`,5)
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb7_thinning30, `[[`,7)
turn.res.all <- do.call(rbind, turn.res.all)

bc.i.all <- lapply(klb7_thinning30, `[[`,9)
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb7_thinning30, `[[`,10)
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb7_thinning30, `[[`,12)
cooc.res.all <- do.call(rbind, cooc.res.all)
klb7_thinned_summary_res30 <- as.data.frame(cbind(klb7_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb7_thinned_summary_res30) <- c("sampled_sites", "coherence",
                                          "turnover", "bc.i", "bc.p", "cooccur")

# klb7 thinning 50 --------------------------------------------------------


klb7_thinning50 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb7_thinning50)[i] <- paste("trial_",i)
  klb7_thinning50[[i]]$klb7_thinned_st <- thin_by_dist(data_coords, 50)
  klb7_thinned_dat <- data.pri[data.pri$station %in% klb7_thinning50[[i]]$klb7_thinned_st$station,]
  rownames(klb7_thinned_dat) <- klb7_thinned_dat$station
  #keep benthic only 
  klb7_thinned_dat <- klb7_thinned_dat[,12:66]
  klb7_thinned_dat <- klb7_thinned_dat[,-c(48, 52)]
  klb7_thinned_dat2 <- modified.remove.zeroes(klb7_thinned_dat)
  klb7_thinned_dat <- remove.sing.zeroes(klb7_thinned_dat)
  klb7_thinning50[[i]]$klb7_thinned_dat <- klb7_thinned_dat
  klb7_thinning50[[i]]$klb7_thinned_sites <- nrow(klb7_thinned_dat)
  #do metacom
  klb7_thinning50[[i]]$coh.value <- metacom::Coherence(klb7_thinned_dat)
  klb7_thinning50[[i]]$coh.z <- klb7_thinning50[[i]]$coh.value[2, 2]
  klb7_thinning50[[i]]$turn.value <- metacom::Turnover(klb7_thinned_dat)
  klb7_thinning50[[i]]$turn.z <- klb7_thinning50[[i]]$turn.value[2,2]
  klb7_thinning50[[i]]$bc.value <- metacom::BoundaryClump(klb7_thinned_dat)
  klb7_thinning50[[i]]$bc.i <- klb7_thinning50[[i]]$bc.value[1,2]
  klb7_thinning50[[i]]$bc.p <- klb7_thinning50[[i]]$bc.value[2,2]
  
  #do cooccur
  temp2 <- as.data.frame(t(klb7_thinned_dat))
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb7_thinning50[[i]]$cooc.mat <- cooc.mat
  klb7_thinning50[[i]]$cooc.percent <- temp.res
  klb7_thinning50[[i]]$sp_dat <- klb7_thinned_dat2
}


klb7_thinned.site.number <- lapply(klb7_thinning50, `[[`,3)
klb7_thinned.site.number <- do.call(rbind, klb7_thinned.site.number)
colnames(klb7_thinned.site.number) <- "sampled_sites"

coh.res.all <- lapply(klb7_thinning50, `[[`,5)
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb7_thinning50, `[[`,7)
turn.res.all <- do.call(rbind, turn.res.all)

bc.i.all <- lapply(klb7_thinning50, `[[`,9)
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb7_thinning50, `[[`,10)
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb7_thinning50, `[[`,12)
cooc.res.all <- do.call(rbind, cooc.res.all)
klb7_thinned_summary_res50 <- as.data.frame(cbind(klb7_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb7_thinned_summary_res50) <- c("sampled_sites", "coherence",
                                          "turnover", "bc.i", "bc.p", "cooccur")






# klb7 thinning 70 --------------------------------------------------------


klb7_thinning70 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb7_thinning70)[i] <- paste("trial_",i)
  klb7_thinning70[[i]]$klb7_thinned_st <- thin_by_dist(data_coords, 70)
  klb7_thinned_dat <- data.pri[data.pri$station %in% klb7_thinning70[[i]]$klb7_thinned_st$station,]
  rownames(klb7_thinned_dat) <- klb7_thinned_dat$station
  #keep benthic only 
  klb7_thinned_dat <- klb7_thinned_dat[,12:66]
  klb7_thinned_dat <- klb7_thinned_dat[,-c(48, 52)]
  klb7_thinned_dat2 <- modified.remove.zeroes(klb7_thinned_dat)
  klb7_thinned_dat <- remove.sing.zeroes(klb7_thinned_dat)
  klb7_thinning70[[i]]$klb7_thinned_dat <- klb7_thinned_dat
  klb7_thinning70[[i]]$klb7_thinned_sites <- nrow(klb7_thinned_dat)
  #do metacom
  klb7_thinning70[[i]]$coh.value <- metacom::Coherence(klb7_thinned_dat)
  klb7_thinning70[[i]]$coh.z <- klb7_thinning70[[i]]$coh.value[2, 2]
  klb7_thinning70[[i]]$turn.value <- metacom::Turnover(klb7_thinned_dat)
  klb7_thinning70[[i]]$turn.z <- klb7_thinning70[[i]]$turn.value[2,2]
  klb7_thinning70[[i]]$bc.value <- metacom::BoundaryClump(klb7_thinned_dat)
  klb7_thinning70[[i]]$bc.i <- klb7_thinning70[[i]]$bc.value[1,2]
  klb7_thinning70[[i]]$bc.p <- klb7_thinning70[[i]]$bc.value[2,2]
  
  #do cooccur
  temp2 <- as.data.frame(t(klb7_thinned_dat))
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb7_thinning70[[i]]$cooc.mat <- cooc.mat
  klb7_thinning70[[i]]$cooc.percent <- temp.res
  klb7_thinning70[[i]]$sp_dat <- klb7_thinned_dat2
}


klb7_thinned.site.number <- lapply(klb7_thinning70, `[[`,3)
klb7_thinned.site.number <- do.call(rbind, klb7_thinned.site.number)
colnames(klb7_thinned.site.number) <- "sampled_sites"

coh.res.all <- lapply(klb7_thinning70, `[[`,5)
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb7_thinning70, `[[`,7)
turn.res.all <- do.call(rbind, turn.res.all)

bc.i.all <- lapply(klb7_thinning70, `[[`,9)
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb7_thinning70, `[[`,10)
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb7_thinning70, `[[`,12)
cooc.res.all <- do.call(rbind, cooc.res.all)
klb7_thinned_summary_res70 <- as.data.frame(cbind(klb7_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb7_thinned_summary_res70) <- c("sampled_sites", "coherence",
                                          "turnover", "bc.i", "bc.p", "cooccur")




# NOW KLB 8 ---------------------------------------------------------------

data.8 <- data.pri[data.pri$KLB == 8,]
rownames(data.8) <- data.8$station
klb8.st <- data.8[1:3]
data.8 <- data.8[,12:66]
data.8 <- data.8[,-c(48, 52)]
data.8 <- modified.remove.zeroes(data.8)

#turn coords into an sf object 
library(sf)
library(tidysdm)

stations <- rownames(data.8)
stations <- as.numeric(stations)
klb8.st <- klb8.st[klb8.st$station %in% stations,]

data_coords <- st_as_sf(klb8.st, coords = c("Longitude", "Latitude"))
st_crs(data_coords) <- 4326

#klb8 thinning 30 ----

klb8_thinning30 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb8_thinning30)[i] <- paste("trial_",i)
  klb8_thinning30[[i]]$klb8_thinned_st <- thin_by_dist(data_coords, 30)
  klb8_thinned_dat <- data.pri[data.pri$station %in% klb8_thinning30[[i]]$klb8_thinned_st$station,]
  rownames(klb8_thinned_dat) <- klb8_thinned_dat$station
  #keep benthic only 
  klb8_thinned_dat <- klb8_thinned_dat[,12:66]
  klb8_thinned_dat <- klb8_thinned_dat[,-c(48, 52)]
  klb8_thinned_dat2 <- modified.remove.zeroes(klb8_thinned_dat)
  klb8_thinned_dat <- remove.sing.zeroes(klb8_thinned_dat)
  klb8_thinning30[[i]]$klb8_thinned_dat <- klb8_thinned_dat
  klb8_thinning30[[i]]$klb8_thinned_sites <- nrow(klb8_thinned_dat)
  #do metacom
  klb8_thinning30[[i]]$coh.value <- metacom::Coherence(klb8_thinned_dat)
  klb8_thinning30[[i]]$coh.z <- klb8_thinning30[[i]]$coh.value[2, 2]
  klb8_thinning30[[i]]$turn.value <- metacom::Turnover(klb8_thinned_dat)
  klb8_thinning30[[i]]$turn.z <- klb8_thinning30[[i]]$turn.value[2,2]
  klb8_thinning30[[i]]$bc.value <- metacom::BoundaryClump(klb8_thinned_dat)
  klb8_thinning30[[i]]$bc.i <- klb8_thinning30[[i]]$bc.value[1,2]
  klb8_thinning30[[i]]$bc.p <- klb8_thinning30[[i]]$bc.value[2,2]
  
  #do cooccur
  temp2 <- as.data.frame(t(klb8_thinned_dat))
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb8_thinning30[[i]]$cooc.mat <- cooc.mat
  klb8_thinning30[[i]]$cooc.percent <- temp.res
  klb8_thinning30[[i]]$sp_dat <- klb8_thinned_dat2
}


klb8_thinned.site.number <- lapply(klb8_thinning30, `[[`,3)
klb8_thinned.site.number <- do.call(rbind, klb8_thinned.site.number)
colnames(klb8_thinned.site.number) <- "sampled_sites"

coh.res.all <- lapply(klb8_thinning30, `[[`,5)
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb8_thinning30, `[[`,7)
turn.res.all <- do.call(rbind, turn.res.all)

bc.i.all <- lapply(klb8_thinning30, `[[`,9)
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb8_thinning30, `[[`,10)
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb8_thinning30, `[[`,12)
cooc.res.all <- do.call(rbind, cooc.res.all)
klb8_thinned_summary_res30 <- as.data.frame(cbind(klb8_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb8_thinned_summary_res30) <- c("sampled_sites", "coherence",
                                          "turnover", "bc.i", "bc.p", "cooccur")

# klb8 thinning 50 --------------------------------------------------------


klb8_thinning50 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb8_thinning50)[i] <- paste("trial_",i)
  klb8_thinning50[[i]]$klb8_thinned_st <- thin_by_dist(data_coords, 50)
  klb8_thinned_dat <- data.pri[data.pri$station %in% klb8_thinning50[[i]]$klb8_thinned_st$station,]
  rownames(klb8_thinned_dat) <- klb8_thinned_dat$station
  #keep benthic only 
  klb8_thinned_dat <- klb8_thinned_dat[,12:66]
  klb8_thinned_dat <- klb8_thinned_dat[,-c(48, 52)]
  klb8_thinned_dat2 <- modified.remove.zeroes(klb8_thinned_dat)
  klb8_thinned_dat <- remove.sing.zeroes(klb8_thinned_dat)
  klb8_thinning50[[i]]$klb8_thinned_dat <- klb8_thinned_dat
  klb8_thinning50[[i]]$klb8_thinned_sites <- nrow(klb8_thinned_dat)
  #do metacom
  klb8_thinning50[[i]]$coh.value <- metacom::Coherence(klb8_thinned_dat)
  klb8_thinning50[[i]]$coh.z <- klb8_thinning50[[i]]$coh.value[2, 2]
  klb8_thinning50[[i]]$turn.value <- metacom::Turnover(klb8_thinned_dat)
  klb8_thinning50[[i]]$turn.z <- klb8_thinning50[[i]]$turn.value[2,2]
  klb8_thinning50[[i]]$bc.value <- metacom::BoundaryClump(klb8_thinned_dat)
  klb8_thinning50[[i]]$bc.i <- klb8_thinning50[[i]]$bc.value[1,2]
  klb8_thinning50[[i]]$bc.p <- klb8_thinning50[[i]]$bc.value[2,2]
  
  #do cooccur
  temp2 <- as.data.frame(t(klb8_thinned_dat))
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb8_thinning50[[i]]$cooc.mat <- cooc.mat
  klb8_thinning50[[i]]$cooc.percent <- temp.res
  klb8_thinning50[[i]]$sp_dat <- klb8_thinned_dat2
}


klb8_thinned.site.number <- lapply(klb8_thinning50, `[[`,3)
klb8_thinned.site.number <- do.call(rbind, klb8_thinned.site.number)
colnames(klb8_thinned.site.number) <- "sampled_sites"

coh.res.all <- lapply(klb8_thinning50, `[[`,5)
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb8_thinning50, `[[`,7)
turn.res.all <- do.call(rbind, turn.res.all)

bc.i.all <- lapply(klb8_thinning50, `[[`,9)
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb8_thinning50, `[[`,10)
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb8_thinning50, `[[`,12)
cooc.res.all <- do.call(rbind, cooc.res.all)
klb8_thinned_summary_res50 <- as.data.frame(cbind(klb8_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb8_thinned_summary_res50) <- c("sampled_sites", "coherence",
                                          "turnover", "bc.i", "bc.p", "cooccur")






# klb8 thinning 70 --------------------------------------------------------


klb8_thinning70 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb8_thinning70)[i] <- paste("trial_",i)
  klb8_thinning70[[i]]$klb8_thinned_st <- thin_by_dist(data_coords, 70)
  klb8_thinned_dat <- data.pri[data.pri$station %in% klb8_thinning70[[i]]$klb8_thinned_st$station,]
  rownames(klb8_thinned_dat) <- klb8_thinned_dat$station
  #keep benthic only 
  klb8_thinned_dat <- klb8_thinned_dat[,12:66]
  klb8_thinned_dat <- klb8_thinned_dat[,-c(48, 52)]
  klb8_thinned_dat2 <- modified.remove.zeroes(klb8_thinned_dat)
  klb8_thinned_dat <- remove.sing.zeroes(klb8_thinned_dat)
  klb8_thinning70[[i]]$klb8_thinned_dat <- klb8_thinned_dat
  klb8_thinning70[[i]]$klb8_thinned_sites <- nrow(klb8_thinned_dat)
  #do metacom
  klb8_thinning70[[i]]$coh.value <- metacom::Coherence(klb8_thinned_dat)
  klb8_thinning70[[i]]$coh.z <- klb8_thinning70[[i]]$coh.value[2, 2]
  klb8_thinning70[[i]]$turn.value <- metacom::Turnover(klb8_thinned_dat)
  klb8_thinning70[[i]]$turn.z <- klb8_thinning70[[i]]$turn.value[2,2]
  klb8_thinning70[[i]]$bc.value <- metacom::BoundaryClump(klb8_thinned_dat)
  klb8_thinning70[[i]]$bc.i <- klb8_thinning70[[i]]$bc.value[1,2]
  klb8_thinning70[[i]]$bc.p <- klb8_thinning70[[i]]$bc.value[2,2]
  
  #do cooccur
  temp2 <- as.data.frame(t(klb8_thinned_dat))
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb8_thinning70[[i]]$cooc.mat <- cooc.mat
  klb8_thinning70[[i]]$cooc.percent <- temp.res
  klb8_thinning70[[i]]$sp_dat <- klb8_thinned_dat2
}


klb8_thinned.site.number <- lapply(klb8_thinning70, `[[`,3)
klb8_thinned.site.number <- do.call(rbind, klb8_thinned.site.number)
colnames(klb8_thinned.site.number) <- "sampled_sites"

coh.res.all <- lapply(klb8_thinning70, `[[`,5)
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb8_thinning70, `[[`,7)
turn.res.all <- do.call(rbind, turn.res.all)

bc.i.all <- lapply(klb8_thinning70, `[[`,9)
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb8_thinning70, `[[`,10)
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb8_thinning70, `[[`,12)
cooc.res.all <- do.call(rbind, cooc.res.all)
klb8_thinned_summary_res70 <- as.data.frame(cbind(klb8_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb8_thinned_summary_res70) <- c("sampled_sites", "coherence",
                                          "turnover", "bc.i", "bc.p", "cooccur")



# KLB 9  ------------------------------------------------------------------


data.9 <- data.pri[data.pri$KLB == 9,]
rownames(data.9) <- data.9$station
klb9.st <- data.9[1:3]
data.9 <- data.9[,12:66]
data.9 <- data.9[,-c(48, 52)]
data.9 <- modified.remove.zeroes(data.9)

#turn coords into an sf object 
library(sf)
library(tidysdm)

stations <- rownames(data.9)
stations <- as.numeric(stations)
klb9.st <- klb9.st[klb9.st$station %in% stations,]

data_coords <- st_as_sf(klb9.st, coords = c("Longitude", "Latitude"))
st_crs(data_coords) <- 4326

#klb9 thinning 30 ----

klb9_thinning30 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb9_thinning30)[i] <- paste("trial_",i)
  klb9_thinning30[[i]]$klb9_thinned_st <- thin_by_dist(data_coords, 30)
  klb9_thinned_dat <- data.pri[data.pri$station %in% klb9_thinning30[[i]]$klb9_thinned_st$station,]
  rownames(klb9_thinned_dat) <- klb9_thinned_dat$station
  #keep benthic only 
  klb9_thinned_dat <- klb9_thinned_dat[,12:66]
  klb9_thinned_dat <- klb9_thinned_dat[,-c(48, 52)]
  klb9_thinned_dat2 <- modified.remove.zeroes(klb9_thinned_dat)
  klb9_thinned_dat <- remove.sing.zeroes(klb9_thinned_dat)
  klb9_thinning30[[i]]$klb9_thinned_dat <- klb9_thinned_dat
  klb9_thinning30[[i]]$klb9_thinned_sites <- nrow(klb9_thinned_dat)
  #do metacom
  klb9_thinning30[[i]]$coh.value <- metacom::Coherence(klb9_thinned_dat)
  klb9_thinning30[[i]]$coh.z <- klb9_thinning30[[i]]$coh.value[2, 2]
  klb9_thinning30[[i]]$turn.value <- metacom::Turnover(klb9_thinned_dat)
  klb9_thinning30[[i]]$turn.z <- klb9_thinning30[[i]]$turn.value[2,2]
  klb9_thinning30[[i]]$bc.value <- metacom::BoundaryClump(klb9_thinned_dat)
  klb9_thinning30[[i]]$bc.i <- klb9_thinning30[[i]]$bc.value[1,2]
  klb9_thinning30[[i]]$bc.p <- klb9_thinning30[[i]]$bc.value[2,2]
  
  #do cooccur
  temp2 <- as.data.frame(t(klb9_thinned_dat))
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb9_thinning30[[i]]$cooc.mat <- cooc.mat
  klb9_thinning30[[i]]$cooc.percent <- temp.res
  klb9_thinning30[[i]]$sp_dat <- klb9_thinned_dat2
}


klb9_thinned.site.number <- lapply(klb9_thinning30, `[[`,3)
klb9_thinned.site.number <- do.call(rbind, klb9_thinned.site.number)
colnames(klb9_thinned.site.number) <- "sampled_sites"

coh.res.all <- lapply(klb9_thinning30, `[[`,5)
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb9_thinning30, `[[`,7)
turn.res.all <- do.call(rbind, turn.res.all)

bc.i.all <- lapply(klb9_thinning30, `[[`,9)
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb9_thinning30, `[[`,10)
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb9_thinning30, `[[`,12)
cooc.res.all <- do.call(rbind, cooc.res.all)
klb9_thinned_summary_res30 <- as.data.frame(cbind(klb9_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb9_thinned_summary_res30) <- c("sampled_sites", "coherence",
                                          "turnover", "bc.i", "bc.p", "cooccur")

# klb9 thinning 50 --------------------------------------------------------


klb9_thinning50 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb9_thinning50)[i] <- paste("trial_",i)
  klb9_thinning50[[i]]$klb9_thinned_st <- thin_by_dist(data_coords, 50)
  klb9_thinned_dat <- data.pri[data.pri$station %in% klb9_thinning50[[i]]$klb9_thinned_st$station,]
  rownames(klb9_thinned_dat) <- klb9_thinned_dat$station
  #keep benthic only 
  klb9_thinned_dat <- klb9_thinned_dat[,12:66]
  klb9_thinned_dat <- klb9_thinned_dat[,-c(48, 52)]
  klb9_thinned_dat2 <- modified.remove.zeroes(klb9_thinned_dat)
  klb9_thinned_dat <- remove.sing.zeroes(klb9_thinned_dat)
  klb9_thinning50[[i]]$klb9_thinned_dat <- klb9_thinned_dat
  klb9_thinning50[[i]]$klb9_thinned_sites <- nrow(klb9_thinned_dat)
  #do metacom
  klb9_thinning50[[i]]$coh.value <- metacom::Coherence(klb9_thinned_dat)
  klb9_thinning50[[i]]$coh.z <- klb9_thinning50[[i]]$coh.value[2, 2]
  klb9_thinning50[[i]]$turn.value <- metacom::Turnover(klb9_thinned_dat)
  klb9_thinning50[[i]]$turn.z <- klb9_thinning50[[i]]$turn.value[2,2]
  klb9_thinning50[[i]]$bc.value <- metacom::BoundaryClump(klb9_thinned_dat)
  klb9_thinning50[[i]]$bc.i <- klb9_thinning50[[i]]$bc.value[1,2]
  klb9_thinning50[[i]]$bc.p <- klb9_thinning50[[i]]$bc.value[2,2]
  
  #do cooccur
  temp2 <- as.data.frame(t(klb9_thinned_dat))
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb9_thinning50[[i]]$cooc.mat <- cooc.mat
  klb9_thinning50[[i]]$cooc.percent <- temp.res
  klb9_thinning50[[i]]$sp_dat <- klb9_thinned_dat2
}


klb9_thinned.site.number <- lapply(klb9_thinning50, `[[`,3)
klb9_thinned.site.number <- do.call(rbind, klb9_thinned.site.number)
colnames(klb9_thinned.site.number) <- "sampled_sites"

coh.res.all <- lapply(klb9_thinning50, `[[`,5)
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb9_thinning50, `[[`,7)
turn.res.all <- do.call(rbind, turn.res.all)

bc.i.all <- lapply(klb9_thinning50, `[[`,9)
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb9_thinning50, `[[`,10)
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb9_thinning50, `[[`,12)
cooc.res.all <- do.call(rbind, cooc.res.all)
klb9_thinned_summary_res50 <- as.data.frame(cbind(klb9_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb9_thinned_summary_res50) <- c("sampled_sites", "coherence",
                                          "turnover", "bc.i", "bc.p", "cooccur")






# klb9 thinning 70 --------------------------------------------------------


klb9_thinning70 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb9_thinning70)[i] <- paste("trial_",i)
  klb9_thinning70[[i]]$klb9_thinned_st <- thin_by_dist(data_coords, 70)
  klb9_thinned_dat <- data.pri[data.pri$station %in% klb9_thinning70[[i]]$klb9_thinned_st$station,]
  rownames(klb9_thinned_dat) <- klb9_thinned_dat$station
  #keep benthic only 
  klb9_thinned_dat <- klb9_thinned_dat[,12:66]
  klb9_thinned_dat <- klb9_thinned_dat[,-c(48, 52)]
  klb9_thinned_dat2 <- modified.remove.zeroes(klb9_thinned_dat)
  klb9_thinned_dat <- remove.sing.zeroes(klb9_thinned_dat)
  klb9_thinning70[[i]]$klb9_thinned_dat <- klb9_thinned_dat
  klb9_thinning70[[i]]$klb9_thinned_sites <- nrow(klb9_thinned_dat)
  #do metacom
  klb9_thinning70[[i]]$coh.value <- metacom::Coherence(klb9_thinned_dat)
  klb9_thinning70[[i]]$coh.z <- klb9_thinning70[[i]]$coh.value[2, 2]
  klb9_thinning70[[i]]$turn.value <- metacom::Turnover(klb9_thinned_dat)
  klb9_thinning70[[i]]$turn.z <- klb9_thinning70[[i]]$turn.value[2,2]
  klb9_thinning70[[i]]$bc.value <- metacom::BoundaryClump(klb9_thinned_dat)
  klb9_thinning70[[i]]$bc.i <- klb9_thinning70[[i]]$bc.value[1,2]
  klb9_thinning70[[i]]$bc.p <- klb9_thinning70[[i]]$bc.value[2,2]
  
  #do cooccur
  temp2 <- as.data.frame(t(klb9_thinned_dat))
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb9_thinning70[[i]]$cooc.mat <- cooc.mat
  klb9_thinning70[[i]]$cooc.percent <- temp.res
  klb9_thinning70[[i]]$sp_dat <- klb9_thinned_dat2
}


klb9_thinned.site.number <- lapply(klb9_thinning70, `[[`,3)
klb9_thinned.site.number <- do.call(rbind, klb9_thinned.site.number)
colnames(klb9_thinned.site.number) <- "sampled_sites"

coh.res.all <- lapply(klb9_thinning70, `[[`,5)
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb9_thinning70, `[[`,7)
turn.res.all <- do.call(rbind, turn.res.all)

bc.i.all <- lapply(klb9_thinning70, `[[`,9)
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb9_thinning70, `[[`,10)
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb9_thinning70, `[[`,12)
cooc.res.all <- do.call(rbind, cooc.res.all)
klb9_thinned_summary_res70 <- as.data.frame(cbind(klb9_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb9_thinned_summary_res70) <- c("sampled_sites", "coherence",
                                          "turnover", "bc.i", "bc.p", "cooccur")


# save all results --------------------------------------------------------

save(klb7_thinned_summary_res30, klb7_thinned_summary_res50, 
     klb7_thinned_summary_res70, klb7_thinning30, klb7_thinning70,
     klb7_thinning50, 
     klb8_thinned_summary_res30, klb8_thinned_summary_res50, 
     klb8_thinned_summary_res70, klb8_thinning30, klb8_thinning70,
     klb8_thinning50, 
     klb9_thinned_summary_res30, klb9_thinned_summary_res50, 
     klb9_thinned_summary_res70, klb9_thinning30, klb9_thinning70,
     klb9_thinning50,
     file = "results/all_thinning_res_08Dec.RData")


#--- species accum ------- goes to end ---------

windows()
svg(file = "figures/specAccum_thinning.svg", h = 14, w = 10)
par(mfrow = c(3,3))

plot(sp.7, main = "Species accumulation KLB 7, thinning by 30m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", 
     xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb7_thinning30[[i]]$sp_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, 
       col = adjustcolor("cornflowerblue", alpha.f = 0.2), 
       lwd = 1.5)
}

plot(sp.7, main = "Species accumulation KLB 7, thinning by 50m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb7_thinning50[[i]]$sp_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, 
       col = adjustcolor("cornflowerblue", alpha.f = 0.2), 
       lwd = 1.5)
}

plot(sp.7, main = "Species accumulation KLB 7, thinning by 70m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb7_thinning70[[i]]$sp_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE,
       col = adjustcolor("cornflowerblue", alpha.f = 0.2),
       lwd = 1.5)
}


plot(sp.8, main = "Species accumulation KLB 8, thinning by 30m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb8_thinning30[[i]]$sp_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, lwd = 1.5, 
       col = adjustcolor("forestgreen", alpha.f = 0.2))
}

plot(sp.8, main = "Species accumulation KLB 8, thinning by 50m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb8_thinning50[[i]]$sp_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, lwd = 1.5, 
       col = adjustcolor("forestgreen", alpha.f = 0.2))
}

plot(sp.8, main = "Species accumulation KLB 8, thinning by 70m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb8_thinning70[[i]]$sp_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, lwd = 1.5, 
       col = adjustcolor("forestgreen", alpha.f = 0.2))
}


plot(sp.9, main = "Species accumulation KLB 9, thinning by 30m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb9_thinning30[[i]]$sp_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, lwd = 1.5, 
       col = adjustcolor("darkorange", alpha.f = 0.2))
}

plot(sp.9, main = "Species accumulation KLB 9, thinning by 50m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb9_thinning50[[i]]$sp_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, lwd = 1.5, 
       col = adjustcolor("darkorange", alpha.f = 0.2))
}

plot(sp.9, main = "Species accumulation KLB 9, thinning by 70m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb9_thinning70[[i]]$sp_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, lwd = 1.5, 
       col = adjustcolor("darkorange", alpha.f = 0.2))
}
dev.off()


# thinning coocur boxplots ------------------------------------------------

load(file = "results/all_thinning_res_08Dec.RData")
load(file = "results/cooccur.sig.res.RData")
cooccur.sig.res <- cooccur.sig.res[-c(4, 9),]

#klb7
klb7_thin30.cooc <- as.data.frame(klb7_thinned_summary_res30$cooccur)
klb7_thin30.cooc$klb <- rep("klb7_thin30", 10)
klb7_thin30.cooc$treatment <- rep("trials", 10)
colnames(klb7_thin30.cooc) <- c("sig.res", "klb", "treatment")

klb7_thin50.cooc <- as.data.frame(klb7_thinned_summary_res50$cooccur)
klb7_thin50.cooc$klb <- rep("klb7_thin50", 10)
klb7_thin50.cooc$treatment <- rep("trials", 10)
colnames(klb7_thin50.cooc) <- c("sig.res", "klb", "treatment")

klb7_thin70.cooc <- as.data.frame(klb7_thinned_summary_res70$cooccur)
klb7_thin70.cooc$klb <- rep("klb7_thin70", 10)
klb7_thin70.cooc$treatment <- rep("trials", 10)
colnames(klb7_thin70.cooc) <- c("sig.res", "klb", "treatment")

#klb8
klb8_thin30.cooc <- as.data.frame(klb8_thinned_summary_res30$cooccur)
klb8_thin30.cooc$klb <- rep("klb8_thin30", 10)
klb8_thin30.cooc$treatment <- rep("trials", 10)
colnames(klb8_thin30.cooc) <- c("sig.res", "klb", "treatment")

klb8_thin50.cooc <- as.data.frame(klb8_thinned_summary_res50$cooccur)
klb8_thin50.cooc$klb <- rep("klb8_thin50", 10)
klb8_thin50.cooc$treatment <- rep("trials", 10)
colnames(klb8_thin50.cooc) <- c("sig.res", "klb", "treatment")

klb8_thin70.cooc <- as.data.frame(klb8_thinned_summary_res70$cooccur)
klb8_thin70.cooc$klb <- rep("klb8_thin70", 10)
klb8_thin70.cooc$treatment <- rep("trials", 10)
colnames(klb8_thin70.cooc) <- c("sig.res", "klb", "treatment")



#klb9
thin30.cooc <- as.data.frame(klb9_thinned_summary_res30$cooccur)
thin30.cooc$klb <- rep("klb9_thin30", 10)
thin30.cooc$treatment <- rep("trials", 10)
colnames(thin30.cooc) <- c("sig.res", "klb", "treatment")

thin50.cooc <- as.data.frame(klb9_thinned_summary_res50$cooccur)
thin50.cooc$klb <- rep("klb9_thin50", 10)
thin50.cooc$treatment <- rep("trials", 10)
colnames(thin50.cooc) <- c("sig.res", "klb", "treatment")

thin70.cooc <- as.data.frame(klb9_thinned_summary_res70$cooccur)
thin70.cooc$klb <- rep("klb9_thin70", 10)
thin70.cooc$treatment <- rep("trials", 10)
colnames(thin70.cooc) <- c("sig.res", "klb", "treatment")


cooccur.sig.res <- rbind(cooccur.sig.res[1:5,], klb7_thin30.cooc, klb7_thin50.cooc,
                         klb7_thin70.cooc, 
                         cooccur.sig.res[6:10,],
                         klb8_thin30.cooc, klb8_thin50.cooc,
                         klb8_thin70.cooc, 
                         cooccur.sig.res[11:13,],
                         thin30.cooc, thin50.cooc,
                         thin70.cooc)
rownames(cooccur.sig.res) <- NULL
cooccur.sig.res$klb <- 
  factor(cooccur.sig.res$klb, 
      levels = c("klb7", "klb7_thin30", 
        "klb7_thin50", "klb7_thin70", 
          "klb8", "klb8_thin30", 
            "klb8_thin50", "klb8_thin70",
                  "klb9", "klb9_thin30", 
                      "klb9_thin50", "klb9_thin70"
      ))

#plot the data

library(ggplot2)

coocc.box <- ggplot(cooccur.sig.res, aes(x = klb, y = sig.res, 
              fill = klb)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.6)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.02),
             aes(shape = treatment, size = 0.001)) + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=
            c("darkblue", "royalblue", "cornflowerblue",  
                "lightblue", 
              "darkgreen", "forestgreen", "olivedrab","yellowgreen",
            "darkorange", "orange", "salmon", "pink")) +
  ggtitle("Cooccur results") + 
  #coord_flip() +
  theme_bw()
coocc.box
ggsave(coocc.box, file = "figures/cooccur.thinned.svg", 
       h = 10, w = 10, scale = 1)
dev.off()

#draw the metacom plots -----------------

#import klb 7 and 8 and 9 metacom results ----

load(file = "results/klb7.metacom.results.RData")
load(file = "results/klb8.metacom.results.RData")
load(file = "results/klb9.metacom.results.RData")

#plot metacom Res 
coh.res.7 <- (klb7.metacom.res.df[2, c(2, 8, 14, 20, 26)])
turn.res.7 <- (klb7.metacom.res.df[2, c(4, 10, 16, 22, 28)])
bc.res.7 <- (klb7.metacom.res.df[1, c(6, 12, 18, 24, 30)])

coh.res.8 <- (klb8.metacom.res.df[2, c(2, 8, 14, 20, 26)])
turn.res.8 <- (klb8.metacom.res.df[2, c(4, 10, 16, 22, 28)])
bc.res.8 <- (klb8.metacom.res.df[1, c(6, 12, 18, 24, 30)])

coh.res.9 <- (klb9.metacom.res.df[2, c(2, 8, 14)])
turn.res.9 <- (klb9.metacom.res.df[2, c(4, 10, 16)])
bc.res.9 <- (klb9.metacom.res.df[1, c(6, 12, 18)])

windows()
svg(file = "figures/S5-metacom_thinning.svg", h = 7, w = 7)
par(mfrow = c(1,1))
plot(xlim = c(-6, 6), ylim = c(-10, 0), 
     xlab = "Turnover Z-score", ylab = "Coherence Z-score",
     x = -3:2, y = -3:2, col = "white", 
     main = "All Metacommunity Results")

abline(v = 0, lwd = 1)
box()
text(-4, 0, "Nested", cex = 1.2)
text(2, -0, "Non-Nested", cex = 1.2, col = "black")

#klb9
points(turn.res.9$raw.turnover_res, coh.res.9$raw.coherence_res, 
       pch = 22, cex = bc.res.9$raw.boundary_clumping_res, 
       bg = "darkorange", lwd = 1, col = "black")
#expand
points(turn.res.9$treat1.turnover_res, coh.res.9$treat1.coherence_res, 
       pch = 22, cex = bc.res.9$treat1.boundary_clumping_res, 
       bg = "darkorange", lwd = 1, col = "black")

#contract
points(turn.res.9$treat2.turnover_res, coh.res.9$treat2.coherence_res, 
       pch = 22, cex = bc.res.9$treat2.boundary_clumping_res,
       bg = "darkorange", lwd = 1, col = "black")

#klb9 thinned points 

points(klb9_thinned_summary_res30$turnover, klb9_thinned_summary_res30$coherence, 
       pch = 16, cex = klb9_thinned_summary_res30$bc.i, 
       col = adjustcolor("darkorange", alpha.f = 0.8) )

points(klb9_thinned_summary_res50$turnover, klb9_thinned_summary_res50$coherence, 
       pch = 16, cex = klb9_thinned_summary_res50$bc.i, 
       col = adjustcolor("salmon", alpha.f = 0.5) )

points(klb9_thinned_summary_res70$turnover, klb9_thinned_summary_res70$coherence, 
       pch = 16, cex = klb9_thinned_summary_res50$bc.i, 
       col = adjustcolor("pink", alpha.f = 0.5) )



#klb8
points(turn.res.8$raw.turnover_res, coh.res.8$raw.coherence_res, 
       pch = 22, cex = bc.res.8$raw.boundary_clumping_res, 
       bg = "forestgreen", lwd = 1, col = "black")


#left
points(turn.res.8$treat1.turnover_res, coh.res.8$treat1.coherence_res, 
       pch = 22, cex = bc.res.8$treat1.boundary_clumping_res,
       bg = "forestgreen", lwd = 1, col = "black")


#right
points(turn.res.8$treat2.turnover_res, coh.res.8$treat2.coherence_res, 
       pch = 22, cex = bc.res.8$treat2.boundary_clumping_res,
       bg = "forestgreen", lwd = 1, col = "black")


#expand
points(turn.res.8$treat3.turnover_res, coh.res.8$treat3.coherence_res, 
       pch = 22, cex = bc.res.8$treat3.boundary_clumping_res, 
       bg = "forestgreen", lwd = 1, col = "black")

#contract
points(turn.res.8$treat4.turnover_res, coh.res.8$treat4.coherence_res, 
       pch = 22, cex = bc.res.8$treat4.boundary_clumping_res,
       bg = "forestgreen", lwd = 1)

#thinned points
points(klb8_thinned_summary_res30$turnover, klb8_thinned_summary_res30$coherence, 
       pch = 16, cex = klb8_thinned_summary_res30$bc.i, 
       col = adjustcolor("forestgreen", alpha.f = 0.8) )

points(klb8_thinned_summary_res50$turnover, klb8_thinned_summary_res50$coherence, 
       pch = 16, cex = klb8_thinned_summary_res50$bc.i, 
       col = adjustcolor("olivedrab", alpha.f = 0.5) )

points(klb8_thinned_summary_res70$turnover, klb8_thinned_summary_res70$coherence, 
       pch = 16, cex = klb8_thinned_summary_res50$bc.i, 
       col = adjustcolor("yellowgreen", alpha.f = 0.5) )


#klb 7
points(turn.res.7$raw.turnover_res, coh.res.7$raw.coherence_res, 
       pch = 22, cex = bc.res.7$raw.boundary_clumping_res, 
       bg = "darkblue", lwd = 1, col = "black")

#left
points(turn.res.7$treat1.turnover_res, coh.res.7$treat1.coherence_res, 
       pch = 22, cex = bc.res.7$treat1.boundary_clumping_res,
       bg = "darkblue", lwd = 1, col = "black")

#right
points(turn.res.7$treat2.turnover_res, coh.res.7$treat2.coherence_res, 
       pch = 22, cex = bc.res.7$treat2.boundary_clumping_res,
       bg = "darkblue", lwd = 1, col = "black")

#expand
points(turn.res.7$treat3.turnover_res, coh.res.7$treat3.coherence_res, 
       pch = 22, cex = bc.res.7$treat3.boundary_clumping_res, 
       bg = "darkblue", lwd = 1, col = "black")
#contract
points(turn.res.7$treat4.turnover_res, coh.res.7$treat4.coherence_res, 
       pch = 22, cex = bc.res.7$treat4.boundary_clumping_res,
       bg = "darkblue", lwd = 1)

points(klb7_thinned_summary_res30$turnover, klb7_thinned_summary_res30$coherence, 
       pch = 16, cex = klb7_thinned_summary_res30$bc.i, 
       col = adjustcolor("royalblue", alpha.f = 0.8) )

points(klb7_thinned_summary_res50$turnover, klb7_thinned_summary_res50$coherence, 
       pch = 16, cex = klb7_thinned_summary_res50$bc.i, 
       col = adjustcolor("cornflowerblue", alpha.f = 0.5) )

points(klb7_thinned_summary_res70$turnover, klb7_thinned_summary_res70$coherence, 
       pch = 16, cex = klb7_thinned_summary_res50$bc.i, 
       col = adjustcolor("lightblue", alpha.f = 0.5) )
dev.off()









