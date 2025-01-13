#this script creates the data tables for further use
# all data classed by station number (270 stations)

library(dplyr)

# Filter by station -------------------------------------------------------

df <- read.csv("data/PRI_Cretaceous_cleaned2_edited_oct24.csv")
pri.cr.data <- df
pri.cr.data$setfam <- paste(pri.cr.data$Set, pri.cr.data$Family, 
                            sep = "_")
pri.cr.data$setfamgen <- paste(pri.cr.data$Set, pri.cr.data$Family, 
                               pri.cr.data$Genus, sep = "_")
pri.cr.data <- pri.cr.data[,c("PRI_Station", "KLB", "setfam", "setfamgen", "Set", "Family", 
                              "Genus", "Species", "Abundance")]

#remove all NAs 
pri.cr.data[pri.cr.data == ""] <- NA
pri.cr.data[pri.cr.data == "unknown"] <- NA
pri.cr.data <- na.omit(pri.cr.data)

klb <- as.numeric(names(table(pri.cr.data$KLB)))

#change klb values in df to be numeric

pri.cr.data$KLB <- as.numeric(pri.cr.data$KLB)


#sort by PRI station numbers 
pri.cr.data <- pri.cr.data[order(pri.cr.data$PRI_Station),]

fam.list <- (unique(pri.cr.data$setfam))

#loop below around all the stations 

list.stations <-list()
list.resgroups<-list()
temp <- list()
temp2 <- list()
KLB <- klb 
station <- as.numeric(names(table(pri.cr.data$PRI_Station)))

for(j in 1:length(station)){
  
  #extract the rows with specified station 
  list.stations[[j]]<-as.data.frame(filter(pri.cr.data,PRI_Station == station[j]))
  
  #gets the setfam IN THE SPECIFIED STATION ONLY 
  fam.names<-names(table(list.stations[[j]][,3])) 
  #new res
  
  list.resgroups[[j]]<-matrix(NA, length(fam.names), 6)
  
  for(i in 1:length(fam.names)) #only ones in station
  {
    line1<-as.data.frame(filter(list.stations[[j]],setfam == fam.names[i]))
    
    #find index number of the family being summed 
    list.resgroups[[j]][i,1] <- station[j]
    list.resgroups[[j]][i,2:5]<-as.matrix(line1[1,c(3, 5, 6, 7)]) 
    #sum the abundance info 
    
    list.resgroups[[j]][i,6]<-sum(as.numeric(line1[,9]))
    list.resgroups[[j]]
    
  }
  
  
}

results <- list()
for (i in 1:length(list.resgroups)){
  temp <- list.resgroups[[i]] #select station 
  y <- nrow(temp)
  station.temp <- matrix(NA, nrow = 1, ncol = y)
  station.temp[1,] <- temp[,6] #abundance
  colnames(station.temp) <- temp[,2] #set fam 
  names(station.temp) <- colnames(station.temp)
  results[[i]] <- station.temp
}

#create site-by-taxa matrix 
abu.station <- bind_rows(results) 
abu.station [is.na(abu.station )] <- "0"
abu.station  <- as.data.frame(sapply(abu.station , as.numeric))

new_order <- sort(colnames(abu.station)) #make taxa alphabetical
abu.station <- abu.station[, new_order]

abu.station <- cbind(station, abu.station)

#bind with lat-long

klb.st <- read.csv(file = "data/station_strat.csv")
coord.st <- df[,c("PRI_Station", "Latitude", "Longitude")]
colnames(coord.st) <- c("station", "Latitude", "Longitude")
coord.st <- unique(coord.st)
klb.st <- merge(klb.st, coord.st, by= "station")

df <- merge(klb.st, abu.station, by = "station") 
#this assigned KLB number to each station
df <- df[order(df$klb),]


# find cold seeps stations ------------------------------------------------

#find Lucina, Solemya, Thyasira or Conchocele 

#Bivalve_Thyasiridae [36]
#Bivalve_Solemyidae [34]
#Bivalve_Lucinidae [22]

seeps <- df[,c(22, 34, 36)]
rownames(seeps) <- df$station
seeps <- seeps[, colSums(seeps != 0) > 0]
seeps <- seeps[rowSums(seeps != 0) > 0,]

seep.station <- rownames(seeps)
seep.station <- as.integer(seep.station)
#find stations to delete
del <- df[df$station %in% seep.station,]
del.ind <- as.numeric(rownames(del))
df <- df[-del.ind,]

#write.csv(df, row.names = FALSE, 
#          file = "data/PRI_Cretaceous_site_by_taxa_Oct30.csv")




