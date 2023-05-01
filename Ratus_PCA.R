################################################################################
require(raster);require(data.table);require(readxl);require(ff);require(parallel);require(purrr)
library(sf);library(spatialRF);require(geodist);library(rgdal);require(tmap);require(FactoMineR)
require(factoextra);require(corrplot);require(tmap)
################################################################################
rasterOptions(tmpdir="D:/DRYAD/ncfd4/raster_tmpdir")
tempdir="D:/DRYAD/ncfd4/raster_tmpdir"
d_dir <- "D:/DRYAD/ncfd4/FINAL"
r_file <- list.files(d_dir,pattern = ".tif")
out_dir <- "D:/DRYAD/RATTUS"
#loading world shapefile
world <- rnaturalearth::ne_countries(
  scale = "medium", 
  returnclass = "sf"
)

################################################################################
x_data <- readxl::read_xlsx(paste0(out_dir,"/","Matris_R_rattus.xlsx"))
x_data$Longitud <- as.numeric(x_data$Longitud)
x_data$Latitud <- as.numeric(x_data$Latitud)
x_data <- x_data[which(!is.na(x_data$Longitud)),]
x_data <- x_data[which(x_data$REMOVE=="F"),]
x_data$REMOVE <- NULL
x_data$species <- "Rattus rattus"
################################################################################
x_data1 <- readxl::read_xlsx(paste0(out_dir,"/","Matris_R_norvegicus.xlsx"))
x_data1$Longitud <- as.numeric(x_data1$Longitud)
x_data1$Latitud <- as.numeric(x_data1$Latitud)
x_data1 <- x_data1[which(!is.na(x_data1$Longitud)),]
x_data1 <- x_data1[which(x_data1$REMOVE=="F"),]
x_data1$REMOVE <- NULL
x_data1$species <- "Rattus norvegicus"

x_data <-as.data.frame(rbind(x_data,x_data1))
x_data$Ho <- as.numeric(x_data$Ho)
x_data$He <- as.numeric(x_data$He)
x_data$FST <- as.numeric(x_data$FST)
x_data$FIS <- as.numeric(x_data$FIS)
x_data$`Riquesa Alellica` <- as.numeric(x_data$`Riquesa Alellica`)

par(mfrow=c(1,3))

boxplot(x_data$Ho~x_data$species,ylab="Observed heterozygosity",xlab="",main="")
boxplot(x_data$He~x_data$species,ylab="Expected heterozygosity",xlab="",main="")
boxplot(x_data$FST~x_data$species,ylab="Fixation index",xlab="",main="")

################################################################################
#read raster files
layers2 <- list.files(paste0(d_dir,"/RESAMPLED"),pattern = ".tif")
layers2 <- layers2[c(1,2,3,4,5,6)] #ayers2 <- layers2[c(1,2,3,4,5,6)] #no NDVI
layers2 <-  lapply(1:length(layers2),function(i){
  x <- raster(paste0(d_dir,"/RESAMPLED/",layers2[[i]]))
  return(x)
})
layers <- stack(layers2)
################################################################################
#creating a mask
sample_p <- layers2[3][[1]]
values(sample_p)[!is.na(values(sample_p))] = 1
# plot(sample_p)
################################################################################
#fitting to poverty raster
layers3 <- lapply(1:length(layers2),function(i){
  x <- layers2[[i]]*sample_p
  return(x)
})
layers3 <- stack(layers3)
names(layers3) <- names(layers)
#plot(layers3)

#layers3a <- scale(layers3)
#https://gis.stackexchange.com/questions/158103/save-raster-to-file-when-using-stack-in-r
#writeRaster(layers3a, filename=paste0(out_dir,"/","multilayer.tif"), options="INTERLEAVE=BAND", overwrite=TRUE)
# ss <- stack(paste0(out_dir,"/","multilayer.tif"))
# names(ss) <-  c("GDP_per_capita_PPP_2015","gpw_v4_population_density_2020","povmap.grdi.v1",
#                 "travel_time_to_cities_12_MOD","travel_time_to_ports_5_MOD","wc2.1_2.5m_bio_1",
#                 "wc2.1_2.5m_bio_12","wc2.1_2.5m_bio_15","wc2.1_2.5m_bio_2","wc2.1_2.5m_bio_3")

################################################################################
#extracting raster values
#plot(layers)
xx_ext <- raster::extract(layers, sp::SpatialPoints(cbind(x_data$Longitud,x_data$Latitud)), sp = T)
#xx_ext <- raster::extract(layers3a, sp::SpatialPoints(cbind(x_data$LONGITUD,x_data$LATITUD)), sp = T)

#View(as.data.frame(xx_ext))
x_data <- cbind(x_data,xx_ext)
#colnames(x_data)
x_data <- x_data[,-c((ncol(x_data)-2):ncol(x_data))]


x_data <- x_data[which(!is.na(x_data$povmap.grdi.v1)),]
#converting charcters to numbers
x_data$He <- as.numeric(x_data$He)
x_data$Ho <- as.numeric(x_data$Ho)
x_data$FIS <- as.numeric(x_data$FIS)
x_data$FST <- as.numeric(x_data$FST)
x_data$`Riquesa Alellica` <- as.numeric(x_data$`Riquesa Alellica`)


length(na.omit(x_data$He))
length(na.omit(x_data$Ho))
length(na.omit(x_data$FST))
#length(na.omit(x_data$FIS))
#length(na.omit(x_data$FIS))
VAR <-c("Ho")#,"FST")
################################################################################
##summaries

tapply(x_data$species, x_data$species,length)

#Ho mean and sd

tapply(x_data[which(!is.na(x_data$Ho)),]$species, x_data[which(!is.na(x_data$Ho)),]$species,length)
tapply(x_data[which(!is.na(x_data$Ho)),]$Ho, x_data[which(!is.na(x_data$Ho)),]$species,function(x){mean(x,na.rm=T)})
tapply(x_data[which(!is.na(x_data$Ho)),]$Ho, x_data[which(!is.na(x_data$Ho)),]$species,function(x){sd(x,na.rm=T)})

#Ho mean and sd
tapply(x_data[which(!is.na(x_data$He)),]$species, x_data[which(!is.na(x_data$He)),]$species,length)
tapply(x_data[which(!is.na(x_data$He)),]$He, x_data[which(!is.na(x_data$He)),]$species,function(x){mean(x,na.rm=T)})
tapply(x_data[which(!is.na(x_data$He)),]$He, x_data[which(!is.na(x_data$He)),]$species,function(x){sd(x,na.rm=T)})


#Fst mean and sd
tapply(x_data[which(!is.na(x_data$FST)),]$species, x_data[which(!is.na(x_data$FST)),]$species,length)
tapply(x_data[which(!is.na(x_data$FST)),]$FST, x_data[which(!is.na(x_data$FST)),]$species,function(x){mean(x,na.rm=T)})
tapply(x_data[which(!is.na(x_data$FST)),]$FST, x_data[which(!is.na(x_data$FST)),]$species,function(x){sd(x,na.rm=T)})

#Fis mean and sd
tapply(x_data[which(!is.na(x_data$FIS)),]$species, x_data[which(!is.na(x_data$FIS)),]$species,length)
tapply(x_data[which(!is.na(x_data$FIS)),]$FIS, x_data[which(!is.na(x_data$FIS)),]$species,function(x){mean(x,na.rm=T)})
tapply(x_data[which(!is.na(x_data$FIS)),]$FIS, x_data[which(!is.na(x_data$FIS)),]$species,function(x){sd(x,na.rm=T)})


#Fis mean and sd
tapply(x_data[which(!is.na(x_data$`Riquesa Alellica`)),]$species, x_data[which(!is.na(x_data$`Riquesa Alellica`)),]$species,length)
tapply(x_data[which(!is.na(x_data$`Riquesa Alellica`)),]$`Riquesa Alellica`, x_data[which(!is.na(x_data$`Riquesa Alellica`)),]$species,function(x){mean(x,na.rm=T)})
tapply(x_data[which(!is.na(x_data$`Riquesa Alellica`)),]$`Riquesa Alellica`, x_data[which(!is.na(x_data$`Riquesa Alellica`)),]$species,function(x){sd(x,na.rm=T)})


#load(paste0("D:/DRYAD/RATTUS/",VAR,".RData"))
#"He" #"FST" #"HETEROCIGOCIDAD HE", "HETEROCIGOCIDAD HO"
################################################################################
#subsetting to Ho

VAR <- c("Ho")

X_FST <- x_data[which(!is.na(x_data[,VAR])),]#`HETEROCIGOCIDAD HO`)),]#`HETEROCIGOCIDAD HO`)),]
#plot(X_FST$LONGITUD,X_FST$LATITUD,pch=15)

X_FST <- X_FST[,c("Longitud","Latitud","species",VAR,#"HETEROCIGOCIDAD HO",#"HETEROCIGOCIDAD HO",#"FST", 
                  "GDP_per_capita_PPP_2015","gpw_v4_population_density_2020",
                  "povmap.grdi.v1","travel_time_to_cities_12_MOD",#,
                  "travel_time_to_ports_5_MOD")]

colnames(X_FST) <- c("Longitud","Latitud","species","observed heterozygosity",
                     "GDP per capita (2015)","population density (2020)",
                     "poverry gridded", "Travel time to cities",
                     "Travel time to ports")#,"NDVI (median 2015-2019)")
#"wc2.1_2.5m_bio_1",
#  "wc2.1_2.5m_bio_12","wc2.1_2.5m_bio_15","wc2.1_2.5m_bio_2"

#colnames(X_FST)[3] <- VAR #FST
X_FST <- X_FST[which(complete.cases(X_FST)),]

tapply(X_FST$species,X_FST$species,length)
res.pca = PCA(X_FST[,-c(1,2,3)],scale.unit = T)
x_H <- HCPC(res.pca, nb.clust = -1, min = 3, max = NULL, graph = TRUE)
# fviz_dend(x_H, 
#           cex = 0.7,                     # Label size
#           palette = "jco",               # Color palette see ?ggpubr::ggpar
#           rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
#           rect_border = "jco",           # Rectangle color
#           labels_track_height = 0.8      # Augment the room for labels
# )


fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
fviz_pca_var(res.pca, col.var = "black")
#library("corrplot")
#corrplot(var$cos2, is.corr=FALSE)


tapply(X_FST$species,X_FST$species,length)
fviz_pca_biplot(res.pca, geom.ind = "point",
               # col.ind = x_H$data.clust$clust, # color by groups
               col.ind = X_FST$species, # color by groups
               col.var = "black",
                #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                addEllipses = TRUE, ellipse.type = "convex",
                legend.title = ""
)

#X_FST <- X_FST[-c(18,41),] #FST

#X_FST <- X_FST[-c(68),] #He
apply(X_FST[,-c(1,2,3)], 2, var) == 0


x_sp_RR <- X_FST[which(X_FST$species=="Rattus rattus"),]
x_sp_RN <- X_FST[which(X_FST$species=="Rattus norvegicus"),]




pdf(paste0("D:/DRYAD/RATTUS/GRAPHICS","/","CORRPLOT.pdf"),height = 8, width = 30)#15


par(mfrow=c(1,2))

corrplot::corrplot(cor(x_sp_RN[,-c(1,2,3)]),  method = 'square', type = 'lower', insig='blank',
                   addCoef.col ='black', number.cex = 0.8, diag=FALSE)#,

corrplot::corrplot(cor(x_sp_RR[,-c(1,2,3)]), method = 'square', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, diag=FALSE)#,

        # main = substitute(paste(italic('Rattus rattus'))))
dev.off()


################################################################################
#MAPS

data("World")
coords <- X_FST[which(!is.na(X_FST$Longitud)),]
coords$Longitud <- as.numeric(as.character(coords$Longitud))
coords$Latitud <- as.numeric(as.character(coords$Latitud))
coordinates(coords)  <- ~Longitud+Latitud
xmin<-round(extent(coords)[1]) ;if(xmin< -180){xmin=-180}
xmax<-round(extent(coords)[2])
ymin<-round(extent(coords)[3])
ymax<-round(extent(coords)[4]) ;if(ymax>90){ymax=90}
sp_NA <-  extent(c(xmin,xmax,ymin,ymax))
coords <- st_as_sf(coords)
coords

sp_NA <- st_bbox(c(xmin=xmin-0.5, xmax=xmax+0.5, ymax=ymax+0.5, ymin=ymin-0.5), crs = st_crs(World))



#st_crs(coords) <- st_crs(World)

tmap_options(check.and.fix = TRUE)
tm_map_points <- 
  tm_shape(World)+#,bbox = sp_NA) +  
  tm_fill(col="gray91",legend.show = F) +
    #Rattus rattus
  tm_shape(coords[which(coords$species=="Rattus rattus"),])+#,bbox = sp_NA)+ 
  tm_dots(legend.show = F,col = "blue",size = 4) + #purple
    #Rattus novergicus
    tm_shape(coords[which(coords$species=="Rattus norvegicus"),])+#,bbox = sp_NA)+ 
    tm_dots(legend.show = F,col = "red",size = 4) + #purple
  tm_shape(World,bbox = sp_NA) +  tm_borders("black") +
  
  tm_facets(nrow = 1, sync = TRUE)+
  tm_layout(inner.margins=0,
            legend.text.size=10,
            legend.title.size=10,
            legend.position = c("left","bottom"),
            legend.bg.color = "white", legend.bg.alpha=.2)#;gc()

tmap_save(filename=paste0(out_dir,"/","POINTS_PCA",".pdf"),tm=tm_map_points,dpi=900,width =100,height=100,units = "cm");gc()
################################################################################
###TOTAL POINTS

coords <- x_data
coords$Longitud <- as.numeric(as.character(coords$Longitud))
coords$Latitud <- as.numeric(as.character(coords$Latitud))
coordinates(coords)  <- ~Longitud+Latitud
xmin<-round(extent(coords)[1]) ;if(xmin< -180){xmin=-180}
xmax<-round(extent(coords)[2])
ymin<-round(extent(coords)[3])
ymax<-round(extent(coords)[4]) ;if(ymax>90){ymax=90}
sp_NA <-  extent(c(xmin,xmax,ymin,ymax))
coords <- st_as_sf(coords)

sp_NA <- st_bbox(c(xmin=xmin-0.5, xmax=xmax+0.5, ymax=ymax+0.5, ymin=ymin-0.5), crs = st_crs(World))



#st_crs(coords) <- st_crs(World)

tmap_options(check.and.fix = TRUE)
tm_map_points <- 
  tm_shape(World)+#,bbox = sp_NA) +  
  tm_fill(col="gray91",legend.show = F) +
  #Rattus rattus
  tm_shape(coords[which(coords$species=="Rattus rattus"),])+#,bbox = sp_NA)+ 
  tm_dots(legend.show = F,col = "blue",size = 4) + #purple
  #Rattus novergicus
  tm_shape(coords[which(coords$species=="Rattus norvegicus"),])+#,bbox = sp_NA)+ 
  tm_dots(legend.show = F,col = "red",size = 4) + #purple
  tm_shape(World,bbox = sp_NA) +  tm_borders("black") +
  
  tm_facets(nrow = 1, sync = TRUE)+
  tm_layout(inner.margins=0,
            legend.text.size=10,
            legend.title.size=10,
            legend.position = c("left","bottom"),
            legend.bg.color = "white", legend.bg.alpha=.2)#;gc()

tmap_save(filename=paste0(out_dir,"/","POINTS_TOTAL",".pdf"),tm=tm_map_points,dpi=900,width =100,height=100,units = "cm");gc()
################################################################################
# require(vegan)
# 
# 
# data_unique <- unique(X_FST[ , c("Longitud", "Latitud")])
# sp_unique <- unique(X_FST$species)
# 
# pres_matrix <- list()
# pres_env <- list()
# for(i in 1:nrow(data_unique)){
# pres_env[[i]] <- X_FST[which(X_FST$Longitud==data_unique$Longitud[[i]] & 
#               X_FST$Latitud==data_unique$Latitud[[i]]),]
# 
# x_j <- data.frame(matrix(nrow = 1,ncol = 2))
# for(j in 1:length(sp_unique)){
#   if(any(pres_env[[i]] %in% sp_unique[j])){
#     x_j[,j] <- 1
#   } else{
#     x_j[,j] <- 0
#   }
# }
# pres_matrix[[i]] <- x_j;rm(x_j)
# }
# 
# pres_matrix <- do.call(rbind,pres_matrix)
# pres_env <- do.call(rbind,pres_env)
# pres_env <- pres_env[,-c(1,2,3)]
# colnames(pres_env)[1] <- "Ho"
# # Model the effect of all environmental variables on fish
# # community composition
# spe.rda <- rda(pres_matrix ~., data = pres_env)
# summary(spe.rda)
# 
# fwd.sel <- ordiR2step(rda(pres_matrix ~ 1, data = pres_env), # lower model limit (simple!)
#                       scope = formula(spe.rda), # upper model limit (the "full" model)
#                       direction = "forward",
#                       R2scope = TRUE, # can't surpass the "full" model's R2
#                       pstep = 1000,
#                       trace = FALSE) # change to TRUE to see the selection process!
# 
# fwd.sel$call
# 
# #https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
# # Write our new model
# spe.rda.signif <- rda(formula = pres_matrix ~ `GDP per capita (2015)` + `Travel time to ports` + 
#                         Ho, data = pres_env)
# # check the adjusted R2 (corrected for the number of
# # explanatory variables)
# RsquareAdj(spe.rda.signif)
# anova.cca(spe.rda.signif, step = 1000)
# anova.cca(spe.rda.signif, step = 1000, by = "term")
# anova.cca(spe.rda.signif, step = 1000, by = "axis")
# 
# # Type 1 scaling
# ordiplot(spe.rda.signif, scaling = 1, type = "text")
# # Type 2 scaling
# ordiplot(spe.rda.signif, scaling = 2, type = "text")
# ################################################################################
# #library(pls)
# library(pls)
# library(mdatools)

##https://www.statology.org/partial-least-squares-in-r/
##make this example reproducible

##https://mdatools.com/docs/pls--models-and-results.html

# set.seed(1)
# X_FST_RN <- X_FST[which(X_FST$species=="Rattus norvegicus"),]
# X_FST_RR <- X_FST[which(X_FST$species=="Rattus rattus"),]
# #fit PCR model
# 
# Xc <- X_FST_RN$`observed heterozygosity`
# Yc <- X_FST_RN[,c(               "GDP per capita (2015)" ,
#                                 "population density (2020)" , 
#                                  "poverry gridded" ,
#                                  "Travel time to cities" , 
#                                  "Travel time to ports")]
# model <- mdatools::pls(as.data.frame(Xc),Yc,scale=TRUE, cv = list("ven", 3))
# print(model)
# res = predict(model, as.data.frame(Xc), Yc)
# print(res$rmse)
# summary(model)
# plot(model)
