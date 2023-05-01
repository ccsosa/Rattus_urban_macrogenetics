################################################################################
require(raster);require(data.table);require(readxl);require(ff);require(parallel);require(purrr)
library(sf);
#library(spatialRF);require(geodist);
library(rgdal);require(tmap);require(FactoMineR)
require(factoextra);library(ranger);library(caret)
library(doParallel);require(tmaptools);require(dismo)
################################################################################
rasterOptions(tmpdir="D:/DRYAD/ncfd4/raster_tmpdir")
tempdir="D:/DRYAD/ncfd4/raster_tmpdir"
d_dir <- "D:/DRYAD/ncfd4/FINAL"
r_file <- list.files(d_dir,pattern = ".tif")
out_dir <- "D:/DRYAD/RATTUS"

################################################################################
#read data
#https://blasbenito.github.io/spatialRF/
x_data <- readxl::read_xlsx(paste0(out_dir,"/","Matris_R_rattus.xlsx"))
x_data$Longitud <- as.numeric(x_data$Longitud)
x_data$Latitud <- as.numeric(x_data$Latitud)
x_data <- x_data[which(!is.na(x_data$Longitud)),]
x_data <- x_data[which(x_data$REMOVE=="F"),]
x_data$REMOVE <- NULL
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
###############################################################################
#extracting raster values
xx_ext <- raster::extract(layers, sp::SpatialPoints(cbind(x_data$Longitud,x_data$Latitud)), sp = T)
x_data <- cbind(x_data,xx_ext)
x_data <- x_data[,-c((ncol(x_data)-2):ncol(x_data))]
################################################################################
x_data <- x_data[which(!is.na(x_data$povmap.grdi.v1)),]
#converting charcters to numbers
x_data$He <- as.numeric(x_data$He)
x_data$Ho <- as.numeric(x_data$Ho)
x_data$FIS <- as.numeric(x_data$FIS)
x_data$FST <- as.numeric(x_data$FST)
x_data$`Riquesa Alellica` <- as.numeric(x_data$`Riquesa Alellica`)

length(na.omit(x_data$He))
length(na.omit(x_data$Ho))
################################################################################
VAR <-"Ho"
################################################################################
#subsetting to He
X_FST <- x_data[which(!is.na(x_data[,VAR])),]#`HETEROCIGOCIDAD HO`)),]#`HETEROCIGOCIDAD HO`)),]

X_FST <- X_FST[,c("Longitud","Latitud",VAR,#"HETEROCIGOCIDAD HO",#"HETEROCIGOCIDAD HO",#"FST", 
                  "GDP_per_capita_PPP_2015","gpw_v4_population_density_2020",
                  "povmap.grdi.v1","travel_time_to_cities_12_MOD",#,
                  "travel_time_to_ports_5_MOD")]

colnames(X_FST)[3] <- VAR #FST
X_FST <- X_FST[which(complete.cases(X_FST)),]
################################################################################
mask1 <- layers3[[1]]
mask1[which(!is.na(mask1[]))] <- 0
plot(mask1)
rr_points <- raster::rasterize(sp::SpatialPoints(cbind(X_FST$Longitud,X_FST$Latitud)),mask1,fun="count", update=TRUE)
plot(rr_points)
length(rr_points[which(rr_points[]>0)])
################################################################################

res.pca = PCA(X_FST[,-c(1,2)],scale.unit = T)
x_H <- HCPC(res.pca, nb.clust = -1, min = 3, max = NULL, graph = TRUE)
fviz_dend(x_H, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
fviz_pca_var(res.pca, col.var = "black")


fviz_pca_biplot(res.pca, geom.ind = "point",
                col.ind = x_H$data.clust$clust, # color by groups
                #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                addEllipses = TRUE, ellipse.type = "convex",
                legend.title = "Groups",col.var = "black",
)

#check variance
apply(X_FST[,-c(1,2)], 2, var) == 0

X_FST$CLUSTER <- x_H$data.clust$clust
write.table(X_FST,paste0(out_dir,"/","Rattus_rattus_data.tsv"),sep = "\t",na = "",row.names = F)
################################################################################
unique_coords <- unique(paste(X_FST$Longitud,"_",X_FST$Latitud))
################################################################################
#plotting points
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
sp_NA <- st_bbox(c(xmin=xmin-0.5, xmax=xmax+0.5, ymax=ymax+0.5, ymin=ymin-0.5), crs = st_crs(World))

sp_NA2 <- sp_NA


#st_crs(coords) <- st_crs(World)

tmap_options(check.and.fix = TRUE)
tm_map_points <- 
  tm_shape(World)+#,bbox = sp_NA) +  
  tm_fill(col="gray91",legend.show = F) +
  #Rattus rattus
  tm_shape(coords)+#,bbox = sp_NA)+ 
  tm_dots(legend.show = F,col = "red",size = 4) + #purple
  #Rattus novergicus
  tm_shape(World,bbox = sp_NA) +  tm_borders("black") +
  
  tm_facets(nrow = 1, sync = TRUE)+
  tm_layout(inner.margins=0,
            legend.text.size=10,
            legend.title.size=10,
            legend.position = c("left","bottom"),
            legend.bg.color = "white", legend.bg.alpha=.2)#;gc()

tmap_save(filename=paste0(out_dir,"/GRAPHICS/","POINTS_RR_RF",".pdf"),tm=tm_map_points,dpi=900,width =100,height=100,units = "cm");gc()

###############################################################################
tm_map_points2 <- 
  tm_shape(World,bbox = sp_NA) +  
  tm_fill(col="gray91",legend.show = F) +
  #Rattus rattus
  tm_shape(coords)+#,bbox = sp_NA)+ 
  tm_dots(legend.show = T,col = "CLUSTER",size = 4,alpha=1) + #purple
  #Rattus novergicus
  tm_shape(World,bbox = sp_NA) +  tm_borders("black") +
  
  tm_facets(nrow = 1, sync = TRUE)+
  tm_layout(inner.margins=0,
            legend.text.size=5,
            legend.title.size=5,
            legend.position = c("right","bottom"),
            legend.bg.color = "white", legend.bg.alpha=1)#;gc()

tmap_save(filename=paste0(out_dir,"/GRAPHICS/","POINTS_RR_RF_CLUSTER",".pdf"),tm=tm_map_points2,dpi=900,width =100,height=100,units = "cm");gc()
################################################################################
##https://rstudio-pubs-static.s3.amazonaws.com/519785_adeb1a4592584939986b1b23103d213e.html
# Get the existing method for the ranger() function
#https://brunaw.com/slides/rladies-dublin/RF/intro-to-rf.html#9
corrplot::corrplot(cor(X_FST[,c(4:8)]), method = 'square', type = 'lower', insig='blank',
                   addCoef.col ='black', number.cex = 0.8, diag=FALSE)#,
set.seed(1)

num_trees = seq(0,10000,50)
num_trees[1] <- 10
#m_try = sqrt(5)

RF_list <- list()
tuning_df <- list()

for(i in 1:length(num_trees)){
 # i <- 1
  x <- as.data.frame(matrix(ncol=5,nrow=1))
  colnames(x) <- c("mtry","ntree","OOB","RMSE","i")
  first_rf <- ranger(Ho ~ GDP_per_capita_PPP_2015 +
                     gpw_v4_population_density_2020 +
                     povmap.grdi.v1 + 
                     travel_time_to_cities_12_MOD +
                     travel_time_to_ports_5_MOD, 
                   num.trees = num_trees[i], #mtry = m_try, #
                   importance = "impurity",
                   data = X_FST,oob.error = T,num.threads = 4,replace = F,seed = 1)

  RF_list[[i]] <- first_rf
  x$mtry <-  first_rf$mtry
  x$ntree <- num_trees[i]
  x$OOB <- first_rf$r.squared 
  x$RMSE <- first_rf$prediction.error
  x$i <- i
  tuning_df[[i]] <- x
  rm(x)
  rm(first_rf)
};rm(i)

tuning_df <- do.call(rbind,tuning_df)
tuning_df <- tuning_df[order(tuning_df$OOB,decreasing = T),]
x <- predict(RF_list[tuning_df$i[1]][1], X_FST[,c(4:8)])
################################################################################
library(ggpmisc);library(ggplot2)
x_o <- data.frame(Observed=X_FST$Ho,predicted=x[[1]]$predictions)
# using default formula, label and methods
ggplot(data = x_o, aes(x = Observed, y = predicted)) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point()
################################################################################
Metrics::rmse(X_FST$Ho,x[[1]]$predictions)
plot(X_FST$Ho,x[[1]]$predictions)
res <- caret::postResample(x[[1]]$predictions, X_FST$Ho)
res

rm(x)

imps <- data.frame(var = colnames(X_FST[,c(4:8)]),
                   imps = RF_list[tuning_df$i[1]][[1]]$variable.importance/max(RF_list[tuning_df$i[1]][[1]]$variable.importance))
imps %>% 
  ggplot(aes(imps, x = reorder(var, imps))) +
  geom_point(size = 10, colour = "#ff6767") +
  coord_flip() +
  labs(x = "Predictors", y = "Importance scores") +
  theme_bw(18)

################################################################################
#subsetting to only cells with values (reduce computing time)
sample_p2 <- sample_p
sample_p2[] <- 1:ncell(sample_p2)
#dajusting to poverty raster
sample_p2 <- sample_p2*sample_p
#plotting cell ids
#plot(sample_p2)
# only getting cells with values to obtain predictions
#
sample_p2_cells <- sample_p2[]
sample_p2_cells <- sample_p2_cells[which(!is.na(sample_p2_cells))]

###############################################################################}
#predict using raster
temp.dt <- sample_p[]
temp.dt2 <- temp.dt
#joining all raster values in a ff object
temp.dt <- ff(1,dim=c(ncell(temp.dt),dim(layers3)[3]-1),vmode="double")

#filling out each column
lapply(2:dim(layers3)[3],function(i){
  t <- getValues(layers3[[i]]) 
  temp.dt[,i-1] <- t[]
  return(cat(names(layers3)[[i]]," done\n"))})

#subsetting to cells with values
temp.dt <- as.data.frame(temp.dt[sample_p2_cells,])
temp.dt$cell <- sample_p2_cells
#subsetting ff object to dataframe and only keeping cells with complete values for the 10 predictors
temp.dt <- temp.dt[which(complete.cases(temp.dt)),]
colnames(temp.dt)[1:(ncol(temp.dt)-1)] <- names(layers3)[2:6]
temp.dt$prediction <- NA

#save.image(paste0("D:/DRYAD/RATTUS/",VAR,"_RN",".RData"))

#predicting by chunks to ommit memory problems (chunks of 100 thousands records)
x_chunks <- chunks(1,  nrow(temp.dt), by=100000) #200000
for(i in 1:length(x_chunks)){
  #i <- 1
  message(paste(round(i/length(x_chunks)*100,2)," %"))
  x_i1 <- as.numeric(as.character(x_chunks[[i]])[1])
  x_i2 <- as.numeric(as.character(x_chunks[[i]])[2])
  
  x <- predict(RF_list[tuning_df$i[1]][1], temp.dt[x_i1:x_i2,c(1:5)])[[1]]$predictions

  temp.dt$prediction[x_i1:x_i2] <- x
  rm(x)  
};rm(i)

#saving an empty raster
sample_p2_cells3 <- sample_p
sample_p2_cells3[which(!is.na(sample_p2_cells3[]))] <- NA
sample_p2_cells3[temp.dt$cell] <- temp.dt$prediction


sp_NA <-  extent(c(xmin=-180,xmax=180,ymin=-60,ymax=80))

  
Map3 <- tm_shape(sample_p2_cells3,bbox = sp_NA,raster.downsample=F) + 
  tm_raster(style= "equal", n=5, palette=rev(get_brewer_pal("PuOr", n = 5, plot=FALSE)), title="Observed heterozygosity")+
  
  tm_layout(legend.outside=TRUE, legend.outside.position="right") +
  tm_shape(World,bbox = sp_NA) + 
  tm_grid(lines=T) +
  tm_borders("black") +
  tm_facets(nrow = 1, sync = TRUE)+
  tm_layout(inner.margins=0,
            legend.text.size=10,
            legend.title.size=10,
            legend.position = c("left","bottom"),
            legend.bg.color = "white", legend.bg.alpha=.2)#;gc()
#Map3

tmap_save(filename=paste0(out_dir,"/GRAPHICS/","Ho_RR",".pdf"),tm=Map3,dpi=900,width =180,height=100,units = "cm");gc()


sp_NA2[1] <- sp_NA2[1] - 2
sp_NA2[2] <- sp_NA2[2] - 2
sp_NA2[3] <- sp_NA2[3] + 2
sp_NA2[4] <- sp_NA2[4] + 2
Map4 <- tm_shape(sample_p2_cells3,bbox = sp_NA2,raster.downsample=F) + 
  tm_raster(style= "equal", n=5, palette=rev(get_brewer_pal("PuOr", n = 5, plot=FALSE)), title="Observed heterozygosity")+
  tm_shape(coords)+#,bbox = sp_NA)+ 
  tm_dots(legend.show = T,col = "black",size = 4,alpha=1) + #purple
  
  tm_layout(legend.outside=TRUE, legend.outside.position="right") +
  tm_shape(World,bbox = sp_NA2) + 
  tm_grid(lines=T) +
  tm_borders("black") +
  tm_facets(nrow = 1, sync = TRUE)+
  tm_layout(inner.margins=0,
            legend.text.size=10,
            legend.title.size=10,
            legend.position = c("left","bottom"),
            legend.bg.color = "white", legend.bg.alpha=.2)#;gc()

tmap_save(filename=paste0(out_dir,"/GRAPHICS/","Ho_RR_POINTS",".pdf"),tm=Map4,dpi=900,width =180,height=100,units = "cm");gc()


writeRaster(sample_p2_cells3,paste0("D:/DRYAD/RATTUS/",VAR,"_RR",".tif"),overwrite=T)

save.image(paste0("D:/DRYAD/RATTUS/",VAR,"_RR",".RData"))

VAR <- "Ho"
load(paste0("D:/DRYAD/RATTUS/",VAR,"_RR",".RData"))


mess_RR <- dismo::mess(x = layers3[[2:6]],X_FST[,c(4:8)])
writeRaster(mess_RR,paste0("D:/DRYAD/RATTUS/","MESS_RR",".tif"),overwrite=T)

mess_RR2 <- mess_RR
mess_RR2[which(mess_RR2[]<=0)] <- NA
#mess_RR2[which(mess_RR2[]>0)] <- 1

qM <- quantile(mess_RR2,0.5)

mess_RR2[which(mess_RR2[]<=qM)] <- NA
mess_RR2[which(mess_RR2[]>qM)] <- 1
#plot(mess_RR2)

sample_p2_cells3_FILT <- sample_p2_cells3
sample_p2_cells3_FILT <- sample_p2_cells3_FILT*mess_RR2
#plot(sample_p2_cells3_FILT)


Map3 <- tm_shape(sample_p2_cells3_FILT,bbox = sp_NA,raster.downsample=F) + 
  tm_raster(style= "equal", n=5, palette=rev(get_brewer_pal("PuOr", n = 5, plot=FALSE)), title="Observed heterozygosity")+
  
  tm_layout(legend.outside=TRUE, legend.outside.position="right") +
  tm_shape(World,bbox = sp_NA) + 
  tm_grid(lines=T) +
  tm_borders("black") +
  tm_facets(nrow = 1, sync = TRUE)+
  tm_layout(inner.margins=0,
            legend.text.size=10,
            legend.title.size=10,
            legend.position = c("left","bottom"),
            legend.bg.color = "white", legend.bg.alpha=.2)#;gc()
#Map3

tmap_save(filename=paste0(out_dir,"/GRAPHICS/","Ho_RR_MESS",".pdf"),tm=Map3,dpi=900,width =180,height=100,units = "cm");gc()


sp_NA2[1] <- sp_NA2[1] - 2
sp_NA2[2] <- sp_NA2[2] - 2
sp_NA2[3] <- sp_NA2[3] + 2
sp_NA2[4] <- sp_NA2[4] + 2

Map4 <- tm_shape(sample_p2_cells3_FILT,bbox = sp_NA2,raster.downsample=F) + 
  tm_raster(style= "equal", n=5, palette=rev(get_brewer_pal("PuOr", n = 5, plot=FALSE)), title="Observed heterozygosity")+
  tm_shape(coords)+#,bbox = sp_NA)+ 
  tm_dots(legend.show = T,col = "black",size = 4,alpha=1) + #purple
  
  tm_layout(legend.outside=TRUE, legend.outside.position="right") +
  tm_shape(World,bbox = sp_NA2) + 
  tm_grid(lines=T) +
  tm_borders("black") +
  tm_facets(nrow = 1, sync = TRUE)+
  tm_layout(inner.margins=0,
            legend.text.size=10,
            legend.title.size=10,
            legend.position = c("left","bottom"),
            legend.bg.color = "white", legend.bg.alpha=.2)#;gc()

tmap_save(filename=paste0(out_dir,"/GRAPHICS/","Ho_RR_POINTS_MESS",".pdf"),tm=Map4,dpi=900,width =180,height=100,units = "cm");gc()

