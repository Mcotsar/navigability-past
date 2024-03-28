library(udunits2)
library(sf)
library(terra)
library(gridExtra)
library(raster)
library(osmdata) 
library(data.table)
library(ggplot2)
library(mcp)
library(changepoint)
library(segmented)
library(patchwork)
library(dplyr)
library(FactoMineR)
library(ggmap)

###############################################################
############# LOADING ALL FILES AND CLEANING THEM ############# 
###############################################################


### empirical points
comp_points <- st_read("Endpoints_empirical.shp")
comp_points <- st_transform(x = comp_points,crs = 4326)
comp_points <- comp_points[order(comp_points$Nom.riv),]


### computing confluence points
conf_points <- st_read("MinervaRiv.shp")
conf_points <- conf_points[order(conf_points$Name),]
conf_points <- st_zm(x = conf_points,drop = TRUE)
X <- st_intersection(x = conf_points[1,],y = conf_points[2,])
for(j in 1:nrow(conf_points))
{
  for(i in 1:nrow(conf_points))
  {
    x <- st_intersection(x = conf_points[i,],y = conf_points[j,])
    if(nrow(x))
    {X <- rbind(X,x)}
  }
}
X <- X[which(X$Name != X$Name.1),]
X <- X[which(duplicated(X$geometry) == FALSE),]
conf_points <- X
conf_points <- st_transform(conf_points,crs = 4326)
conf_points$Name[which(conf_points$Name == "Morava West")] <- "Morava west"
rm(x,X,i,j)

# loading rivers
raster = rast("SRTM romanFull.tif")
rivers = st_read("MinervaRiv_FullSections.shp")
rivers$Name[is.na(rivers$Name)] <- "Aude"
# aude = st_read("aude_review.shp")
# aude <- st_transform(x = aude,crs = 3035)
# colnames(aude)[1] <- "Name"
# aude$Name <- "Aude"
# rivers <- rbind(rivers,aude)
rivers <- rivers[order(rivers$Name),]
rivers$Name[1] <- "Rhône"
rivers <- rivers[order(rivers$Name),]
rivers <- st_zm(x = rivers,drop = TRUE)
rivers <- st_cast(x = rivers,to = "LINESTRING",do_split = TRUE)
if(any(duplicated(rivers$Name)) == TRUE)
{
  cat(paste("There seems to be problem in the input shapefile.",
  "The conversion from MULTILINESTRING to LINESTRING produced sub-segments of a river.",
  "Pleaser review your input shapefile with QGIS or R.",
  paste("Here is the list of rivers presenting problems :",paste(unique(rivers$Name[which(duplicated(rivers$Name))]),collapse = " ")),
  sep = "\n"))
}else
{
  cat("There seems to be no problem in the input shapefile.")
}

# computing a point every 100m in all rivers
distance <- units::set_units(100, "m")
points <- st_line_sample(rivers, density = distance)
points <- st_transform(points,4326)


#retrieving elevation from the raster for each points
elevation <- extract(x = raster,y = st_coordinates(points)[,1:2])


# creation of the table with the point coordinates and their elevation
table_points <- st_coordinates(points)
table_points <- data.table(table_points)
how_much_points_per_river <- as.numeric(summary(as.factor(table_points$L1)))
table_points$name <- NA
table_points$name <- rep(rivers$Name,how_much_points_per_river) 
table_points$Z <- elevation
correction <- table_points[table_points$name == "Aude",]
correction <- correction[order(correction$Z,decreasing = TRUE),]
table_points[table_points$name == "Aude",] <- correction
rm(correction)


##############################################################
################## PLOTTING ELEVATION CURVE ################## 
################## WITH ELBOW METHOD RESULT ##################
##################    WITHOUT CONFLUENCE    ##################
##############################################################


# number of rivers we are working on
n_rivers <- max(table_points$L1)
# creating table RES in which we will store all results
# we will store the coordinates and the altitude of the 
# estimated change points and their estimation bounds. 
RES <- as.data.frame(matrix(0,nrow = 3*3*n_rivers,ncol = 3))
RES <- st_as_sf(RES,coords = c("V1","V2","V3"))
# We store for each estimation the name of the river
RES$rivers <- rep(rivers$Name,each = 3*3)
# We store the model type (was it a single or two CP model)
RES$model_type <- rep(rep(c("M1","M2","M2"),n_rivers),each = 3)
# We indicated the line corresponds to the estimated value
# the lower bound value or the upper bound value
RES$stat <- rep(c("Est","Low","Up"),n_rivers*3)
st_crs(RES) <- 4326
for(i in 1:n_rivers)
{
  # storing the number of points in the river
  n <- sum(table_points$L1 == i)
  # finding all points belonging to the river i in table_points
  pos <- which(table_points$L1 == i)
  # creation of a table with two columns storing 
  # the elevation of the river and the distance from the source
  tab <- cbind(table_points[pos,"Z"],0:(length(pos)-1)*100)
  tab <- data.table(tab)
  colnames(tab) <- c("Elevation","Distance")
  # performing change point detection with
  my.lm <- lm(formula = Elevation~Distance,data = tab)
  ## 1 changepoint (model 1)
  seg.lm1 <- segmented(obj = my.lm,seg.Z = ~Distance,npsi = 1)
  ## 2 changepoints (model 2)
  seg.lm2 <- segmented(obj = my.lm,seg.Z = ~Distance,npsi = 2)
  # retrieving the estimated changepoint in model 1
  res <- as.numeric(round(seg.lm1$psi[,2]/100)+1)
  res_low <- as.numeric(round(confint(seg.lm1)[,2]/100)+1)
  res_up <- as.numeric(round(confint(seg.lm1)[,3]/100)+1)
  # retrieving the coordinates of this changepoint in table_points
  res_p <- st_as_sf(table_points[pos[res],],coords = c("X","Y","Z"))
  res_up <- st_as_sf(table_points[pos[res_up],],coords = c("X","Y","Z"))
  res_low <- st_as_sf(table_points[pos[res_low],],coords = c("X","Y","Z"))
  # storing it into RES
  RES$geometry[which(RES$rivers == res_p$name & RES$model_type == "M1" & RES$stat == "Est")] <- res_p$geometry
  RES$geometry[which(RES$rivers == res_p$name & RES$model_type == "M1" & RES$stat == "Low")] <- res_low$geometry
  RES$geometry[which(RES$rivers == res_p$name & RES$model_type == "M1" & RES$stat == "Up")] <- res_up$geometry
  # retrieving the estimated changepoints in model 2
  res <- as.numeric(round(seg.lm2$psi[,2]/100)+1)
  res_low <- as.numeric(round(confint(seg.lm2)[,2]/100)+1)
  res_up <- as.numeric(round(confint(seg.lm2)[,3]/100)+1)
  # retrieving the coordinates of this changepoint in table_points
  res_p <- st_as_sf(table_points[pos[res],],coords = c("X","Y","Z"))
  res_up <- st_as_sf(table_points[pos[res_up],],coords = c("X","Y","Z"))
  res_low <- st_as_sf(table_points[pos[res_low],],coords = c("X","Y","Z"))
  # storing it into RES
  RES$geometry[which(RES$rivers == res_p$name & RES$model_type == "M2" & RES$stat == "Est")] <- res_p$geometry
  RES$geometry[which(RES$rivers == res_p$name & RES$model_type == "M2" & RES$stat == "Low")] <- res_low$geometry
  RES$geometry[which(RES$rivers == res_p$name & RES$model_type == "M2" & RES$stat == "Up")] <- res_up$geometry
}


################################################################
###################### PLOT GENERATION :  ######################
###################### WITHOUT CONFLUENCE ######################
################################################################

# making the names of the rivers match between the empirical
# points and the imported shapefile 
comp_points$Nom.riv[comp_points$Nom.riv == "Meuse"] <- "Maas"
comp_points$Nom.riv[comp_points$Nom.riv == "Rhin"] <- "Rhine"
comp_points$Nom.riv[comp_points$Nom.riv == "Moselle"] <- "Mosel"
# creating a spatial object from table_points
table_points2 <- st_as_sf(x = table_points,coords = c("X","Y"))
st_crs(table_points2) <- 4326
# cleaning the work environment a little
rm(elevation,raster,rivers,i,n,pos,res,my.lm,res_p,seg.lm1,seg.lm2,tab,distance,res_low,res_up,aude)
# creating a list that will hold several information on each river
L<-as.list(rep(NA,n_rivers))
# for each river
for(i in 1:n_rivers)
{
  # we store 13 information on each river
  L[[i]] <- as.list(rep(NA,14))
  # pos <- which(table_points$L1 == i)
  # we keep the position of points belonging to that river in table_points
  L[[i]][[1]] <- which(table_points$L1 == i)
  name_river <- unique(table_points$name[L[[i]][[1]]])
  n <- sum(table_points$L1 == i)
  tab <- cbind(table_points[L[[i]][[1]],"Z"],0:(length(L[[i]][[1]])-1)*100)
  tab <- data.table(tab)
  tab2 <- cbind(table_points[L[[i]][[1]],"X"],table_points[L[[i]][[1]],"Y"])
  tab2 <- data.table(tab2)
  colnames(tab) <- c("Elevation","Distance")
  # we store the table with the elevation VS the distance in the river (tab)
  # we store the table with the GPS coordinates of of the points in the river (tab2)
  L[[i]][[2]] <- tab
  L[[i]][[3]] <- tab2
  my.lm <- lm(formula = Elevation~Distance,data = L[[i]][[2]])
  seg.lm1 <- segmented(obj = my.lm,seg.Z = ~Distance,npsi = 1)
  seg.lm2 <- segmented(obj = my.lm,seg.Z = ~Distance,npsi = 2)
  my.fitted1 <- fitted(seg.lm1)
  my.fitted2 <- fitted(seg.lm2)
  my.model1 <- data.frame(Distance = tab$Distance, Elevation = my.fitted1)
  my.model2 <- data.frame(Distance = tab$Distance, Elevation = my.fitted2)
  # we store the results of the changepoint models (with one and two changepoints)
  L[[i]][[4]] <- my.model1
  L[[i]][[5]] <- my.model2
  # we store the elevation and distance from the source of the changepoint found in model 1
  L[[i]][[6]] <- which(apply(X = L[[i]][[3]],MARGIN = 1,FUN = "paste",collapse = "_") %in%
                         paste(st_coordinates(RES)[RES$rivers == name_river & RES$model_type=="M1" & RES$stat=="Est",1:2],collapse = "_"))
  L[[i]][[7]] <- which(apply(X = L[[i]][[3]],MARGIN = 1,FUN = "paste",collapse = "_") %in%
                          apply(X = st_coordinates(RES)[RES$rivers == name_river & RES$model_type=="M1" & RES$stat %in% c("Low","Up"),1:2],MARGIN = 1,FUN = "paste",collapse = "_"))
  # we store the elevation and distance from the source of the changepoints found in model 2
  L[[i]][[8]] <- which(apply(X = L[[i]][[3]],MARGIN = 1,FUN = "paste",collapse = "_") %in%
                         apply(X = st_coordinates(RES)[RES$rivers == name_river & RES$model_type=="M2" & RES$stat == "Est",1:2],MARGIN = 1,FUN = "paste",collapse = "_"))
  L[[i]][[9]] <- which(apply(X = L[[i]][[3]],MARGIN = 1,FUN = "paste",collapse = "_") %in%
                          apply(X = st_coordinates(RES)[RES$rivers == name_river & RES$model_type=="M2" & RES$stat %in% c("Low","Up"),1:2],MARGIN = 1,FUN = "paste",collapse = "_"))
  
  if(name_river %in% unique(c(conf_points$Name,conf_points$Name.1))  & name_river %in% comp_points$Nom.riv)
  {
    # computing the position of the confluence in the elevation profile
    RW <- apply(X = conf_points,MARGIN = 1,FUN = function(x) which(x == name_river))
    RW <- as.numeric(which(lapply(RW,"length")!=0))
    for(p in 1:length(RW))
    {
      dist_max <- st_distance(x = st_cast(x = points[i],to = "POINT"),y = conf_points$geometry[RW[p]])
      dist_max <- which.min(dist_max)
      RW[p] <- dist_max
    }
    L[[i]][[13]] <- min(tab$Distance[RW])
    # computing the position of the empirical points in the elevation profile
    L[[i]][[14]] <- apply(X = st_distance(x = comp_points[comp_points$Nom.riv %in% name_river,],
                                     y = table_points2[L[[i]][[1]],]),
                     MARGIN = 1,
                     FUN = "which.min")
    # graphic of the 1 CP model
    L[[i]][[10]] <- eval(substitute(
      ggplot()+
        geom_point(data = L[[i]][[2]],aes(x = Distance,y = Elevation),size = 0.5,alpha = 0.1)+
        geom_line(data = L[[i]][[4]], aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
        geom_vline(aes(xintercept = (L[[i]][[6]]-1)*100,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)+
        geom_vline(aes(xintercept = L[[i]][[13]],color = "Confluence"),linetype = "solid")+
        geom_vline(aes(xintercept = L[[i]][[14]]*100,color = "Empirical CP"),linetype = "dashed")+
        geom_rect(aes(xmin = L[[i]][[7]][1]*100,xmax = L[[i]][[7]][2]*100,ymin = -Inf,ymax = Inf), col = "red",alpha = 0.2)+
        scale_color_manual("Information",breaks = c("Empirical CP","Confluence","Estimated CP","Linear model"),values = c("skyblue","black","red","orange"))+
        ggtitle("Model with a single change point")+theme_bw()+theme(legend.key.size = unit(1, 'cm'),legend.position = "bottom")
      ,list(i = i)))
    # graphic of the 2 CP model
    L[[i]][[11]] <- eval(substitute(
      ggplot()+
        geom_point(data = L[[i]][[2]],aes(x = Distance,y = Elevation),size = 0.5,alpha = 0.1)+
        geom_line(data = L[[i]][[5]], aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
        geom_vline(aes(xintercept = (L[[i]][[8]]-1)*100,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)+
        geom_vline(aes(xintercept = L[[i]][[13]],color = "Confluence"),linetype = "solid")+
        geom_vline(aes(xintercept = L[[i]][[14]]*100,color = "Empirical CP"),linetype = "dashed")+
        geom_rect(aes(xmin = L[[i]][[9]][1]*100,xmax = L[[i]][[9]][2]*100,ymin = -Inf,ymax = Inf), col = "red",alpha = 0.2)+
        geom_rect(aes(xmin = L[[i]][[9]][3]*100,xmax = L[[i]][[9]][4]*100,ymin = -Inf,ymax = Inf), col = "red",alpha = 0.2)+
        scale_color_manual("Information",breaks = c("Empirical CP","Confluence","Estimated CP","Linear model"),values = c("skyblue","black","red","orange"))+
        ggtitle("Model with two changepoints")+theme_bw()+theme(legend.key.size = unit(1, 'cm'),legend.position = "bottom")
      ,list(i = i)))
    # both graphics combined 
    L[[i]][[12]] <- eval(substitute((L[[i]][[10]]/L[[i]][[11]])+plot_annotation(paste(table_points$name[L[[i]][[1]][1]],"river"),theme = theme(plot.title=element_text(hjust=0.5)))+plot_layout(guides = "collect")&theme(legend.position = "bottom",text = element_text(size = 7))
                                    ,list(i=i)))
    
  }else if(name_river %in% unique(c(conf_points$Name,conf_points$Name.1)) & !name_river %in% comp_points$Nom.riv)
  {
    # computing the position of the confluence in the elevation profile
    RW <- apply(X = conf_points,MARGIN = 1,FUN = function(x) which(x == name_river))
    RW <- as.numeric(which(lapply(RW,"length")!=0))
    for(p in 1:length(RW))
    {
      dist_max <- st_distance(x = st_cast(x = points[i],to = "POINT"),y = conf_points$geometry[RW[p]])
      dist_max <- which.min(dist_max)
      RW[p] <- dist_max
    }
    L[[i]][[13]] <- min(tab$Distance[RW])
    # graphic of the 1 CP model
    L[[i]][[10]] <- eval(substitute(
      ggplot()+
        geom_point(data = L[[i]][[2]],aes(x = Distance,y = Elevation),size = 0.5,alpha = 0.1)+
        geom_line(data = L[[i]][[4]], aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
        geom_vline(aes(xintercept = (L[[i]][[6]]-1)*100,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)+
        geom_vline(aes(xintercept = L[[i]][[13]],color = "Confluence"),linetype = "solid")+
        geom_rect(aes(xmin = L[[i]][[7]][1]*100,xmax = L[[i]][[7]][2]*100,ymin = -Inf,ymax = Inf), col = "red",alpha = 0.2)+
        scale_color_manual("Information",breaks = c("Confluence","Estimated CP","Linear model"),values = c("black","red","orange"))+
        ggtitle("Model with a single change point")+theme_bw()+theme(legend.key.size = unit(1, 'cm'),legend.position = "bottom")
      ,list(i = i)))
    # graphic of the 2 CP model 
    L[[i]][[11]] <- eval(substitute(
      ggplot()+
        geom_point(data = L[[i]][[2]],aes(x = Distance,y = Elevation),size = 0.5,alpha = 0.1)+
        geom_line(data = L[[i]][[5]], aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
        geom_vline(aes(xintercept = (L[[i]][[8]]-1)*100,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)+
        geom_vline(aes(xintercept = L[[i]][[13]],color = "Confluence"),linetype = "solid")+
        geom_rect(aes(xmin = L[[i]][[9]][1]*100,xmax = L[[i]][[9]][2]*100,ymin = -Inf,ymax = Inf), col = "red",alpha = 0.2)+
        geom_rect(aes(xmin = L[[i]][[9]][3]*100,xmax = L[[i]][[9]][4]*100,ymin = -Inf,ymax = Inf), col = "red",alpha = 0.2)+
        scale_color_manual("Information",breaks = c("Confluence","Estimated CP","Linear model"),values = c("black","red","orange"))+
        ggtitle("Model with two changepoints")+theme_bw()+theme(legend.key.size = unit(1, 'cm'),legend.position = "bottom")
      ,list(i = i)))
    # both graphics combined 
    L[[i]][[12]] <- eval(substitute((L[[i]][[10]]/L[[i]][[11]])+plot_annotation(paste(table_points$name[L[[i]][[1]][1]],"river"),theme = theme(plot.title=element_text(hjust=0.5)))+plot_layout(guides = "collect")&theme(legend.position = "bottom",text = element_text(size = 7))
                                    ,list(i=i)))
  }else if(!name_river %in% unique(c(conf_points$Name,conf_points$Name.1)) & name_river %in% comp_points$Nom.riv)
  {
    # computing the position of the empirical points in the elevation profile
    L[[i]][[14]] <- apply(X = st_distance(x = comp_points[comp_points$Nom.riv %in% name_river,],
                                     y = table_points2[L[[i]][[1]],]),
                     MARGIN = 1,
                     FUN = "which.min")
    # graphic of the 1 CP model
    L[[i]][[10]] <- eval(substitute(
      ggplot()+
        geom_point(data = L[[i]][[2]],aes(x = Distance,y = Elevation),size = 0.5,alpha = 0.1)+
        geom_line(data = L[[i]][[4]], aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
        geom_vline(aes(xintercept = (L[[i]][[6]]-1)*100,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)+
        geom_vline(aes(xintercept = L[[i]][[14]]*100,color = "Empirical CP"),linetype = "dashed")+
        geom_rect(aes(xmin = L[[i]][[7]][1]*100,xmax = L[[i]][[7]][2]*100,ymin = -Inf,ymax = Inf), col = "red",alpha = 0.2)+
        scale_color_manual("Information",breaks = c("Empirical CP","Estimated CP","Linear model"),values = c("skyblue","red","orange"))+
        ggtitle("Model with a single change point")+theme_bw()+theme(legend.key.size = unit(1, 'cm'),legend.position = "bottom")
      ,list(i = i)))
    # graphic of the 2 CP model
    L[[i]][[11]] <- eval(substitute(
      ggplot()+
        geom_point(data = L[[i]][[2]],aes(x = Distance,y = Elevation),size = 0.5,alpha = 0.1)+
        geom_line(data = L[[i]][[5]], aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
        geom_vline(aes(xintercept = (L[[i]][[8]]-1)*100,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)+
        geom_vline(aes(xintercept = L[[i]][[14]]*100,color = "Empirical CP"),linetype = "dashed")+
        geom_rect(aes(xmin = L[[i]][[9]][1]*100,xmax = L[[i]][[9]][2]*100,ymin = -Inf,ymax = Inf), col = "red",alpha = 0.2)+
        geom_rect(aes(xmin = L[[i]][[9]][3]*100,xmax = L[[i]][[9]][4]*100,ymin = -Inf,ymax = Inf), col = "red",alpha = 0.2)+
        scale_color_manual("Information",breaks = c("Empirical CP","Estimated CP","Linear model"),values = c("skyblue","red","orange"))+
        ggtitle("Model with two changepoints")+theme_bw()+theme(legend.key.size = unit(1, 'cm'),legend.position = "bottom")
      ,list(i = i)))
    # both graphics combined
    L[[i]][[12]] <- eval(substitute((L[[i]][[10]]/L[[i]][[11]])+plot_annotation(paste(table_points$name[L[[i]][[1]][1]],"river"),theme = theme(plot.title=element_text(hjust=0.5)))+plot_layout(guides = "collect")&theme(legend.position = "bottom",text = element_text(size = 7))
                                    ,list(i=i)))
  }else
  {
    # graphic of the 1 CP model
    L[[i]][[10]] <- eval(substitute(
      ggplot()+
        geom_point(data = L[[i]][[2]],aes(x = Distance,y = Elevation),size = 0.5,alpha = 0.1)+
        geom_line(data = L[[i]][[4]], aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
        geom_vline(aes(xintercept = (L[[i]][[6]]-1)*100,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)+
        geom_rect(aes(xmin = L[[i]][[7]][1]*100,xmax = L[[i]][[7]][2]*100,ymin = -Inf,ymax = Inf), col = "red",alpha = 0.2)+
        scale_color_manual("Information",breaks = c("Estimated CP","Linear model"),values = c("red","orange"))+
        ggtitle("Model with a single change point")+theme_bw()+theme(legend.key.size = unit(1, 'cm'),legend.position = "bottom")
      ,list(i = i)))
    # graphic of the 2 CP model
    L[[i]][[11]] <- eval(substitute(
      ggplot()+
        geom_point(data = L[[i]][[2]],aes(x = Distance,y = Elevation),size = 0.5,alpha = 0.1)+
        geom_line(data = L[[i]][[5]], aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
        geom_vline(aes(xintercept = (L[[i]][[8]]-1)*100,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)+
        geom_rect(aes(xmin = L[[i]][[9]][1]*100,xmax = L[[i]][[9]][2]*100,ymin = -Inf,ymax = Inf), col = "red",alpha = 0.2)+
        geom_rect(aes(xmin = L[[i]][[9]][3]*100,xmax = L[[i]][[9]][4]*100,ymin = -Inf,ymax = Inf), col = "red",alpha = 0.2)+
        scale_color_manual("Information",breaks = c("Estimated CP","Linear model"),values = c("red","orange"))+
        ggtitle("Model with two changepoints")+theme_bw()+theme(legend.key.size = unit(1, 'cm'),legend.position = "bottom")
      ,list(i = i)))
    # both graphics combined
    L[[i]][[12]] <- eval(substitute((L[[i]][[10]]/L[[i]][[11]])+plot_annotation(paste(table_points$name[L[[i]][[1]][1]],"river"),theme = theme(plot.title=element_text(hjust=0.5)))+plot_layout(guides = "collect")&theme(legend.position = "bottom",text = element_text(size = 7))
                                    ,list(i=i)))
  }
}


# Saving the plots in a single PDF
pdf(file = "fullsection_wo_constraint.pdf",onefile = TRUE,paper = "a4")
for(i in 1:n_rivers)
{
  print(L[[i]][[12]])
}
dev.off()


# Saving the table containing the results
RES <- st_transform(x = RES,crs = 3035)
RES$X <- st_coordinates(x = RES$geometry)[,1]
RES$Y <- st_coordinates(x = RES$geometry)[,2]
RES$alt <- st_coordinates(x = RES$geometry)[,3]
RES$id <- NA
for(i in 1:nrow(RES))
{
  if(RES$model_type[i] == "M1")
  {
    RES$id[i] <- "Only_point"
  }else if(is.na(RES$id[i]) & RES$model_type[i] == "M2")
  {
    if(RES$alt[i] > RES$alt[i+3])
    {
      RES$id[i] <- "Upper_point"
      RES$id[i+3] <- "Lower_point"
    }else
    {
      RES$id[i] <- "Lower_point"
      RES$id[i+3] <- "Upper_point"
    }
  }
}
save(RES,L,table_points,comp_points,conf_points,file = "fullsection_wo_constraint.Rdata")
# save(RES,file = "fullsection_wo_constraint.Rdata")
write.table(x = RES,file = "fullsection_wo_constraint.csv",sep = ";")


################################################################
###################### INDICATORS OF PERF ######################
######################    COMPUTATION     ######################
################################################################


# these indicators can only be computed on rivers 
# that have empirical data points of navigabilty at disposal
# keeping these rivers
L_key <- which(unique(RES$rivers) %in% comp_points$Nom.riv)
river_names <- unique(RES$rivers)[L_key]
# constructing the table
tab_indic <- as.data.table(matrix(NA,nrow = length(unique(river_names)),ncol = 18)) 
colnames(tab_indic) <- c("Name",
                         "Prec_2up_M1","Prec_2low_M1","Pts_order_M1",
                         "Prec_up2up_M2","Prec_low2low_M2","Overlap_est_over_emp_M2","Overlap_emp_over_est_M2","Pts_order_M2",
                         "pos_up_emp","pos_low_emp","pos_M1","pos_up_M2","pos_low_M2",
                         "char_pos_conf","char_max_height","char_max_length","char_AUC")
tab_indic$Name <- river_names
for(i in 1:nrow(tab_indic))
{
  # temp_res <- RES[which(RES$rivers == tab_indic$Name[i]),]
  # temp_tab <- table_points[which(table_points$name == tab_indic$Name[i]),]
  # temp_comp <- comp_points[which(comp_points$Nom.riv == tab_indic$Name[i]),]
  tab_indic$Prec_2up_M1[i] <- abs(L[[L_key[i]]][[14]][1]-L[[L_key[i]]][[6]])*100
  tab_indic$Prec_2low_M1[i] <- abs(L[[L_key[i]]][[14]][2]-L[[L_key[i]]][[6]])*100
  tab_indic$Pts_order_M1[i] <- paste(c("Up_emp","Low_emp","CP1")[order(c(L[[L_key[i]]][[14]][1],L[[L_key[i]]][[14]][2],L[[L_key[i]]][[6]]))],collapse = "->")
  tab_indic$Prec_up2up_M2[i] <- abs(L[[L_key[i]]][[14]][1]-L[[L_key[i]]][[8]][1])*100
  tab_indic$Prec_low2low_M2[i] <- abs(L[[L_key[i]]][[14]][2]-L[[L_key[i]]][[8]][2])*100
  tab_indic$Overlap_est_over_emp_M2[i] <- round(length(which(L[[L_key[i]]][[14]][1]:L[[L_key[i]]][[14]][2] %in% L[[L_key[i]]][[8]][1]:L[[L_key[i]]][[8]][2]))/length(L[[L_key[i]]][[14]][1]:L[[L_key[i]]][[14]][2])*100,2)
  tab_indic$Overlap_emp_over_est_M2[i] <- round(length(which(L[[L_key[i]]][[8]][1]:L[[L_key[i]]][[8]][2] %in% L[[L_key[i]]][[14]][1]:L[[L_key[i]]][[14]][2]))/length(L[[L_key[i]]][[8]][1]:L[[L_key[i]]][[8]][2])*100,2)
  tab_indic$Pts_order_M2[i] <- paste(c("Up_emp","Low_emp","Up_CP2","Low_CP2")[order(c(L[[L_key[i]]][[14]][1],L[[L_key[i]]][[14]][2],L[[L_key[i]]][[8]][1],L[[L_key[i]]][[8]][2]))],collapse = "->")
  tab_indic$pos_up_emp[i] <- L[[L_key[i]]][[14]][1]
  tab_indic$pos_low_emp[i] <- L[[L_key[i]]][[14]][2]
  tab_indic$pos_M1[i] <- L[[L_key[i]]][[6]]
  tab_indic$pos_up_M2[i] <- L[[L_key[i]]][[8]][1]
  tab_indic$pos_low_M2[i] <- L[[L_key[i]]][[8]][2]
  tab_indic$char_pos_conf[i] <- L[[L_key[i]]][[13]]/100
  tab_indic$char_max_height[i] <- max(L[[L_key[i]]][[2]][,1])
  tab_indic$char_max_length[i] <- L[[L_key[i]]][[2]][nrow(L[[L_key[i]]][[2]]),2]
  temp <- L[[L_key[i]]][[2]]$Elevation/tab_indic$char_max_height[i]
  res <- 0.5*(temp[-1]+temp[-length(temp)])
  n <- nrow(L[[i]][[2]])
  tab_indic$char_AUC[i] <- 1/n*sum(res)
}
# saving the results
save(tab_indic,file = "Indic_fullsection_wo_constraint.Rdata")

