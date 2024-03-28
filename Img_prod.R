library(udunits2)
library(sf)
library(terra)
library(gridExtra)
library(dtw)
library(ggmap)
library(raster)
library(ggpubr)
library(osmdata) 
library(transport)
library(data.table)
library(fda)
library(fdacluster)
library(ggplot2)
library(mcp)
library(tidyterra)
library(changepoint)
library(segmented)
library(grid)
library(patchwork)
library(dplyr)
library(FactoMineR)
library(ggmap)
library(longitudinalData)
library(rPref)
library(maotai)
library(ggmap)


# production of Figure 2

# ggmap::register_stadiamaps(key = "insert_your_registration_key_to_the_stamen_maps_API")
rivers = st_read("MinervaRiv_FullSections.shp")
rivers <- rivers[which(is.na(rivers$Name)),]
distance <- units::set_units(2500, "m")
points <- st_line_sample(rivers, density = distance)
rivers <- st_transform(rivers,crs = 4326)
points <- st_transform(points,crs = 4326)
rivers <- rivers %>% st_zm(drop = TRUE,what =  "ZM")
mylocation <- c(min(st_coordinates(rivers)[,1]),
                min(st_coordinates(rivers)[,2]),
                max(st_coordinates(rivers)[,1]),
                max(st_coordinates(rivers)[,2]))
map1 <- get_stadiamap(bbox = mylocation,
                color="bw",maptype = "stamen_terrain")
sitemap <- ggmap(map1, extent = 'device')



G1 <- sitemap+geom_sf(inherit.aes = FALSE,
                data = rivers)+ylab("Latitude")+xlab("Longitude")+
  ggtitle("(a)")

G2 <- sitemap+geom_sf(inherit.aes = FALSE,
                      data = points,
                      shape = 19, 
                      stroke = 0.3,size = 1)+
  ylab("Latitude")+
  xlab("Longitude")+
  ggtitle("(b)")

raster = rast("SRTM romanFull.tif")
subR <- crop(raster,extent(mylocation[c(1,3,2,4)]))
G3 <- G2 + geom_spatraster(data = subR,alpha = 0.2) + 
  scale_fill_gradientn("Elevation (m)", colors = topo.colors(255))+
  theme(legend.position = "bottom")+ggtitle("(c)")

load("fullsection_w_constraint.Rdata")
table_points <- table_points %>% filter(name == "Aude")
G4 <- ggplot(data = table_points[seq(1,nrow(table_points),by = 25),])+
  geom_point(aes(x = D,y = Z))+xlab("Distance in the river (m)")+ylab("Elevation (m)")+
  ggtitle("(d)")

G_FIN <- ((G1|G2|G3)/(G4)) & theme(text =  element_text(size = 7))

ggsave(filename = "Img_profiles.pdf",plot = G_FIN,device = "pdf",width = 8,height = 5,units = "in")


# production of Figure 3

load("partsection_wo_constraint.Rdata")
# pos_part <- L[[which(unique(table_points$name) == "Ill")]][[8]]
table_points_part <- table_points[table_points$name == "Ill",]
# load("fullsection_wo_constraint.Rdata")
# pos_fullwo <- L[[which(unique(table_points$name) == "Ill")]][[8]]
load("fullsection_w_constraint.Rdata")
# pos_fullw <- L[[which(unique(table_points$name) == "Ill")]][[8]]
pos_conf <- L[[which(unique(table_points$name) == "Ill")]][[13]]
table_points_full <- table_points[table_points$name == "Ill",]
rm(L,comp_points,RES,table_points)
table_points_part$D <- (1:nrow(table_points_part)-1)*100
my.lm <- lm(formula = Z~D,data = table_points_full)
seg.lm2 <- segmented(obj = my.lm,seg.Z = ~D,npsi = 2)
seg.lm7 <- segmented(obj = my.lm,seg.Z = ~D,npsi = 3)
my.lm <- lm(formula = Z~D,data = table_points_part)
seg.lm2_part <- segmented(obj = my.lm,seg.Z = ~D,npsi = 2)
my.fitted1 <- fitted(seg.lm2_part)
my.fitted2 <- fitted(seg.lm2)
my.fitted3 <- fitted(seg.lm7)
my.model1 <- data.frame(Distance = table_points_part$D, Elevation = my.fitted1)
my.model2 <- data.frame(Distance = table_points_full$D, Elevation = my.fitted2)
my.model3 <- data.frame(Distance = table_points_full$D, Elevation = my.fitted3)

G1 <- ggplot()+
  geom_point(data = table_points_full,aes(x = D,y = Z),size = 0.6,alpha = 0.1)+
  geom_line(data = my.model1, aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.75)+
  geom_vline(aes(xintercept = seg.lm2_part$psi[,2],color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.75)+
  geom_vline(aes(xintercept = pos_conf,color = "Confluence"),linetype = "solid")+
  scale_color_manual("Information",breaks = c("Empirical CP","Confluence","Estimated CP","Linear model"),values = c("skyblue","black","red","orange"))+
  ggtitle("Method 1")+theme_bw()+theme(legend.key.size = unit(1, 'cm'),legend.position = "bottom")+
  xlab("Distance")+ylab("Elevation")

G2 <- ggplot()+
  geom_point(data = table_points_full,aes(x = D,y = Z),size = 0.6,alpha = 0.1)+
  geom_line(data = my.model2, aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.75)+
  geom_vline(aes(xintercept = seg.lm2$psi[,2],color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.75)+
  geom_vline(aes(xintercept = pos_conf,color = "Confluence"),linetype = "solid")+
  scale_color_manual("Information",breaks = c("Empirical CP","Confluence","Estimated CP","Linear model"),values = c("skyblue","black","red","orange"))+
  ggtitle("Method 2")+theme_bw()+theme(legend.key.size = unit(1, 'cm'),legend.position = "bottom")+
  xlab("Distance")+ylab("Elevation")

G3 <- ggplot()+
  geom_point(data = table_points_full,aes(x = D,y = Z),size = 0.6,alpha = 0.1)+
  geom_line(data = my.model3, aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.75)+
  geom_vline(aes(xintercept = seg.lm7$psi[1:2,2],color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.75)+
  geom_vline(aes(xintercept = seg.lm7$psi[3,2]),color = "darkgrey",linetype="dashed", show.legend = F,lwd = 0.75)+
  geom_vline(aes(xintercept = pos_conf,color = "Confluence"),linetype = "solid")+
  scale_color_manual("Information",breaks = c("Empirical CP","Confluence","Estimated CP","Linear model"),values = c("skyblue","black","red","orange"))+
  ggtitle("Method 3")+theme_bw()+theme(legend.key.size = unit(1, 'cm'),legend.position = "bottom")+
  xlab("Distance")+ylab("Elevation")

G_FIN <- (G1/G2/G3)+plot_layout(guides = "collect")&theme(text = element_text(size = 7),legend.position = "bottom")

ggsave(filename = "Method_illustration.pdf",plot = G_FIN,device = "pdf",width = 8,height = 5,units = "in")



# best_of_1 <- tab[tab$meth == 1,]
# best_of_1$mean_precM1 <- (best_of_1$Prec_2up_M1+best_of_1$Prec_2low_M1)/2
# best_of_1$mean_precM2 <- (best_of_1$Prec_up2up_M2+best_of_1$Prec_low2low_M2)/2
# best_of_2 <- tab[tab$meth == 2,]
# best_of_2$mean_precM1 <- (best_of_2$Prec_2up_M1+best_of_2$Prec_2low_M1)/2
# best_of_2$mean_precM2 <- (best_of_2$Prec_up2up_M2+best_of_2$Prec_low2low_M2)/2
# best_of_3 <- tab[tab$meth == 3,]
# best_of_3$mean_precM1 <- (best_of_3$Prec_2up_M1+best_of_3$Prec_2low_M1)/2
# best_of_3$mean_precM2 <- (best_of_3$Prec_up2up_M2+best_of_3$Prec_low2low_M2)/2
# Cher M1 > M2 == M3
# Doubs M3 > M2 > M1
# Ill M2 > M3 > M1 

rm(best_of_1,best_of_2,best_of_3)

# production of Figure 4

load(file = "fullsection_w_constraint.Rdata")
pos_ill <- which(unique(RES$rivers) == "Ill")
pos_doubs <- which(unique(RES$rivers) == "Doubs")
pos_cher <- which(unique(RES$rivers) == "Cher")
tab_ill_Z <- L[[pos_ill]][[2]]
tab_doubs_Z <- L[[pos_doubs]][[2]]
tab_cher_Z <- L[[pos_cher]][[2]]
conf_ill <- L[[pos_ill]][[13]]
conf_doubs <-  L[[pos_doubs]][[13]]
conf_cher <- L[[pos_cher]][[13]]
comp_ill <- L[[pos_ill]][[14]]
comp_doubs <-  L[[pos_doubs]][[14]]
comp_cher <- L[[pos_cher]][[14]]
G_ill <- ggplot()+ 
  geom_point(data = tab_ill_Z,aes(Distance,Elevation),size = 0.5)+
  geom_vline(aes(xintercept = comp_ill*100,color = "Empirical"),linetype = "dashed",lwd = 0.7)+
  geom_vline(aes(xintercept = conf_ill,color = "Confluence"),lwd = 0.2)+
  scale_color_manual("Information",breaks = c("Empirical","Confluence","Estimated CP","Linear model"),values = c("skyblue","black","red","orange"))
G_doubs <- ggplot()+ 
  geom_point(data = tab_doubs_Z,aes(Distance,Elevation))+
  geom_vline(aes(xintercept = comp_doubs*100,color = "Empirical"),linetype = "dashed",lwd = 0.7)+
  geom_vline(aes(xintercept = conf_doubs,color = "Confluence"),lwd = 0.2)+
  scale_color_manual("Information",breaks = c("Empirical","Confluence","Estimated CP","Linear model"),values = c("skyblue","black","red","orange"))
G_cher <- ggplot()+ 
  geom_point(data = tab_cher_Z,aes(Distance,Elevation))+
  geom_vline(aes(xintercept = comp_cher*100,color = "Empirical"),linetype = "dashed",lwd = 0.7)+
  geom_vline(aes(xintercept = conf_cher,color = "Confluence"),lwd = 0.2)+
  scale_color_manual("Information",breaks = c("Empirical","Confluence","Estimated CP","Linear model"),values = c("skyblue","black","red","orange"))


# M3 plots

M3_ill_CP <- (L[[pos_ill]][[8]]-1)*100
M3_ill <- L[[pos_ill]][[5]]
M3_doubs_CP <- (L[[pos_doubs]][[8]]-1)*100
M3_doubs <- L[[pos_doubs]][[5]]
M3_cher_CP <- (L[[pos_cher]][[8]]-1)*100
M3_cher <- L[[pos_cher]][[5]]
G_ill_M3 <- G_ill+
  geom_line(data = M3_ill, aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
  geom_vline(aes(xintercept = M3_ill_CP,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)
G_doubs_M3 <- G_doubs+
  geom_line(data = M3_doubs, aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
  geom_vline(aes(xintercept = M3_doubs_CP,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)
G_cher_M3 <- G_cher+
  geom_line(data = M3_cher, aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
  geom_vline(aes(xintercept = M3_cher_CP,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)

# M2 plots
load(file = "fullsection_wo_constraint.Rdata")
M2_ill_CP <- (L[[pos_ill]][[8]]-1)*100
M2_ill <- L[[pos_ill]][[5]]
M2_doubs_CP <- (L[[pos_doubs]][[8]]-1)*100
M2_doubs <- L[[pos_doubs]][[5]]
M2_cher_CP <- (L[[pos_cher]][[8]]-1)*100
M2_cher <- L[[pos_cher]][[5]]
G_ill_M2 <- G_ill+
  geom_line(data = M2_ill, aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
  geom_vline(aes(xintercept = M2_ill_CP,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)
G_doubs_M2 <- G_doubs+
  geom_line(data = M2_doubs, aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
  geom_vline(aes(xintercept = M2_doubs_CP,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)
G_cher_M2 <- G_cher+
  geom_line(data = M2_cher, aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
  geom_vline(aes(xintercept = M2_cher_CP,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)


# M1 plots
load(file = "partsection_wo_constraint.Rdata")
M1_ill_CP <- (L[[pos_ill]][[8]]-1)*100
M1_ill <- L[[pos_ill]][[5]]
M1_doubs_CP <- (L[[pos_doubs]][[8]]-1)*100
M1_doubs <- L[[pos_doubs]][[5]]
M1_cher_CP <- (L[[pos_cher]][[8]]-1)*100
M1_cher <- L[[pos_cher]][[5]]
G_ill_M1 <- G_ill+
  geom_line(data = M1_ill, aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
  geom_vline(aes(xintercept = M1_ill_CP,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)
G_doubs_M1 <- G_doubs+
  geom_line(data = M1_doubs, aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
  geom_vline(aes(xintercept = M1_doubs_CP,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)
G_cher_M1 <- G_cher+
  geom_line(data = M1_cher, aes(x = Distance, y = Elevation,color = "Linear model"),linetype = "solid",lwd = 0.5)+
  geom_vline(aes(xintercept = M1_cher_CP,color = "Estimated CP"),linetype="dashed", show.legend = F,lwd = 0.5)



row1 <- ggplot() + annotate(geom = "text", x=1, y=1, label= 'Method 1', angle = 90) + theme_void() 
row2 <- ggplot() + annotate(geom = "text", x=1, y=1, label= 'Method 2', angle = 90) + theme_void() 
row3 <- ggplot() + annotate(geom = "text", x=1, y=1, label= 'Method 3', angle = 90) + theme_void() 
col1 <- ggplot() + annotate(geom = "text", x=1, y=1, label= 'Cher') + theme_void() 
col2 <- ggplot() + annotate(geom = "text", x=1, y=1, label= 'Doubs') + theme_void() 
col3 <- ggplot() + annotate(geom = "text", x=1, y=1, label= 'Ill') + theme_void() 
layoutplot <- "
#dddeeefff
aggghhhiii
aggghhhiii
aggghhhiii
bjjjkkklll
bjjjkkklll
bjjjkkklll
cmmmnnnooo
cmmmnnnooo
cmmmnnnooo
"
plotlist <- list(a = row1, b = row2, c = row3, 
                 d = col1, e = col2, f = col3, 
                 g = G_cher_M1,h = G_doubs_M1,i = G_ill_M1,
                 j = G_cher_M2,k = G_doubs_M2,l = G_ill_M2,
                 m = G_cher_M3,n = G_doubs_M3,o = G_ill_M3)

G_FIN <- wrap_plots(plotlist, guides = 'collect', design = layoutplot)&theme(text = element_text(size = 7),legend.position = "bottom")

ggsave(filename = "Method_results.pdf",plot = G_FIN,device = "pdf",width = 8,height = 7,units = "in")

# production of Figure 1

load("fullsection_w_constraint.Rdata")
table_eure <- table_points %>% filter(name %in% c("Eure"))
table_rhine <- table_points %>% filter(name %in% c("Rhine"))
G1 <- ggplot(data = table_eure[seq(1,nrow(table_eure),length.out = 1000),])+
  geom_point(aes(x = D,y = Z),size = 1)+
  xlab("Distance in the river (m)")+
  ylab("Elevation (m)")+
  ggtitle("Eure river")
G2 <- ggplot(data = table_rhine[seq(1,nrow(table_rhine),length.out = 1000),])+
  geom_point(aes(x = D,y = Z),size = 1)+
  xlab("Distance in the river (m)")+
  ylab("Elevation (m)")+
  ggtitle("Rhine river")

GG <- (G1|G2)+plot_layout(guides = "collect")
ggsave(filename = "topology.pdf",plot = GG,device = "pdf",width = 8,height = 4,units = "in")


# Production of the results for a single river

load(file = "fullsection_w_constraint.Rdata")
pdf(file = paste(unique(table_points$name)[i],"pdf",sep = "."),width = 6,height = 3)
L[[i]][[11]]+
  ggtitle(unique(table_points$name)[i])+
  theme(text = element_text(size = 7))
dev.off()
