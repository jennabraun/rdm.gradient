#get site climate data from worldclim

library(raster)
library(sp)

sites <- read.csv("clean_data/sites.csv")
site <- sites$ID
coords <- data.frame(x=sites$Longitude,y=sites$Latitude)

points <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
#panoche
r <- getData("worldclim",var="bio",res=0.5, lon = -120.8120, lat = 36.70654)
r <- r[[c(1,12, 5)]]
names(r) <- c("Temp","Prec", "Max")
values <- extract(r,points)
pan <- cbind.data.frame(coordinates(points),values)
pan

r <- getData("worldclim",var="bio",res=0.5, lon = -116.8349	, lat = 35.09405)
r <- r[[c(1,12,5)]]
names(r) <- c("Temp","Prec", "Max")
values <- extract(r,points)
moj <- cbind.data.frame(coordinates(points),values)
moj

#put together
pan <- filter(pan, Temp != "NA")
moj <- filter(moj, Temp != "NA")
sites <- rbind(pan, moj)
sites$site <- site
sites <- mutate(sites, Temp = Temp/10)
sites <- mutate(sites, Max = Max/10)
sites <- mutate(sites, arid = Prec/(Temp+10))
write.csv(sites, "clean_data/sites_worldclim.csv")
