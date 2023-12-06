plotting = FALSE

## Reading in the strata system in the mackerel world and add it on top of what we use 
## Note that area 8 (norway coast) is missing so be careful about the numbering
## And other manually defined strata

files <- list.files("data/")
files <- files[grep("stratum", files)]
files <- files[order(as.numeric(gsub("_", "", substr(files, start=16, stop=17))))]

strata <- list()
All_strata <- c()
for (fl in seq_along(files)){
  bounds <- read.csv(paste0(getwd(), "/data/", files[fl]), sep=" ")
  ID = substr(files[fl], start=16, stop=17)
  ID = gsub("_", "", ID)
  bounds$ID <- as.numeric(ID)
  bounds <- bounds[,-3]
  All_strata <- rbind(All_strata, bounds)
  my_df_sf <- st_as_sf(x=bounds, coords = c("lon", "lat")) %>% st_set_crs(4326) %>% group_by(ID) %>% summarize(geometry = st_combine(geometry)) %>%
    st_cast("POLYGON")
  strata[[fl]] <- my_df_sf
}

my_df_sf <- st_as_sf(x=All_strata, coords = c("lon", "lat")) %>% st_set_crs(4326) %>% group_by(ID) %>% 
  summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
if (plotting == TRUE) ggplot(Atlantic) + geom_sf() + theme_bw() + geom_sf(data = my_df_sf, aes(fill=ID), size = 0.2) + 
  geom_sf_label(data = my_df_sf, aes(label = ID)) + scale_fill_viridis_c()

 all_strata <- st_read("Shapefiles/all_strata.shp")
if (plotting == TRUE) ggplot(Atlantic) + geom_sf() + theme_bw() + geom_sf(data = all_strata, size = 0.2, fill="lightblue1") 

## Load in the triangle area in between area 2 and 3
Triangle_area <- read.table(file="data/yellow triangle.txt",header=T)
Triangle_area$ID <- "Triangle"
Triangle_area_sf <- st_as_sf(x=Triangle_area, coords = c("lon", "lat")) %>% st_set_crs(4326) %>% 
  summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
if (plotting == TRUE) ggplot(Atlantic) + geom_sf() + theme_bw() + geom_sf(data = Triangle_area_sf, size = 0.2, fill="lightblue1") 

## Also modify the area 9 to include only north of 74 and 75
# Area9_74 <- read.table(file="data/Area9_74.txt",header=T)
Area9_74 <- read.table(file="data/Area9_74_broad.txt",header=T)
Area9_74$ID <- "9_74"
Area9_74_sf <- st_as_sf(x=Area9_74, coords = c("lon", "lat")) %>% st_set_crs(4326) %>% group_by(ID) %>% summarize(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
if (plotting == TRUE) ggplot(Atlantic) + geom_sf() + theme_bw() + geom_sf(data = Area9_74_sf, size = 0.2, fill="lightblue1") 
# Area9_75 <- read.table(file="data/Area9_75.txt",header=T)
Area9_75 <- read.table(file="data/Area9_75_broad.txt",header=T)
Area9_75$ID <- "9_75"
Area9_75_sf <- st_as_sf(x=Area9_75, coords = c("lon", "lat")) %>% st_set_crs(4326) %>% group_by(ID) %>% summarize(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
if (plotting == TRUE) ggplot(Atlantic) + geom_sf() + theme_bw() + geom_sf(data = Area9_75_sf, size = 0.2, fill="lightblue1") 

# bla <- st_intersection(Area9_75, Area9_74)
if (plotting == TRUE) ggplot(Atlantic) + geom_sf() + theme_bw() + geom_sf(data = bla, size = 0.2, fill="lightblue1") 
# head(sf_to_df(bla))

## Create a fake land masses south of 60N
South_end <- data.frame(lon=c(-65,-15,-15,-65,-65), lat=c(60,60,46,46,60))
South_end$ID <- "South_end"
South_end_sf <- st_as_sf(x=South_end, coords = c("lon", "lat")) %>% st_set_crs(4326) %>% 
  summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
if (plotting == TRUE) ggplot(Atlantic) + geom_sf() + theme_bw() + geom_sf(data = South_end_sf, size = 0.2, fill="lightblue1") 

## Create a fake land masses west of area 10
West_10 <- data.frame(lon=c(-65,-65,-29,-29,-40,-40.5,-40.79,-65), lat=c(62,66.2,66.2,66.2,63.7,62.5,62,62))
West_10$ID <- "West_10"
West_10_sf <- st_as_sf(x=West_10, coords = c("lon", "lat")) %>% st_set_crs(4326) %>% 
  summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
if (plotting == TRUE) ggplot(Atlantic) + geom_sf() + theme_bw() + geom_sf(data = West_10_sf, size = 0.2, fill="lightblue1") 

## Create a fake land masses west of area 11
West_11 <- data.frame(lon=c(-65,-40.79,-40.79,-65,-65), lat=c(62,62,60,60,62))
West_11$ID <- "West_11"
West_11_sf <- st_as_sf(x=West_11, coords = c("lon", "lat")) %>% st_set_crs(4326) %>% 
  summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
if (plotting == TRUE) ggplot(Atlantic) + geom_sf() + theme_bw() + geom_sf(data = West_11_sf, size = 0.2, fill="lightblue1") 

## Create a fake land masses west of greenland
West_Greenland <- data.frame(lon=c(-65,-43,-43,-65,-65), lat=c(80,80,60,60,80))
West_Greenland$ID <- "West_Greenland"
West_Greenland_sf <- st_as_sf(x=West_Greenland, coords = c("lon", "lat")) %>% st_set_crs(4326) %>% 
  summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
if (plotting == TRUE) ggplot(Atlantic) + geom_sf() + theme_bw() + geom_sf(data = West_Greenland_sf, size = 0.2, fill="lightblue1") 


## Update All_bounds to include the triange area
## Update All_bounds to include the triange area
All_bounds <- rbind(All_strata, Triangle_area, Area9_74, Area9_75, South_end, West_10, West_11, West_Greenland)
my_df_sf <- st_as_sf(x=All_bounds, coords = c("lon", "lat")) %>% st_set_crs(4326) %>% group_by(ID) %>% 
  summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
if (plotting == TRUE) ggplot(Atlantic) + geom_sf() + theme_bw() + geom_sf(data = my_df_sf, aes(fill=as.factor(ID)), size = 0.2) + 
  # geom_sf_label(data = my_df_sf, aes(label = ID)) + scale_fill_viridis_d()

if (plotting == TRUE) {
 p1 <- ggplot(Atlantic) + geom_sf() + theme_bw() + geom_sf(data = all_strata, size = 0.2, fill="lightblue1") +
   geom_sf(data = my_df_sf, aes(col=as.factor(ID)), fill=NA, size = 0.2) + 
   geom_sf_label(data = my_df_sf, aes(label = ID)) + scale_color_viridis_d(name="Strata ID") +
   xlab("Longitude") + ylab("Latitude") + 
   theme(axis.title = element_text(size=15), axis.text = element_text(size=13))
 ggsave(p1, filename = paste0(getwd(), "/plots/Strata.pdf"), dpi ="retina", width = 10, height = 8, device = "pdf")
}




asd <- my_df_sf %>% filter(ID %in% c(1:7,9:12))
asd <- asd %>% mutate(col = ifelse(ID %in% c(1,2,3,5,6,7,10,11), "Permanent", "Dynamic"))

Fig1 <- ggplot(Atlantic) + geom_sf() + theme_bw() + 
  geom_sf(data = asd, aes(fill=col), col=1, size = 0.2) + 
  geom_sf_label(data = asd, aes(label = ID)) + 
  scale_fill_viridis_d(name="Strata type") +
  xlab("Longitude") + ylab("Latitude") + 
  theme(axis.title = element_text(size=15), axis.text = element_text(size=13))
ggsave(Fig1, filename = paste0(getwd(), "/MS/Figs/Strata.png"), dpi ="retina", width = 8, height = 6, device = "png")

