library(tidync)
library(tidyverse)
library(lubridate)
library(RNetCDF)
library(sf)
library(marmap)

#  Extract the grid data from one of the netcdf files (this creates an array)
#  Note that hyper_array here pulls the full lat-lon grid, which is important since no one day will have 
#  the full grid because of missing chlorophyll data. 
df1 <- attr(tidync("Data/chlorophyll/modis_8day/chl_annual_2005.nc") %>% 
              hyper_filter(time=time==time[1]) %>% 
              hyper_array(),
            "transforms")

#  Extract bathymetry for the negative longitude regions of Alaskan waters
r.ak <- marmap::as.raster(getNOAA.bathy(lon1=-180,lon2=-129,lat1=47.5,lat2=69, resolution=1))

#  Extract the latitudes and negative longitudes and pull the matching bathymetry data.
grid1 <- expand.grid(longitude=df1$longitude$longitude,
                     latitude=df1$latitude$latitude) %>% 
  mutate(depth=round(raster::extract(r.ak,cbind(longitude,latitude),method="bilinear"),0))

#  Now repeat for the positive longitude grids
dfpos <- attr(tidync("Data/chlorophyll/modis_8day/chl_annual_positive2005.nc") %>% 
              hyper_filter(time=time==time[1]) %>% 
              hyper_array(),
            "transforms")

#  Extract bathymetry for the negative longitude regions of Alaskan waters
r.ak <- marmap::as.raster(getNOAA.bathy(lon1=166,lon2=180,lat1=46,lat2=61, resolution=1))

#  Extract the lat/long data from the array (above) and truncate the decimals to four places
gridpos <- expand.grid(longitude=dfpos$longitude$longitude,
                     latitude=dfpos$latitude$latitude) %>% 
  mutate(depth=round(raster::extract(r.ak,cbind(longitude,latitude),method="bilinear"),0))

#  Merge the positive and negative grids. 
mygrid <- grid1 %>% 
  bind_rows(gridpos) %>% 
  data.frame %>% 
  distinct() %>% 
  mutate(lon=longitude, # Add a separate lat/lon field for the spatial transformation to avoid unintended transformation of these primary keys.
         lat=latitude)
  
#  Clean up.
rm(grid1,gridpos);gc()

#  Read in the shapefile with all of our polygons
esr_shp <- st_read('Data/Shapefiles/Alaska_Marine_Management_Areas.gdb',layer="Alaska_Marine_Areas_dd")

# Convert the lookup grid to a sf object with CRS and then transform to that of the shapefile and 
# perform the point-in-polygon operation, removing unmatched coordinates.
esr_pts = st_join(
  st_as_sf(mygrid, coords = c("lon", "lat"), crs = 4326, agr = "constant")  %>% # Use the duplicated lat/lon columns for matching to avoid rounding issues.
    st_transform(st_crs(esr_shp)$proj4string),
  esr_shp) %>% 
  #data.frame %>% # THis will drop the latitude and longitude point geometry column
  # dplyr::mutate(longitude = sf::st_coordinates(.)[,1],  #  Extract the lat-lon columns from spatial object
  #               latitude = sf::st_coordinates(.)[,2]) %>% 
  data.frame %>% 
  dplyr::select(-geometry)

#  We don't want duplicated rows for lat-long coordinates so we need to join coordinates that are matched to individual strata (e.g., nmfs areas versus ESR areas)
lkp.grid <- esr_pts %>% # select the points that fall within ESR polygons
  dplyr::select(longitude,latitude,Ecosystem=Ecosystem_Area,Ecosystem_sub=Ecosystem_Subarea) %>% 
  filter(!is.na(Ecosystem_sub)) %>% 
  full_join(esr_pts %>% # select the points that fall within nmfs area polygons
              dplyr::select(longitude,latitude,nmfsarea=NMFS_REP_AREA) %>% 
              filter(!is.na(nmfsarea))) %>% 
  full_join(esr_pts %>% # select the points that fall within stat area polygons
              dplyr::select(longitude,latitude,stat_area=STAT_AREA,statefed=WATERS_COD) %>% 
              filter(!is.na(stat_area))) %>% 
  full_join(esr_pts %>% # select the points that fall within BSIERP polygons
              dplyr::select(longitude,latitude,bsierp_name=BSIERP_Region_Name,bsierp_id=BSIERP_ID) %>% 
              filter(!is.na(bsierp_name))) %>% 
  left_join(mygrid %>% dplyr::select(-c(lon,lat))) # Now get the depth data

#  Plot it to make sure if worked.
# lkp.grid %>%
#   mutate(longitude=ifelse(longitude<0,longitude+360,longitude)) %>%
#   ggplot() +
#   geom_point(aes(longitude,latitude))

# lkp.grid %>%
#   saveRDS("Data/chlorophyll/Updated_modis_lookup_table.RDS")

lkp.grid <- readRDS("Data/chlorophyll/Updated_modis_lookup_table.RDS")

#  Do a test case
# test <- tidync("Data/chlorophyll/modis_8day/chl_annual_2005.nc") %>%
#   hyper_filter(time=time==time[10]) %>%
#   hyper_tibble() %>%
#   right_join(lkp.grid %>% filter(longitude<0)) %>%
#   bind_rows(tidync("Data/chlorophyll/modis_8day/chl_annual_positive2005.nc") %>%
#               hyper_filter(time=time==time[10]) %>%
#               hyper_tibble() %>%
#               right_join(lkp.grid %>% filter(longitude>0)))
# 
# test %>%
#   mutate(longitude=ifelse(longitude<0,longitude+360,longitude)) %>%
#   ggplot(aes(longitude,latitude,color=chlorophyll)) +
#   geom_point()

#  Load the annual files, perform the spatial join and clip to our strata boundaries.
for(i in 2003:2020){
  tidync(paste0("Data/chlorophyll/modis_8day/chl_annual_",i,".nc")) %>%
    #hyper_filter(time=time==time[10]) %>%
    hyper_tibble() %>%
    right_join(lkp.grid %>% filter(longitude<0)) %>%
    bind_rows(tidync(paste0("Data/chlorophyll/modis_8day/chl_annual_positive",i,".nc")) %>%
                #hyper_filter(time=time==time[10]) %>%
                hyper_tibble() %>%
                right_join(lkp.grid %>% filter(longitude>0))) %>% 
    saveRDS(paste0("Data/chlorophyll/modis_8day/merged/chl_8day_merged_",i,".RDS"))
}

#  Deal with the unfinished current year.
tidync(paste0("Data/chlorophyll/modis_8day/chl_2021_through_07272021.nc")) %>%
  #hyper_filter(time=time==time[10]) %>%
  hyper_tibble() %>%
  right_join(lkp.grid %>% filter(longitude<0)) %>%
  bind_rows(tidync(paste0("Data/chlorophyll/modis_8day/chl_positive_2021_through_07272021.nc")) %>%
              #hyper_filter(time=time==time[10]) %>%
              hyper_tibble() %>%
              right_join(lkp.grid %>% filter(longitude>0))) %>% 
  saveRDS(paste0("Data/chlorophyll/modis_8day/merged/chl_8day_merged_2021_through_07272021.RDS"))

# Make a big ol' merged file.
newdat <- lapply(list.files("Data/chlorophyll/modis_8day/merged/"),
       function(x) readRDS(paste0("Data/chlorophyll/modis_8day/merged/",x)) %>% 
         mutate(date=as_date(as_datetime(time))) %>% 
         dplyr::select(-time)) %>% 
  bind_rows()

# Create a merged MODIS chlorophyll dataset.
newdat %>% 
  saveRDS("Data/chlorophyll/modis_8day/merged/merged_8day_2003_2021.RDS")

# Save smaller versions, including a csv version for Noel Pelland.
newdat %>% 
  filter(Ecosystem=="Aleutian Islands") %>% 
  write.csv("Data/chlorophyll/modis_8day/merged/merged_8day_2003_2021_Aleutians.csv")

newdat %>% 
  filter(Ecosystem=="Aleutian Islands") %>% 
  saveRDS("Data/chlorophyll/modis_8day/merged/merged_8day_2003_2021_Aleutians.RDS")

newdat %>% 
  filter(Ecosystem=="Eastern Bering Sea") %>% 
  saveRDS("Data/chlorophyll/modis_8day/merged/merged_8day_2003_2021_EBS.RDS")

newdat %>% 
  filter(Ecosystem=="Gulf of Alaska") %>% 
  saveRDS("Data/chlorophyll/modis_8day/merged/merged_8day_2003_2021_GOA.RDS")
