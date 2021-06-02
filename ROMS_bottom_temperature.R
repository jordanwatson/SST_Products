library(RNetCDF)
library(tidync)
require(tidyverse)
require(lubridate)
library(sf)

#Redone with ROMS data

#  Access netcdfs from THREDDS server (https://data.pmel.noaa.gov/aclim/thredds/catalog.html)

#  We are interested in bottom temperature for the Bering Sea ecosystem regions, which are much smaller than the ROMS extent. 
#  Read in our lookup table for the Bering Sea to find some spatial bounds.
lkp <- readRDS("Data/crwsst/crwsst_spatial_lookup_table.RDS") %>% 
  filter(Ecosystem=="Eastern Bering Sea") %>% 
  dplyr::select(longitude,latitude) %>% 
  summarise(maxlat=max(latitude),
            minlat=min(latitude),
            maxlon=max(longitude),
            minlon=min(longitude))

#  We need to access the file containing the extended spatial grids, which contains different transformations of the spatial coordinates.
grid<-tidync("https://data.pmel.noaa.gov/aclim/thredds/dodsC/extended_grid/Bering10K_extended_grid.nc")

#  One of the native grids contains the columns we need to match the bottom temperature data to latitudes and longitudes.
grid_lkp <- grid %>% 
  activate("D8,D2") %>% 
  hyper_tibble() %>% 
  dplyr::select(lat_rho,lon_rho,xi_rho,eta_rho) %>% 
  filter(lon_rho>180) %>% 
  mutate(lon_rho=lon_rho-360,
         latitude=lat_rho,      # Make a duplicate column for the point-in-polygon operation
         longitude=lon_rho) %>%       # Make a duplicate column for the point-in-polygon operation
  filter(lon_rho>=lkp$minlon & 
           lon_rho<=lkp$maxlon & 
           lat_rho>=lkp$minlat & 
           lat_rho<=lkp$maxlat)

# Right now it's just a rectangular grid.
grid_lkp %>% 
  ggplot(aes(lon_rho,lat_rho)) + 
  geom_point()

#  Now load load the ROMS data file (note that this does not extract the data so it is quick).
ROMS <- tidync("https://data.pmel.noaa.gov/aclim/thredds/dodsC/Level2/B10K-K20_CORECFS_bottom5m.nc")

# We can create a data frame of dates. To match our SST data, we only need bottom temperatures since 1985-01-01.
# So we identify which dates to filter our from our subsequent query. 
datevec <- ROMS %>% 
  activate("D4") %>% 
  hyper_tibble() %>% 
  mutate(date=as_date(as_datetime(ocean_time,origin="1900-01-01 00:00:00", tz = "UTC")),
         date_index=1:n()) %>% 
  filter(date>=as.Date("1985-01-01"))

#  This is a beast. Extract all the temperature across the lookup grid since 1985. May need to come back and do this by year.
#  Ends up being about 27 million rows. Save directly to RDS.
ROMS %>% 
  activate("D8,D1,D4") %>% 
  hyper_filter(ocean_time=ocean_time>=ocean_time[min(datevec$date_index)],    # Filter dates based on the datevec index for 1985-01-01
               xi_rho=xi_rho>=min(grid_lkp$xi_rho) & xi_rho<=max(grid_lkp$xi_rho), # Filter xi_rho based on the spatial lookup
               eta_rho=eta_rho>=min(grid_lkp$eta_rho) & eta_rho<=max(grid_lkp$eta_rho)) %>%  # Filter eta_rho based on the spatial lookup
  hyper_tibble(select_var="temp") %>% 
  saveRDS("ROMS_bottom_temp_EBS_since_1985.RDS")# Only extract the temperature variable



#  Read in the ESR shapefile and subset for the Bering areas
esr_shp <- st_read('Data/Shapefiles/Alaska_Marine_Management_Areas.gdb',layer="Alaska_Marine_Areas_dd") %>% 
  filter(Ecosystem_Subarea%in%c("Northern Bering Sea","Southeastern Bering Sea"))

#  Plot it.
#ggplot() + 
#  geom_sf(data=esr_shp)

# Convert the lookup grid to a sf object with CRS and then transform to that of the shapefile and 
# perform the point-in-polygon operation, removing unmatched coordinates.
esr_pts = st_join(
  st_as_sf(grid_lkp, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")  %>% # Use the duplicated lat/lon columns for matching to avoid rounding issues.
  st_transform(st_crs(esr_shp)$proj4string),
  esr_shp) %>% 
  filter(!is.na(Ecosystem_Subarea)) %>% 
  data.frame %>% # THis will drop the latitude and longitude point geometry column
  dplyr::select(Ecosystem_Subarea,lat_rho,lon_rho,xi_rho,eta_rho,BSIERP_ID)

#  Make sure it looks alright
esr_pts %>% 
  ggplot(aes(lon_rho,lat_rho)) + 
  geom_point()


#  Join and save!
readRDS("ROMS_bottom_temp_EBS_since_1985.RDS") %>% 
  inner_join(esr_pts) %>% 
  mutate(date=as_date(as_datetime(ocean_time,origin="1900-01-01 00:00:00", tz = "UTC"))) %>% 
  saveRDS("ROMS_bottom_temp_since_1985_merged_ESR.RDS")


data <- readRDS("ROMS_bottom_temp_since_1985_merged_ESR.RDS") %>% 
  rename(Ecosystem_sub=Ecosystem_Subarea) %>% 
  group_by(date,Ecosystem_sub) %>% 
  summarise(meanbt=mean(temp)) %>% 
  inner_join(httr::content(httr::GET('https://apex.psmfc.org/akfin/data_marts/akmp/ecosystem_sub_crw_avg_sst?ecosystem_sub=Southeastern%20Bering%20Sea,Northern%20Bering%20Sea&start_date=19850101&end_date=20211231'), type = "application/json") %>% 
               bind_rows %>% 
               mutate(date=as_date(READ_DATE)) %>% 
               data.frame %>% 
               dplyr::select(date,meansst=MEANSST,Ecosystem_sub=ECOSYSTEM_SUB)) %>% 
  gather(index,temperature,-c(date,Ecosystem_sub))


data %>% 
  ggplot(aes(date,temperature,linetype=index)) + 
  geom_line() + 
  facet_wrap(~Ecosystem_sub,ncol=1) + 
  theme_bw()
