#  Extract ROMS bottom temperature from the weekly (level-2) ROMS data and clip those data to a shapefile.
#  In this case, we clip to the ESR Bering Sea regions.
#  The final output shows a figure that merges the SST (web service) and BT (ROMS)


library(RNetCDF)
library(tidync)
require(tidyverse)
require(lubridate)
library(sf)

#  Access netcdfs from THREDDS server (https://data.pmel.noaa.gov/aclim/thredds/catalog.html)

#  We are interested in bottom temperature for the Bering Sea ecosystem regions, which are much smaller than the ROMS extent. 
#  Read in our shapefile (which we'll use again later) and extract the coordinates of the bounding box. 
esr_shp <- st_read('Data/Shapefiles/Alaska_Marine_Management_Areas.gdb',layer="Alaska_Marine_Areas_dd") %>% 
  filter(Ecosystem_Subarea%in%c("Northern Bering Sea","Southeastern Bering Sea"))

lkp <- esr_shp %>% 
  st_bbox()%>% 
  st_as_sfc %>% 
  st_cast("MULTIPOINT") %>% 
  st_coordinates() %>% 
  data.frame %>% 
  summarise(maxlat=max(Y),
            minlat=min(Y),
            maxlon=max(X),
            minlon=min(X))

#  Extract bathymetry for the negative longitude regions of Alaskan waters
#  The ROMS file includes bathymetry but I can't figure out how to extract it for these purposes so I use marmap instead.
r.ak <- marmap::as.raster(getNOAA.bathy(lon1=lkp$minlon,lon2=lkp$maxlon,lat1=lkp$minlat,lat2=lkp$maxlat, resolution=1))

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
           lat_rho<=lkp$maxlat) %>% 
  mutate(depth=round(raster::extract(r.ak,cbind(lon_rho,lat_rho),method="bilinear"),0))


# Right now it's just a rectangular grid.
grid_lkp %>% 
  ggplot(aes(lon_rho,lat_rho)) + 
  geom_point()


#  Now load the ROMS data file (note that this does not extract the data so it is quick).
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


#  In case you removed the shapefile that was read in previously, reload it here. 
# esr_shp <- st_read('Data/Shapefiles/Alaska_Marine_Management_Areas.gdb',layer="Alaska_Marine_Areas_dd") %>% 
#   filter(Ecosystem_Subarea%in%c("Northern Bering Sea","Southeastern Bering Sea"))


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
  dplyr::select(Ecosystem_Subarea,lat_rho,lon_rho,xi_rho,eta_rho,BSIERP_ID,depth)

#  Make sure it looks alright
esr_pts %>% 
  ggplot(aes(lon_rho,lat_rho)) + 
  geom_point()


#  Join and save!
readRDS("ROMS_bottom_temp_EBS_since_1985.RDS") %>% 
  inner_join(esr_pts) %>% 
  mutate(date=as_date(as_datetime(ocean_time,origin="1900-01-01 00:00:00", tz = "UTC"))) %>% 
  saveRDS("ROMS_bottom_temp_since_1985_merged_ESR.RDS")




#----------------------------------------------------------------------------------------------
#  Now generate a weekly temperature plot that shows SST and BT for the NBS and SEBS.

#  Read in the data from both the bottom temperature and the surface temperature and plot them together. 
data <- readRDS("ROMS_bottom_temp_since_1985_merged_ESR.RDS") %>% 
  rename(Ecosystem_sub=Ecosystem_Subarea) %>% 
  filter(between(depth,-200,-10)) %>% 
  group_by(date,Ecosystem_sub) %>% 
  summarise(meanbt=mean(temp)) %>% 
  inner_join(httr::content(httr::GET('https://apex.psmfc.org/akfin/data_marts/akmp/ecosystem_sub_crw_avg_sst?ecosystem_sub=Southeastern%20Bering%20Sea,Northern%20Bering%20Sea&start_date=19850101&end_date=20211231'), type = "application/json") %>% 
               bind_rows %>% 
               mutate(date=as_date(READ_DATE)) %>% 
               data.frame %>% 
               dplyr::select(date,meansst=MEANSST,Ecosystem_sub=ECOSYSTEM_SUB))


#  Plot sst and bt
data %>%
  gather(index,temperature,-c(date,Ecosystem_sub)) %>% 
  ggplot(aes(date,temperature,linetype=index)) + 
  geom_line() + 
  facet_wrap(~Ecosystem_sub,ncol=1) + 
  theme_bw()
  

  # Plot the difference between SST and BT  
  data %>%
    mutate(diff=meansst-meanbt) %>% 
    ggplot(aes(date,diff)) + 
    geom_line() + 
    geom_smooth() +
    facet_wrap(~Ecosystem_sub,ncol=1) + 
    theme_bw()
  
temp <-  data %>%
    mutate(diff=meansst-meanbt,
           doy=yday(date),
           year=year(date),
           month=month(date),
           day=day(date),
           newdate=as.Date(ifelse(month>=9,as.character(as.Date(paste("1999",month,day,sep="-"),format="%Y-%m-%d")),
                                  as.character(as.Date(paste("2000",month,day,sep="-"),format="%Y-%m-%d"))),format("%Y-%m-%d")),
           year2=ifelse(month>=9,year+1,year),
           doy2=(as.vector(newdate-as.Date("1999-09-01"))+1)) %>% 
    arrange(date) %>% 
  filter(year2>=1986)

OceansBlue2='#0055A4' # NOAA dark blue

# What was the earliest date during which the surface temperatures were less than the bottom temperatures?
png("ROMS_ESR/Cold_surface_timing.png",width=7,height=5,units="in",res=300)
temp %>% 
  group_by(Ecosystem_sub,year2) %>% 
  summarise(crossdate=min(doy2[diff<0])) %>% 
  mutate(crossdate=ifelse(is.infinite(crossdate),NA,crossdate)) %>% 
  ggplot(aes(year2,crossdate)) + 
  geom_bar(stat="identity",fill=OceansBlue2) + 
  facet_wrap(~Ecosystem_sub,ncol=1) + 
  xlab("Year") + 
  ylab("Date") + 
  theme(strip.text = element_text(size=10,color="white",family="sans",face="bold"),
        strip.background = element_rect(fill=OceansBlue2),
        axis.title = element_text(size=10,family="sans"),
        axis.text = element_text(size=10,family="sans"),
        panel.background = element_blank(),
        panel.border=element_rect(fill=NA,colour="black",size=0.5),
        plot.margin=unit(c(0.65,0,0.65,0),"cm"))
dev.off()

#  How many days had SST less than bottom temperatures?
png("ROMS_ESR/Cold_days.png",width=7,height=5,units="in",res=300)
temp %>% 
  group_by(Ecosystem_sub,year2) %>% 
  summarise(colddays=length(newdate[diff<0])) %>% 
  ggplot(aes(year2,colddays)) + 
  geom_bar(stat="identity",fill=OceansBlue2) + 
  facet_wrap(~Ecosystem_sub,ncol=1) + 
  xlab("Year") + 
  ylab("Number of weeks") + 
  theme(strip.text = element_text(size=10,color="white",family="sans",face="bold"),
        strip.background = element_rect(fill=OceansBlue2),
        axis.title = element_text(size=10,family="sans"),
        axis.text = element_text(size=10,family="sans"),
        panel.background = element_blank(),
        panel.border=element_rect(fill=NA,colour="black",size=0.5),
        plot.margin=unit(c(0.65,0,0.65,0),"cm"))
dev.off()

#  Explore the data based on depth bins, 0-950 and 51-200 within each Bering ESR region
roms <- readRDS("ROMS_bottom_temp_since_1985_merged_ESR.RDS") %>% 
  rename(Ecosystem_sub=Ecosystem_Subarea) %>% 
  filter(between(depth,-200,0)) %>% 
  mutate(eco2=ifelse(depth<(-50),paste0(Ecosystem_sub," (51-200m)"),paste0(Ecosystem_sub," (0-50m)"))) %>% 
  group_by(date,eco2) %>% 
  summarise(meanbt=mean(temp)) %>% 
  ungroup

#  The ESR_sst_depthbinds.RDS file was created using the "SQL_query_sst_depthbins.R" file in the ROMS_ESR folder
#  THis was done via the Oracle database.
mergedat <- roms %>% 
  inner_join(readRDS("ROMS_ESR/ESR_sst_depthbins.RDS") %>% 
               mutate(date=as_date(READ_DATE)) %>% 
               dplyr::select(sst=SST,eco2,date))

#  Spatial illustration of ROMS bottom temperature only by depth stratified regions
png("ROMS_ESR/ROMS_bt_by_depth_region.png",width=7,height=5,units="in",res=300)
mergedat %>% 
  ggplot(aes(date,meanbt),size=0.35) + 
  geom_line() + 
  geom_smooth(size=0.35) + 
  facet_wrap(~eco2,ncol=2) + 
  theme(strip.text = element_text(size=10,color="white",family="sans",face="bold"),
        strip.background = element_rect(fill=OceansBlue2),
        axis.title = element_text(size=10,family="sans"),
        axis.text = element_text(size=10,family="sans"),
        panel.background = element_blank(),
        panel.border=element_rect(fill=NA,colour="black",size=0.5),
        plot.margin=unit(c(0.65,0,0.65,0),"cm")) + 
  ylab("Bottom temperature (C)") + 
  xlab("Date")
dev.off() 

#  Spatial illustration of ROMS bottom temperature only by depth stratified regions
png("ROMS_ESR/Diss_sst_bt_depth_region.png",width=7,height=5,units="in",res=300)
mergedat %>% 
  mutate(diff=sst-meanbt) %>% 
  ggplot(aes(date,diff)) + 
  geom_line(size=0.35) + 
  geom_smooth(size=0.35) + 
  facet_wrap(~eco2,ncol=2) + 
  theme(strip.text = element_text(size=10,color="white",family="sans",face="bold"),
        strip.background = element_rect(fill=OceansBlue2),
        axis.title = element_text(size=10,family="sans"),
        axis.text = element_text(size=10,family="sans"),
        panel.background = element_blank(),
        panel.border=element_rect(fill=NA,colour="black",size=0.5),
        plot.margin=unit(c(0.65,0,0.65,0),"cm")) + 
  ylab("Surface temperature - bottom temperature (C)") +
  xlab("Date")
dev.off()