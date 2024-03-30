library(sf)
library(terra)
library(grid)
library(rgbif)
library(maptools)
library(spData)
library(deeptime)
library(patchwork)
library(viridis)
#library(ggtern)
library(stringr)
#devtools::install_github("r-barnes/dggridR", vignette=TRUE)
library(dggridR)
library(spdep)
library(moments)
library(quantreg)
library(ggforce)
library(relaimpo)
library(MuMIn)
library(parallel)
library(tidyverse)

#from https://github.com/skgallagher/EpiCompare/blob/master/R/aaa.R#L61
imports_hidden_from <- function(pkg, name){
  pkg <- as.character(substitute(pkg))
  name <- as.character(substitute(name))
  get(name, envir = asNamespace(pkg), inherits = FALSE)
}

#Fix the print method for ggtern
print.ggplot <- function(x, newpage = is.null(vp), vp = NULL, ...){
  if(inherits(x$coordinates, "CoordTern")){
    #ggtern:::print.ggplot(x, newpage = newpage, vp = vp, ...)
    print.ggplot <- imports_hidden_from("ggtern", "print.ggplot")
    print.ggplot(x, newpage = newpage, vp = vp, ...)
  } else {
    #ggplot2:::print.ggplot(x, newpage = newpage, vp = vp, ...)
    print.ggplot <- imports_hidden_from("ggplot2", "print.ggplot")
    print.ggplot(x, newpage = newpage, vp = vp, ...)
  }
}

colsea = "#ebf7ff"
colland = "#ededed"

#try different cut-offs of area of overlap of grid cells and ranges
#try different grid cell sizes
#remove cosmopolitan species? top 10%?

#IUCN ranges ####
#terrestrial mammals using IUCN ranges
#download from here: https://www.iucnredlist.org/resources/spatial-data-download
#then unzip all of the files into a folder named "MAMMALS_TERRESTRIAL_ONLY"
mammal_IUCN <- st_read(dsn = "MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")

#a map of the world
data(wrld_simpl)
wrld_sf <- st_as_sf(wrld_simpl)

#get IUCN species, including synonyms
# "search summary" data retrieved from:
# https://www.iucnredlist.org/search?dl=true&permalink=6da36cdc-5f19-46e7-b980-aeb06b2b6208
mammal_IUCN_synonyms <- subset(read.csv("synonyms.csv", stringsAsFactors = FALSE), infraType != "subspecies") %>%
  unite(species, genusName, speciesName, sep = " ") %>% dplyr::select(IUCN_species = scientificName, species, speciesAuthor)
mammal_IUCN_synonyms$species <- trimws(gsub("<i>|</i>", "", mammal_IUCN_synonyms$species))
mammal_IUCN_synonyms$speciesAuthor <- mammal_IUCN_synonyms$speciesAuthor %>% str_match("[0-9]+") %>% unlist %>% as.numeric
mammal_IUCN_synonyms <- mammal_IUCN_synonyms[order(mammal_IUCN_synonyms$species, -mammal_IUCN_synonyms$speciesAuthor),]
mammal_IUCN_synonyms <- mammal_IUCN_synonyms[!duplicated(mammal_IUCN_synonyms$species),1:2]

mammal_IUCN_taxonomy <- read.csv("taxonomy.csv", stringsAsFactors = FALSE) %>%
  dplyr::select(IUCN_species = scientificName)
mammal_IUCN_taxonomy$species <- mammal_IUCN_taxonomy$IUCN_species
mammal_IUCN_synonyms <- unique(rbind(mammal_IUCN_synonyms[!mammal_IUCN_synonyms$species %in% mammal_IUCN_taxonomy$species,], mammal_IUCN_taxonomy))

#Trait data####
#get trait data from EltonTraits https://doi.org/10.6084/m9.figshare.c.3306933.v1
mammal_traits <- read.table("EltonTraits.txt", sep = "\t", header = T, stringsAsFactors = FALSE) #Wilman et al 2014
mammal_traits$Plant.Per <- rowSums(mammal_traits[c("Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO")])
mammal_traits$Trophic.Level <- cut(mammal_traits$Plant.Per, c(-1,5,95,101), c("Carnivore", "Omnivore", "Herbivore"))
mammal_traits$BodyMass.log10 <- log10(mammal_traits$BodyMass.Value)
# small: 1g - 1000g (1kg)
# medium: 1kg - 10kg
# large: 10kg+
mammal_traits$size_cat <- cut(mammal_traits$BodyMass.log10, c(0,3,4,10), labels = c("small", "medium", "large"), ordered_result = TRUE)
# large: 10kg - 100kg
# xlarge: 100kg - 1000kg (bovids, big deer, big cats, bears, big boars, giraffes, camels, gorillas)
# xxlarge: 1000kg+ (elephants and rhinos)
mammal_traits$size_cat2 <- cut(mammal_traits$BodyMass.log10, c(0,3,4,5,6,10),
                               labels = c("small", "medium", "large", "xlarge", "xxlarge"), ordered_result = TRUE)

#merge in trait data
mammal_traits_IUCN <- merge(mammal_IUCN_synonyms, mammal_traits, by.x = "IUCN_species", by.y = "Scientific", all.x = TRUE)

#fill in traits based on synonyms if still missing
mammal_traits_IUCN[is.na(mammal_traits_IUCN$Trophic.Level), colnames(mammal_traits)[colnames(mammal_traits) %in% colnames(mammal_traits_IUCN)]] <- mammal_traits[match(mammal_traits_IUCN$species[is.na(mammal_traits_IUCN$Trophic.Level)], mammal_traits$Scientific), colnames(mammal_traits)[colnames(mammal_traits) %in% colnames(mammal_traits_IUCN)]]

#remove duplicates due to synonyms
mammal_traits_IUCN <- mammal_traits_IUCN %>% dplyr::select(-species) %>% unique()
mammal_traits_IUCN <- mammal_traits_IUCN[!is.na(mammal_traits_IUCN$Trophic.Level),]
mammal_traits_IUCN <- mammal_traits_IUCN[order(mammal_traits_IUCN$IUCN_species, mammal_traits_IUCN$Diet.Certainty, -(mammal_traits_IUCN$BodyMass.SpecLevel == 1)),]
mammal_traits_IUCN <- mammal_traits_IUCN[!duplicated(mammal_traits_IUCN$IUCN_species),]

# remove bats
mammal_IUCN_filt <- mammal_IUCN %>%
  filter(order_ != "CHIROPTERA")

#Setup Spatial Cells####
#setup equal area grid cells
#https://cran.microsoft.com/snapshot/2017-12-10/web/packages/dggridR/vignettes/dggridR.html
#resolution based on this: http://www.pnas.org/content/104/33/13384
#and this: https://www.jstor.org/stable/pdf/23064563.pdf
#~250km resolution (~70,000km^2 per grid cell)
dggs <- dgconstruct(res = 6)
grid <- st_wrap_dateline(dgearthgrid(dggs))

#get coordinates of midpoints of grid cells
grid_coord <- st_centroid(grid)

#remove grid cells that aren't over land
sf_use_s2(FALSE) #st_make_valid doesn't work well with s2 right now
grid_overlap <- cbind(grid_coord, over_land = sapply(st_within(grid_coord, st_make_valid(wrld_sf)), function(x) length(x) > 0))
sf_use_s2(TRUE)
grid_filt <- grid[grid_overlap$over_land,]
grid_filt_vect <- vect(grid_filt)

#Get geographic data####
#human footprint index
#https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-2009-human-footprint
#https://www.nature.com/articles/sdata201667
#https://academic.oup.com/bioscience/article/52/10/891/354831
#~1km resolution
footprint.2009 <- app(rast("wildareas-v3-2009-human-footprint.tif"), fun=function(x) { x[x>50 | x<0] <- NA; return(x) })
footprint.1993 <- app(rast("wildareas-v3-1993-human-footprint.tif"), fun=function(x) { x[x>50 | x<0] <- NA; return(x) })

footprint.2009.means <- terra::extract(footprint.2009, project(grid_filt_vect, crs(footprint.2009)), fun = mean, na.rm = T)
colnames(footprint.2009.means) <- c("seqnum", "footprint.2009")
footprint.1993.means <- terra::extract(footprint.1993, project(grid_filt_vect, crs(footprint.2009)), fun = mean, na.rm = TRUE)
colnames(footprint.1993.means) <- c("seqnum", "footprint.1993")
footprint.2009.means$seqnum <- footprint.1993.means$seqnum <- grid_filt$seqnum

#deforestation data
#https://millenniumassessment.org/documents/document.356.aspx.pdf
#https://sedac.ciesin.columbia.edu/data/set/ma-rapid-land-cover-change
#.1 degree = ~10km resolution
deforestation <- rast("def_hot_for/w001001.adf")
#make other land pixels 0
wrld_rast <- subst(rasterize(wrld_sf, deforestation, fun=min), from = 1, to = 0)
deforestation <- merge(deforestation, wrld_rast)
deforestation.means <- terra::extract(deforestation, project(grid_filt_vect, crs(deforestation)), fun=mean, na.rm = TRUE)
colnames(deforestation.means) <- c("seqnum", "deforestation")
deforestation.means$seqnum <- grid_filt$seqnum

#WorldClim data
#download from https://www.worldclim.org/data/worldclim21.html
#https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/joc.5086
#30 seconds = ~1km resolution
#unzip the files to a folder named "WorldClim"
WorldClim_raster <- do.call(c, lapply(1:19, function(i) {
  rast(paste0("WorldClim/wc2.1_30s_bio_", i,".tif"))
}))
WorldClim_means <- lapply(1:19, function(i) {
  terra::extract(WorldClim_raster[[i]], project(grid_filt_vect, crs(WorldClim_raster)), fun=mean, na.rm = TRUE)
}) %>% reduce(full_join, by="ID")
colnames(WorldClim_means) <- c("seqnum", paste0("bio",1:19,"_mean"))

WorldClim_vars <- lapply(1:19, function(i) {
  terra::extract(WorldClim_raster[[i]], project(grid_filt_vect, crs(WorldClim_raster)), fun=var, na.rm = TRUE)
}) %>% reduce(full_join, by="ID")
colnames(WorldClim_vars) <- c("seqnum", paste0("bio",1:19,"_var"))
WorldClim_means$seqnum <- WorldClim_vars$seqnum <- grid_filt$seqnum

#elevation data
#https://pubs.er.usgs.gov/publication/ofr20111073
#https://www.usgs.gov/land-resources/eros/coastal-changes-and-impacts/gmted2010
#https://topotools.cr.usgs.gov/gmted_viewer/gmted2010_global_grids.php
#mean statistic, 30 second resolution = ~1km resolution
elev_raster <- rast("gmted2010/mn30_grd/w001001.adf")
wrld_vect <- vect(wrld_sf)
#set oceans to NA
elev_raster <- mask(elev_raster, wrld_vect)
elev_means <- terra::extract(elev_raster, project(grid_filt_vect, crs(elev_raster)), fun=mean, na.rm = TRUE)
colnames(elev_means) <- c("seqnum", "elev_mean")
elev_vars <- terra::extract(elev_raster, project(grid_filt_vect, crs(elev_raster)), fun=var, na.rm = TRUE)
colnames(elev_vars) <- c("seqnum", "elev_var")
elev_means$seqnum <- elev_vars$seqnum <- grid_filt$seqnum

#Get intersections####
# get range size for each species
# area in square meters
sf_use_s2(FALSE)
mammal_IUCN_clean <- st_make_valid(mammal_IUCN_filt) %>%
  filter(grepl("Extant", legend, fixed = TRUE) & category != "EX") %>%
  rename(binomial = "sci_name")
mammal_ranges <- mammal_IUCN_clean %>%
  group_by(binomial) %>%
  summarize(n_shapes = n(), geometry = st_union(geometry)) %>%
  mutate(range_size = st_area(geometry))
sf_use_s2(TRUE)
mammal_ranges <- left_join(mammal_ranges, mammal_traits_IUCN, by = c("binomial" = "IUCN_species"))

# get IUCN red list category
mammal_IUCN_cat <- mammal_IUCN_clean %>%
  st_drop_geometry() %>%
  dplyr::select(binomial, category) %>%
  unique() %>%
  mutate(category = factor(category, ordered = TRUE,
                           levels = c("DD","LC","NT","VU","EN","CR","EW","EX")))
mammal_ranges <- left_join(mammal_ranges, mammal_IUCN_cat, by = "binomial")

# find which species intersect with which grid cells
# this can require a bit of time and RAM
intrsct <- terra::intersect(grid_filt_vect, vect(mammal_ranges))

# Extract areas from polygon objects then attach as attribute
# area in square meters
mammal_grids <- as.data.frame(intrsct)
mammal_grids$area <- expanse(intrsct)
mammal_grids$range_size.log10 <- log10(as.numeric(mammal_grids$range_size))

#determine continent of each seqnum
# world comes from spData
conts <- world %>%
  st_wrap_dateline() %>%
  st_make_valid() %>%
  dplyr::select(continent) %>%
  drop_na() %>%
  group_by(continent) %>%
  summarize(geometry = st_union(geom)) %>%
  st_intersection(grid_filt, .)
land <- world %>%
  st_union()
conts$area <- st_area(conts)
mammal_grids <- merge(mammal_grids,
                      conts %>%
                        st_drop_geometry() %>%
                        group_by(seqnum) %>%
                        summarize(continent = continent[which.max(area)]) %>%
                        mutate(continent = factor(continent)),
                      by = "seqnum", all.x = TRUE)

#Grid cell data####
#calculate number of species overlapping with grid cell, size stats, diet stats, range stats
grid_data <- mammal_grids %>% group_by(seqnum, continent) %>%
  summarize(n = n(),
            mean_mass = mean(BodyMass.log10, na.rm = TRUE), median_mass = median(BodyMass.log10, na.rm = TRUE),
            min_mass = min(BodyMass.log10, na.rm = TRUE), max_mass = max(BodyMass.log10, na.rm = TRUE),
            n_mass = sum(!is.na(BodyMass.Value)), var_mass = var(BodyMass.log10, na.rm = TRUE),
            mass_disp = mean(dist(BodyMass.log10, method = "manhattan"), na.rm = TRUE), #body mass dispersion https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12667
            mass_disp_s = mean(dist(BodyMass.log10[size_cat == "small"], method = "manhattan"), na.rm = TRUE),
            mass_disp_m = mean(dist(BodyMass.log10[size_cat == "medium"], method = "manhattan"), na.rm = TRUE),
            mass_disp_l = mean(dist(BodyMass.log10[size_cat == "large"], method = "manhattan"), na.rm = TRUE),
            mass_disp_sm = mean(dist(BodyMass.log10[size_cat != "large"], method = "manhattan"), na.rm = TRUE),
            mass_disp_no_xl = mean(dist(BodyMass.log10[BodyMass.log10 < 5], method = "manhattan"), na.rm = TRUE),
            skew_mass = skewness(BodyMass.log10, na.rm = TRUE), kurt_mass = kurtosis(BodyMass.log10, na.rm = TRUE),
            n_small = sum(size_cat == "small", na.rm = TRUE), n_med = sum(size_cat == "medium", na.rm = TRUE),
            n_large = sum(size_cat == "large", na.rm = TRUE, na.rm = TRUE),
            n_xl = sum(size_cat2 == "xlarge", na.rm = TRUE), n_xxl = sum(size_cat2 == "xxlarge", na.rm = TRUE),
            n_xl_xxl = sum(BodyMass.log10 > 5, na.rm = TRUE),
            carn = sum(Trophic.Level == "Carnivore", na.rm = TRUE), omni = sum(Trophic.Level == "Omnivore", na.rm = TRUE),
            herb = sum(Trophic.Level == "Herbivore", na.rm = TRUE), n_diet = sum(!is.na(Plant.Per)),
            mean_plant = mean(Plant.Per, na.rm = TRUE), var_plant = var(Plant.Per, na.rm = TRUE),
            plant_disp = mean(dist(Plant.Per, method = "manhattan"), na.rm = TRUE),
            mean_range = mean(range_size.log10, na.rm = TRUE), var_range = var(range_size.log10, na.rm = TRUE),
            skew_range = skewness(range_size.log10, na.rm = TRUE), kurt_range = kurtosis(range_size.log10, na.rm = TRUE),
            #calculate lower quantile of log_range~log_mass relationship for each grid cell
            range_mass_slope = lm(range_size.log10~BodyMass.log10)$coefficients[[2]],
            range_mass_10 = tryCatch(rq(range_size.log10~BodyMass.log10, tau = .1)$coefficients[[2]], error = function(e) NA),
            range_mass_90 = tryCatch(rq(range_size.log10~BodyMass.log10, tau = .9)$coefficients[[2]], error = function(e) NA),
            range_disp = mean(dist(range_size.log10, method = "manhattan"), na.rm = TRUE),
            max_range = max(range_size.log10, na.rm = TRUE), min_range = min(range_size.log10, na.rm = TRUE),
            .groups = "drop")

#merge cell coordinates
grid_data <- grid_data %>%
  left_join(grid_coord, by = "seqnum")
grid_data <- cbind(grid_data, st_coordinates(grid_data$geometry))
# merge polygons
grid_data <- grid_data %>%
  left_join(grid_filt %>% rename(polygon = geometry), by = "seqnum")
#merge footprint data
grid_data <- grid_data %>%
  left_join(footprint.2009.means, by = "seqnum") %>%
  left_join(footprint.1993.means, by = "seqnum")
#merge deforestation data
grid_data <- grid_data %>%
  left_join(deforestation.means, by = "seqnum")
#merge climate data
grid_data <- grid_data %>%
  left_join(WorldClim_means, by = "seqnum") %>%
  left_join(WorldClim_vars, by = "seqnum")
#merge elevation data
grid_data <- grid_data %>%
  left_join(elev_means, by = "seqnum") %>%
  left_join(elev_vars, by = "seqnum")

#calculate proportions
grid_data[c("pcarn","pomni","pherb")] <- grid_data[, c("carn","omni","herb")]/grid_data$n_diet
grid_data[c("psmall","pmed","plarge")] <- grid_data[, c("n_small","n_med","n_large")]/grid_data$n_mass

# Filter out communities with fewer than 5 species
grid_data_sub <- grid_data %>% filter(n >= 5)

# convert geometry and coordinates to projected coordinates
# Interrupted Goode homolosine
xyz <- sf_project(grid_data_sub[, c("X", "Y")], to = "+proj=igh", from = "+proj=longlat +datum=WGS84")
grid_data_sub$X_igh <- xyz[, 1]
grid_data_sub$Y_igh <- xyz[, 2]
grid_data_sub$polygon_igh <- st_transform(grid_data_sub$polygon, crs = "+proj=igh")
# Robinson
xyz <- sf_project(grid_data_sub[, c("X", "Y")], to = "+proj=robin", from = "+proj=longlat +datum=WGS84")
grid_data_sub$X_robin <- xyz[, 1]
grid_data_sub$Y_robin <- xyz[, 2]
grid_data_sub$polygon_robin <- st_transform(grid_data_sub$polygon, crs = "+proj=robin")

# Setup for projected plots ####
## Goode homolosine ####
# projection outline in long-lat coordinates
crs_goode <- "+proj=igh"
lats <- c(
  90:-90, # right side down
  -90:0, 0:-90, # third cut bottom
  -90:0, 0:-90, # second cut bottom
  -90:0, 0:-90, # first cut bottom
  -90:90, # left side up
  90:0, 0:90, # cut top
  90 # close
)
longs <- c(
  rep(180, 181), # right side down
  rep(c(80.01, 79.99), each = 91), # third cut bottom
  rep(c(-19.99, -20.01), each = 91), # second cut bottom
  rep(c(-99.99, -100.01), each = 91), # first cut bottom
  rep(-180, 181), # left side up
  rep(c(-40.01, -39.99), each = 91), # cut top
  180 # close
)

goode_outline <-
  list(cbind(longs, lats)) %>%
  st_polygon() %>%
  st_sfc(
    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  )
# now we need to work in transformed coordinates, not in long-lat coordinates
goode_outline <- st_transform(goode_outline, crs = crs_goode)

# get the bounding box in transformed coordinates and expand by 10%
xlim <- st_bbox(goode_outline)[c("xmin", "xmax")] * 1.05
ylim <- st_bbox(goode_outline)[c("ymin", "ymax")] * 1.05

# turn into enclosing rectangle
goode_encl_rect <-
  list(
    cbind(
      c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]),
      c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1])
    )
  ) %>%
  st_polygon() %>%
  st_sfc(crs = crs_goode)

# calculate the area outside the earth outline as the difference
# between the enclosing rectangle and the earth outline
goode_without <- st_difference(goode_encl_rect, goode_outline)

## Robinson ####
lats <- c(
  90:-90, # right side down
  -90:90, # left side up
  90 # close
)
longs <- c(
  rep(180, 181), # right side down
  rep(-180, 181), # left side up
  180 # close
)
robin_outline <-
  list(cbind(longs, lats)) %>%
  st_polygon() %>%
  st_sfc(
    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  )
# now we need to work in transformed coordinates, not in long-lat coordinates
robin_outline <- st_transform(robin_outline, crs = "+proj=robin")

# get the bounding box in transformed coordinates and expand by 10%
xlim <- st_bbox(robin_outline)[c("xmin", "xmax")] * 1.05
ylim <- st_bbox(robin_outline)[c("ymin", "ymax")] * 1.05

# turn into enclosing rectangle
robin_encl_rect <-
  list(
    cbind(
      c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]),
      c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1])
    )
  ) %>%
  st_polygon() %>%
  st_sfc(crs = "+proj=robin")

# calculate the area outside the earth outline as the difference
# between the enclosing rectangle and the earth outline
robin_without <- st_difference(robin_encl_rect, robin_outline)

#Figure 2####
lmsumm1 <- summary(lm(mass_disp~n, data = grid_data_sub))
lmsumm2 <- summary(lm(kurt_mass~n, data = grid_data_sub))
lmsumm3 <- summary(lm(skew_mass~n, data = grid_data_sub))

g2_a <- ggplot(grid_data_sub, aes(x = n, y = mass_disp)) +
  geom_point(shape = 21, fill = "grey90", color = "grey30") +
  geom_smooth(method = lm, se = FALSE, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = expression("Mass Dispersion (log"[10]*"g)")) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines")) +
  annotate(geom = "text", x = 200, y = 1.95, hjust = 1, size = 10,
           label = deparse(bquote("p:"~.(format.pval(lmsumm1$coefficients[2,4], 1))*";"~R^2*":"~.(round(lmsumm1$r.squared, 3)))), parse = T)
g2_b <- ggplot(grid_data_sub, aes(x = n, y = kurt_mass)) +
  geom_point(shape = 21, fill = "grey90", color = "grey30") +
  geom_smooth(method = lm, se = FALSE, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = "Mass Kurtosis") +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines")) +
  annotate(geom = "text", x = 200, y = 4.8, hjust = 1, size = 10,
           label = deparse(bquote("p:"~.(format.pval(lmsumm2$coefficients[2,4], 1))*";"~R^2*":"~.(round(lmsumm2$r.squared, 3)))), parse = T)
g2_c <- ggplot(grid_data_sub, aes(x = n, y = skew_mass)) +
  geom_point(shape = 21, fill = "grey90", color = "grey30") +
  geom_smooth(method = lm, se = FALSE, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_x_continuous(name = "Community Species Richness") +
  scale_y_continuous(name = "Mass Skewness") +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines")) +
  annotate(geom = "text", x = 200, y = 1.4, hjust = 1, size = 10,
           label = deparse(bquote("p:"~.(format.pval(lmsumm3$coefficients[2,4], 1))*";"~R^2*":"~.(round(lmsumm3$r.squared, 3)))), parse = T)
gg <- ggarrange2(g2_a, g2_b, g2_c, ncol = 1, labels = c("(a)", "(b)", "(c)"),
                 label.args = list(gp = grid::gpar(font = 2, cex = 2)),  draw = FALSE)
ggsave("Mammal Ranges Stats by N.pdf", gg, width = 10, height = 15)

#Figure 1####
#mass moments and species richness
g1 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = st_graticule(crs = "+proj=robin", lon = seq(-180, 180, 30),
                              lat = seq(-90, 90, 30)), color = "gray70") +
  geom_sf(data = grid_data_sub, aes(geometry = polygon_robin, fill = mass_disp), color = NA) +
  geom_sf(data = robin_without, fill = "white", color = "NA") +
  geom_sf(data = robin_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = bquote(bold("(b)")~" Mass Dispersion")) +
  theme(panel.background = element_rect(fill = colsea, color = "white", linewidth = 1),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(.5,.5,0,.5), "cm"), panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = -.05)) +
  scale_fill_viridis(option = "plasma", guide = "none")
g1_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(mass_disp), fill = viridis(30, option = "plasma")) +
  labs(x = expression("Mass Dispersion (log"[10]*"g)"), y = "# of Communities") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,0,.5), "cm"))
g2 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = st_graticule(crs = "+proj=robin", lon = seq(-180, 180, 30),
                              lat = seq(-90, 90, 30)), color = "gray70") +
  geom_sf(data = grid_data_sub, aes(geometry = polygon_robin, fill = kurt_mass), color = NA) +
  geom_sf(data = robin_without, fill = "white", color = "NA") +
  geom_sf(data = robin_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = bquote(bold("(c)")~" Mass Kurtosis")) +
  theme(panel.background = element_rect(fill = colsea, color = "white", linewidth = 1),
        panel.grid.major = element_blank(), panel.border = element_blank(),
        plot.margin = unit(c(.5,.5,0,.5), "cm"),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = -.025)) +
  scale_fill_viridis(option = "plasma", guide = "none")
g2_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(kurt_mass), fill = viridis(30, option = "plasma")) +
  labs(x = "Mass Kurtosis", y = "# of Communities") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,0,.5), "cm"))
g3 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = st_graticule(crs = "+proj=robin", lon = seq(-180, 180, 30),
                              lat = seq(-90, 90, 30)), color = "gray70") +
  geom_sf(data = grid_data_sub, aes(geometry = polygon_robin, fill = skew_mass), color = NA) +
  geom_sf(data = robin_without, fill = "white", color = "NA") +
  geom_sf(data = robin_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = bquote(bold("(d)")~" Mass Skewness")) +
  theme(panel.background = element_rect(fill = colsea, color = "white", linewidth = 1),
        panel.grid.major = element_blank(), panel.border = element_blank(),
        plot.margin = unit(c(.5,.5,0,.5), "cm"),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = -.05)) +
  scale_fill_viridis(option = "plasma", guide = "none")
g3_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(skew_mass), fill = viridis(30, option = "plasma")) +
  labs(x = "Mass Skewness", y = "# of Communities") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,0,.5), "cm"))
g4 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = st_graticule(crs = "+proj=robin", lon = seq(-180, 180, 30),
                              lat = seq(-90, 90, 30)), color = "gray70") +
  geom_sf(data = grid_data_sub, aes(geometry = polygon_robin, fill = n), color = NA) +
  geom_sf(data = robin_without, fill = "white", color = "NA") +
  geom_sf(data = robin_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = bquote(bold("(a)")~" Species Richness")) +
  theme(panel.background = element_rect(fill = colsea, color = "white", linewidth = 1),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(.5,.5,0,.5), "cm"), panel.border = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = -.05)) +
  scale_fill_viridis(option = "plasma", guide = "none")
g4_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(n), fill = viridis(30, option = "plasma")) +
  labs(x = "Species Richness", y = "# of Communities") +
  coord_cartesian(expand = FALSE, xlim = c(0, NA)) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,0,.5), "cm"))
gg <- ggarrange2(g4,g4_hist,g1,g1_hist,g2,g2_hist,g3,g3_hist, ncol = 2, widths = c(2,1),  draw = FALSE)
ggsave("Mammal Ranges Mass Moments Maps.pdf", gg, width = 16, height = 25)

#Figure 1/2 combined####
gg <- ggarrange2(g1 + ggtitle(bquote(bold("(a)")~" Mass Dispersion")),
                 g1_hist,
                 g2_a + scale_x_continuous("Community Species Richness"),
                 g2 + ggtitle(bquote(bold("(b)")~" Mass Kurtosis")),
                 g2_hist,
                 g2_b + scale_x_continuous("Community Species Richness"),
                 g3 + ggtitle(bquote(bold("(c)")~" Mass Skewness")),
                 g3_hist,
                 g2_c,
                 ncol = 3, widths = c(2, 1, 1.5), draw = FALSE)
ggsave("Mammal Ranges Mass Moments Maps and Regressions.pdf", gg, width = 23, height = 18.75)
# then plot species richness separately
gg <- ggarrange2(g4 + ggtitle("Species Richness") + theme(plot.title = element_text(hjust = 0)),
                 g4_hist, ncol = 2, widths = c(2,1),  draw = FALSE)
ggsave("Mammal Ranges Species Richness Map.pdf", gg, width = 15.5, height = 6.25)

#diet and range
g1 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = mean_plant)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Mean of Plant Percentage of Diet") +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = "none")
g1_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(mean_plant), fill = viridis(30, option = "plasma")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
g2 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = plant_disp)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Dispersion of Plant Percentage of Diet") +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = "none")
g2_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(plant_disp), fill = viridis(30, option = "plasma")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
g3 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = mean_range)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = expression("Mean Range Size (log"[10]*"m"^2*")")) +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = "none")
g3_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(mean_range), fill = viridis(30, option = "plasma")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
g4 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = range_disp)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = expression("Range Size Dispersion (log"[10]*"m"^2*")")) +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = "none")
g4_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(range_disp), fill = viridis(30, option = "plasma")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
gg <- ggarrange2(g1,g1_hist,g2,g2_hist,g3,g3_hist,g4,g4_hist, ncol = 2, widths = c(2,1), draw = FALSE)
ggsave("Mammal Ranges Diet and Range Maps.pdf", gg, width = 17.5, height = 25)

#climate
g1 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = bio1_mean)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = expression("Mean Annual Temperature ("*degree*"C)")) +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = "none")
g1_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(bio1_mean), fill = viridis(30, option = "plasma")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
g2 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = bio4_mean)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = expression("Mean Temperature Seasonality ("*degree*"C * 100)")) +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = "none")
g2_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(bio4_mean), fill = viridis(30, option = "plasma")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
g3 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = bio12_mean)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Mean Annual Precipitation (mm)") +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = "none")
g3_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(bio12_mean), fill = viridis(30, option = "plasma")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
g4 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = bio15_mean)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Mean Precipitation Seasonality (%)") +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = "none")
g4_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(bio15_mean), fill = viridis(30, option = "plasma")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
gg <- ggarrange2(g1,g1_hist,g2,g2_hist,g3,g3_hist,g4,g4_hist, ncol = 2, widths = c(2,1), draw = FALSE)
ggsave("Mammal Ranges Climate Maps.pdf", gg, width = 17.5, height = 25)

#elevation and climate variability ("Microclimateness")
g1 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = elev_mean)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Mean Elevation (m)") +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = "none")
g1_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(elev_mean), fill = viridis(30, option = "plasma")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
g2 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = elev_var)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = expression("Elevation Variance (m"^2*")")) +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = "none")
g2_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(elev_var), fill = viridis(30, option = "plasma")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
g3 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = bio1_var)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = expression("Mean Annual Temperature Variance ("*degree*"C"^2*")")) +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = "none")
g3_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(bio1_var), fill = viridis(30, option = "plasma")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
g4 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = bio12_var)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = expression("Annual Precipitation Variance (mm"^2*")")) +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = "none")
g4_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(bio12_var), fill = viridis(30, option = "plasma")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
gg <- ggarrange2(g1,g1_hist,g2,g2_hist,g3,g3_hist,g4,g4_hist, ncol = 2, widths = c(2,1), draw = FALSE)
ggsave("Mammal Ranges Elevation and Microclimate Maps.pdf", gg, width = 17.5, height = 25)

#diet and body size pie charts on maps
diet_gather <- gather(grid_data_sub, "type", "value", c("carn","omni","herb"))
diet_gather$type <- factor(diet_gather$type, levels = c("herb","omni","carn"))
g1 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_arc_bar(data = diet_gather, aes(x0 = X, y0 = Y, r0 = 0, r = 1.1, amount = value,
                                       group = seqnum, fill = type), stat='pie') +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Dietary Composition of Communities") +
  theme(axis.text = element_text(colour = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_fill_viridis(name = "Diet", discrete = TRUE, guide = "none")
pdiet_gather <- gather(grid_data_sub, "type", "value", c("pcarn","pomni","pherb"))
pdiet_gather$type <- factor(pdiet_gather$type, levels = c("pherb","pomni","pcarn"))
g1_hist <- ggplot(pdiet_gather) +
  geom_histogram(aes(value * 100, fill = type), binwidth = 2.5, boundary = 25, show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_x_continuous(name = NULL, limits = c(0, 100)) +
  scale_y_continuous(name = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  facet_wrap(~type, ncol = 1, strip.position = "right", scales = "free_x",
             labeller = as_labeller(c("pcarn" = "% Carnivore","pherb" = "% Herbivore", "pomni" = "% Omnivore")))
gg <- ggarrange2(g1, g1_hist, nrow = 1, widths = c(7,1), draw = FALSE)
ggsave("Mammal Ranges Diet Pie Chart Map.pdf", gg, width = 30, height = 13.5)

#Figure 3####
pie_colors <- c("#21908CFF", "#ED7953FF", "#F0F921FF")
size_gather <- gather(grid_data_sub, "type", "value", c("n_small","n_med","n_large"))
size_gather$type <- factor(size_gather$type, levels = c("n_small","n_med","n_large"))
g2 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = st_graticule(crs = "+proj=igh", lon = seq(-180, 180, 30),
                              lat = seq(-90, 90, 30)), color = "gray70") +
  geom_sf(data = goode_without, fill = "white", color = "NA") +
  geom_sf(data = goode_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
  coord_sf(crs = crs_goode, expand = FALSE) +
  geom_arc_bar(data = size_gather, aes(x0 = X_igh, y0 = Y_igh, r0 = 0, r = 110000, amount = value,
                                       group = seqnum, fill = type), stat='pie', linewidth = .125) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Size Composition of Communities") +
  theme(panel.background = element_rect(fill = colsea, color = "white", linewidth = 1),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"),
        axis.ticks = element_blank(), axis.text = element_blank()) +
  scale_fill_manual(name = "Diet", values = pie_colors, guide = "none")
psize_gather <- gather(grid_data_sub, "type", "value", c("psmall","pmed","plarge"))
psize_gather$type <- factor(psize_gather$type, levels = c("psmall","pmed","plarge"))
g2_hist <- ggplot(psize_gather) +
  geom_histogram(aes(value * 100, fill = type), binwidth = 2.5, boundary = 25, show.legend = FALSE) +
  scale_fill_manual(values = pie_colors) +
  scale_x_continuous(name = NULL, limits = c(0, 100)) +
  scale_y_continuous(name = "# of Communities") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  facet_wrap(~type, ncol = 1, strip.position = "right", scales = "free_x",
             labeller = as_labeller(c("psmall" = "% <1kg","pmed" = "% 1kg - 10kg", "plarge" = "% >10kg")))
g2_hist_horiz <- ggplot(psize_gather) +
  geom_histogram(aes(value * 100, fill = type), binwidth = 2.5, boundary = 25, show.legend = FALSE) +
  scale_fill_manual(values = pie_colors) +
  scale_x_continuous(name = NULL, limits = c(0, 100)) +
  scale_y_continuous(name = "# of Communities") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"), panel.spacing.x = unit(.05, "npc")) +
  facet_wrap(~type, nrow = 1, strip.position = "top", scales = "free_x",
             labeller = as_labeller(c("psmall" = "% <1kg","pmed" = "% 1kg - 10kg", "plarge" = "% >10kg")))
gg <- ggarrange2(g2, g2_hist, nrow = 1, widths = c(7,1), draw = FALSE)
ggsave("Mammal Ranges Size Pie Chart Map_Goode.pdf", gg, width = 30, height = 12)
gg <- g2 / g2_hist_horiz + plot_layout(heights = c(3, 1))
ggsave("Mammal Ranges Size Pie Chart Map_Goode2.pdf", gg, width = 30, height = 21)

# Faceted map
g2_facet <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = st_graticule(crs = "+proj=igh", lon = seq(-180, 180, 30),
                              lat = seq(-90, 90, 30)) %>% dplyr::select(-type), color = "gray70") +
  geom_sf(data = psize_gather, aes(geometry = polygon_igh, fill = value * 100), linewidth = .125) +
  geom_sf(data = goode_without, fill = "white", color = "NA") +
  geom_sf(data = goode_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
  coord_sf(crs = crs_goode, expand = FALSE) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Size Composition of Communities") +
  theme(panel.background = element_rect(fill = colsea, color = "white", linewidth = 1),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(.5,0,.5,.5), "cm"),
        axis.ticks = element_blank(), axis.text = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = rel(2)),
        panel.spacing.y = unit(0.04, "npc")) +
  scale_fill_viridis(name = "% of community", guide = "none", option = "plasma") +
  facet_wrap(~type, ncol = 1, strip.position = "left",
             labeller = as_labeller(c("psmall" = "% <1kg","pmed" = "% 1kg - 10kg", "plarge" = "% >10kg")))
g2_hist_facet <- ggplot(psize_gather) +
  geom_histogram(aes(value * 100, fill = after_stat(x)), binwidth = 2.5, boundary = 25, show.legend = FALSE) +
  scale_fill_viridis(option = "plasma") +
  scale_x_continuous(name = "% of Community", limits = c(0, 100)) +
  scale_y_continuous(name = "# of Communities") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,0), "cm")) +
  facet_wrap(~type, ncol = 1, scales = "free_x", strip.position = "right",
             labeller = as_labeller(c("psmall" = "% <1kg","pmed" = "% 1kg - 10kg", "plarge" = "% >10kg")))
gg <- g2_facet + g2_hist_facet + plot_layout(widths = c(2.5,1))
ggsave("Mammal Ranges Size Pie Chart Map_Goode_Faceted.pdf", gg, width = 18, height = 16)

# With Robinson projection
g2 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = st_graticule(crs = "+proj=robin", lon = seq(-180, 180, 30),
                              lat = seq(-90, 90, 30)), color = "gray70") +
  geom_sf(data = robin_without, fill = "white", color = NA) +
  geom_sf(data = robin_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +
  geom_arc_bar(data = size_gather, aes(x0 = X_robin, y0 = Y_robin, r0 = 0, r = 110000, amount = value,
                                       group = seqnum, fill = type), stat='pie', linewidth = .125) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Size Composition of Communities") +
  theme(panel.background = element_rect(fill = colsea, color = "white", linewidth = 1),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"),
        axis.ticks = element_blank(), axis.text = element_blank()) +
  scale_fill_manual(name = "Diet", values = pie_colors, guide = "none")
gg <- ggarrange2(g2, g2_hist, nrow = 1, widths = c(7,1), draw = FALSE)
ggsave("Mammal Ranges Size Pie Chart Map_Robin.pdf", gg, width = 30, height = 13.5)
gg <- g2 / g2_hist_horiz + plot_layout(heights = c(3, 1))
ggsave("Mammal Ranges Size Pie Chart Map_Robin2.pdf", gg, width = 30, height = 21)

# Faceted plot
g2_facet <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = st_graticule(crs = "+proj=robin", lon = seq(-180, 180, 30),
                              lat = seq(-90, 90, 30)) %>% dplyr::select(-type), color = "gray70") +
  geom_sf(data = psize_gather, aes(geometry = polygon_robin, fill = value * 100), color = NA) +
  geom_sf(data = robin_without, fill = "white", color = "NA") +
  geom_sf(data = robin_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Size Composition of Communities") +
  theme(panel.background = element_rect(fill = colsea, color = "white", linewidth = 1),
        panel.grid.major = element_blank(), panel.border = element_blank(),
        plot.margin = unit(c(.5,0,.5,.5), "cm"),
        axis.ticks = element_blank(), axis.text = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = rel(2)),
        panel.spacing.y = unit(0.04, "npc")) +
  scale_fill_viridis(name = "% of community", guide = "none", option = "plasma") +
  facet_wrap(~factor(type, levels = c("plarge", "pmed", "psmall")), ncol = 1, strip.position = "left",
             labeller = as_labeller(c("psmall" = "% <1kg","pmed" = "% 1kg - 10kg", "plarge" = "% >10kg")))
g2_hist_facet <- ggplot(psize_gather) +
  geom_histogram(aes(value * 100, fill = after_stat(x)), binwidth = 2.5, boundary = 25, show.legend = FALSE) +
  scale_fill_viridis(option = "plasma") +
  scale_x_continuous(name = "% of Community", limits = c(0, 100)) +
  scale_y_continuous(name = "# of Communities") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,0), "cm")) +
  facet_wrap(~factor(type, levels = c("plarge", "pmed", "psmall")), ncol = 1, scales = "free_x", strip.position = "right",
             labeller = as_labeller(c("psmall" = "% <1kg","pmed" = "% 1kg - 10kg", "plarge" = "% >10kg")))
gg <- g2_facet + g2_hist_facet + plot_layout(widths = c(2.5,1))
ggsave("Mammal Ranges Size Pie Chart Map_Robin_Faceted.pdf", gg, width = 15.5, height = 16)

#Null model setup####
mammals_in_grids <- merge(unique(mammal_grids["binomial"]) %>% rename(IUCN_species = binomial), mammal_traits_IUCN, all.x = TRUE)

colors <- c(rgb(0,0,0), rgb(230/255, 159/255, 0/255), rgb(86/255, 180/255, 233/255), rgb(0/255, 158/255, 115/255), rgb(240/255,228/255,66/255), rgb(0/255, 114/255, 178/255), rgb(213/255,94/255,0/255), rgb(204/255,121/255,167/255))
# Ordered as in http://mkweb.bcgsc.ca/colorblind/img/colorblindness.palettes.trivial.png

#Null model with Diet, Ranges, and Continents Maintained####
mammals_in_grids_cont <- merge(unique(mammal_grids[c("binomial", "continent")]) %>% rename(IUCN_species = binomial), mammal_traits_IUCN, all.x = TRUE)
num_sim <- 100
sim_data_cont <- matrix(NA, nrow = nrow(grid_data_sub), ncol = num_sim)
pb <- txtProgressBar(min = 0, max = num_sim, style = 3)
for(i in 1:num_sim){
  setTxtProgressBar(pb, i)
  mammal_size_sim <- mammals_in_grids_cont %>% group_by(continent, Trophic.Level) %>% mutate(mass = sample(BodyMass.log10)) %>% rename(binomial = IUCN_species)
  mammal_grids_sim <- merge(mammal_grids %>% dplyr::select(-BodyMass.Value, -Trophic.Level), mammal_size_sim, by = c("binomial", "continent"))
  grid_data_sim <- mammal_grids_sim %>% group_by(seqnum) %>% summarise(n = n(), n_mass = sum(!is.na(mass)), mean = mean(mass, na.rm = TRUE), var = var(mass, na.rm = TRUE), mass_disp = mean(dist(mass, method = "manhattan"), na.rm = TRUE))
  grid_data_sim_sub <- grid_data_sim %>% filter(n >= 5)
  sim_data_cont[,i] <- grid_data_sim_sub$mass_disp
}
close(pb)

sim_results_cont <- cbind(data.frame(sim = rep(seq(1:num_sim), each = nrow(grid_data_sub)), sim_disp = as.vector(sim_data_cont), diff = as.vector(grid_data_sub$mass_disp - sim_data_cont)),
                          grid_data_sub[rep(seq_len(nrow(grid_data_sub)), num_sim), ])

grid_data_sub_cont <- subset(sim_results_cont, !is.na(mass_disp)) %>%
  group_by(across(c(-sim, -sim_disp, -diff))) %>%
  summarize(mean_sim_disp = mean(sim_disp, na.rm = TRUE), sd_sim_disp = sd(sim_disp, na.rm = TRUE),
            mean_diff = mean(diff, na.rm = TRUE),
            var_diff = var(diff, na.rm = TRUE),
            p_two = wilcox.test(sim_disp, mu = unique(mass_disp))$p.value,
            p_less = wilcox.test(sim_disp, mu = unique(mass_disp), alternative = "less")$p.value,
            p_greater = wilcox.test(sim_disp, mu = unique(mass_disp), alternative = "greater")$p.value,
            .groups = "drop")
grid_data_sub_cont$p_two_adjust <- p.adjust(grid_data_sub_cont$p_two, method = "fdr", n = nrow(grid_data_sub_cont))
grid_data_sub_cont$p_less_adjust <- p.adjust(grid_data_sub_cont$p_less, method = "fdr", n = nrow(grid_data_sub_cont))
grid_data_sub_cont$p_greater_adjust <- p.adjust(grid_data_sub_cont$p_greater, method = "fdr", n = nrow(grid_data_sub_cont))

table(cut(grid_data_sub_cont$p_less_adjust, c(0,0.05,1)), exclude = NULL)/nrow(grid_data_sub_cont)

#Null model 2 with Diet, Size Category, Ranges, and Continents Maintained####
num_sim <- 100
sim_data_cont2 <- matrix(NA, nrow = nrow(grid_data_sub), ncol = num_sim)
pb <- txtProgressBar(min = 0, max = num_sim, style = 3)
for(i in 1:num_sim){
  setTxtProgressBar(pb, i)
  mammal_size_sim <- mammals_in_grids_cont %>% group_by(continent, size_cat, Trophic.Level) %>% mutate(mass = sample(BodyMass.log10)) %>% rename(binomial = IUCN_species)
  mammal_grids_sim <- merge(mammal_grids %>% dplyr::select(-BodyMass.Value, -size_cat, -Trophic.Level), mammal_size_sim, by = c("binomial", "continent"))
  grid_data_sim <- mammal_grids_sim %>% group_by(seqnum) %>% summarise(n = n(), n_mass = sum(!is.na(mass)), mean = mean(mass, na.rm = TRUE), var = var(mass, na.rm = TRUE), mass_disp = mean(dist(mass, method = "manhattan"), na.rm = TRUE))
  grid_data_sim_sub <- grid_data_sim %>% filter(n >= 5)
  sim_data_cont2[,i] <- grid_data_sim_sub$mass_disp
}
close(pb)

sim_results_cont2 <- cbind(data.frame(sim = rep(seq(1:num_sim), each = nrow(grid_data_sub)), sim_disp = as.vector(sim_data_cont2), diff = as.vector(grid_data_sub$mass_disp - sim_data_cont2)),
                           grid_data_sub[rep(seq_len(nrow(grid_data_sub)), num_sim), ])

grid_data_sub_cont2 <- subset(sim_results_cont2, !is.na(mass_disp)) %>%
  group_by(across(c(-sim, -sim_disp, -diff))) %>%
  summarize(mean_sim_disp = mean(sim_disp, na.rm = TRUE), sd_sim_disp = sd(sim_disp, na.rm = TRUE),
            mean_diff = mean(diff, na.rm = TRUE),
            var_diff = var(diff, na.rm = TRUE),
            p_two = wilcox.test(sim_disp, mu = unique(mass_disp))$p.value,
            p_less = wilcox.test(sim_disp, mu = unique(mass_disp), alternative = "less")$p.value,
            p_greater = wilcox.test(sim_disp, mu = unique(mass_disp), alternative = "greater")$p.value,
            .groups = "drop")
grid_data_sub_cont2$p_two_adjust <- p.adjust(grid_data_sub_cont2$p_two, method = "fdr", n = nrow(grid_data_sub_cont2))
grid_data_sub_cont2$p_less_adjust <- p.adjust(grid_data_sub_cont2$p_less, method = "fdr", n = nrow(grid_data_sub_cont2))
grid_data_sub_cont2$p_greater_adjust <- p.adjust(grid_data_sub_cont2$p_greater, method = "fdr", n = nrow(grid_data_sub_cont2))

table(cut(grid_data_sub_cont2$p_less_adjust, c(0,0.05,1)), exclude = NULL)/nrow(grid_data_sub_cont2)

#Figure 5####
g1 <- ggplot(grid_data_sub_cont2) +
  geom_errorbar(aes(x = n, ymin = mean_sim_disp - 2 * sd_sim_disp,
                    ymax = mean_sim_disp + 2 * sd_sim_disp), width = 0, color = "grey70") +
  geom_point(aes(x = n, y = mean_sim_disp), shape = 21, fill = "grey100", show.legend = TRUE) +
  geom_point(aes(x = n, y = mass_disp, fill = cut(plarge, quantile(plarge, c(0, .5, 1)), include.lowest = TRUE)), shape = 21) +
  coord_cartesian(ylim = c(.5, 2.25), xlim = c(0, max(grid_data_sub_cont$n) + 2)) +
  scale_x_continuous(name = "Community Species Richness", expand = c(0.001, 0.001)) +
  scale_y_continuous(name = expression("Mass Dispersion (log"[10]*"g)")) +
  scale_fill_manual(name = NULL, labels = c("< 18.5% Large Species", "> 18.5% Large Species", "Null Model"),
                    values = viridis_pal()(5)[c(2,4)]) +
  theme_classic(base_size = 30) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position = c(.73, .92), legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0), strip.background = element_blank(),
        legend.background = element_rect(fill = NA)) +
  guides(fill = guide_legend(override.aes = list("size" = 3)))
g2 <- ggplot(grid_data_sub_cont2) +
  geom_abline(slope = 1, intercept = 0) +
  geom_errorbar(aes(x = mass_disp, ymin = mean_sim_disp - 2 * sd_sim_disp,
                    ymax = mean_sim_disp + 2 * sd_sim_disp), width = 0, color = "grey70") +
  geom_point(aes(x = mass_disp, y = mean_sim_disp, fill = cut(plarge, quantile(plarge, c(0, .5, 1)), include.lowest = TRUE)), shape = 21) +
  scale_x_continuous(name = expression("Observed Mass Dispersion (log"[10]*"g)")) +
  scale_y_continuous(name = expression("Null Mass Dispersion (log"[10]*"g)")) +
  scale_fill_manual(name = NULL, labels = c("< 18.5% Large Species", "> 18.5% Large Species", "Null Model"),
                    values = viridis_pal()(5)[c(2,4)]) +
  coord_cartesian(xlim = c(.5, 2.7), ylim = c(.5, 2.7)) +
  theme_classic(base_size = 30) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position = c(.27, .92), legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0), strip.background = element_blank(),
        legend.background = element_rect(fill = NA)) +
  guides(fill = guide_legend(override.aes = list("size" = 3)))
gg <- ggarrange2(g1, g2, ncol = 2, widths = c(1,1), labels = c("(a)", "(b)"),
                 label.args = list(gp = grid::gpar(font = 2, cex = 2.5)),  draw = FALSE)
ggsave("Mammal Ranges Simulation Variance - Ranges and Diet Maintained.pdf", gg, width = 20, height = 10)

#Figure 4####
#plot simulation value and simulation difference on map
g1 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = st_graticule(crs = "+proj=robin", lon = seq(-180, 180, 30),
                              lat = seq(-90, 90, 30)), color = "gray70") +
  geom_sf(data = robin_without, fill = "white", color = NA) +
  geom_sf(data = robin_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
  geom_sf(data = grid_data_sub_cont2, aes(geometry = polygon_robin, fill = mean_sim_disp), color = NA) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = bquote(bold("(a)")~" Null Model Mean Mass Dispersion")) +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea, color = "white", linewidth = 1),
        panel.grid.major = element_blank(), plot.margin = unit(c(.5,.5,.5,.5), "cm"),
        panel.border = element_blank()) +
  scale_fill_viridis(option = "plasma", guide = "none")
g1_hist <- ggplot(data = grid_data_sub_cont2) +
  geom_histogram(aes(mean_sim_disp), fill = viridis(30, option = "plasma")) +
  labs(x = expression("Null Mass Dispersion (log"[10]*"g)"), y = "# of Communities") +
  coord_cartesian(expand =FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))

g2 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = st_graticule(crs = "+proj=robin", lon = seq(-180, 180, 30),
                              lat = seq(-90, 90, 30)), color = "gray70") +
  geom_sf(data = robin_without, fill = "white", color = NA) +
  geom_sf(data = robin_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
  geom_sf(data = grid_data_sub_cont2, aes(geometry = polygon_robin, fill = mean_diff), color = NA) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = bquote(bold("(b)")~" Observed Mass Dispersion - Null Model Mean Mass Dispersion")) +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea, color = "white", linewidth = 1),
        panel.grid.major = element_blank(), plot.margin = unit(c(.5,.5,.5,.5), "cm"),
        panel.border = element_blank(), ) +
  scale_fill_viridis(option = "plasma", guide = "none")
g2_hist <- ggplot(data = grid_data_sub_cont2) +
  geom_histogram(aes(mean_diff), fill = viridis(30, option = "plasma")) +
  labs(x = expression("Deviation from the Null (log"[10]*"g)"), y = "# of Communities") +
  coord_cartesian(expand =FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
gg <- ggarrange2(g1, g1_hist, g2, g2_hist, ncol = 2, widths = c(2,1),  draw = FALSE)
ggsave("Mammal Ranges Simulation Values and Deviations - Ranges and Diet Maintained.pdf", gg, width = 19, height = 15)


#Null model 3 without big things####
num_sim <- 100
grid_data_sub_sm <- grid_data_sub %>% filter((n - n_large) >= 5)
sim_data_cont3 <- matrix(NA, nrow = nrow(grid_data_sub_sm), ncol = num_sim)
pb <- txtProgressBar(min = 0, max = num_sim, style = 3)
for(i in 1:num_sim){
  setTxtProgressBar(pb, i)
  mammal_size_sim <- mammals_in_grids_cont %>% filter(size_cat != "large") %>% group_by(continent, size_cat, Trophic.Level) %>% mutate(mass = sample(BodyMass.log10)) %>% rename(binomial = IUCN_species)
  mammal_grids_sim <- merge(mammal_grids %>% dplyr::select(-BodyMass.Value, -size_cat, -Trophic.Level), mammal_size_sim, by = c("binomial", "continent"))
  grid_data_sim <- mammal_grids_sim %>% group_by(seqnum) %>% summarise(n = n(), n_mass = sum(!is.na(mass)), mean = mean(mass, na.rm = TRUE), var = var(mass, na.rm = TRUE), mass_disp = mean(dist(mass, method = "manhattan"), na.rm = TRUE))
  grid_data_sim_sub <- grid_data_sim %>% filter(n >= 5)
  sim_data_cont3[,i] <- grid_data_sim_sub$mass_disp
}
close(pb)

sim_results_cont3 <- cbind(data.frame(sim = rep(seq(1:num_sim), each = nrow(grid_data_sub_sm)), sim_disp = as.vector(sim_data_cont3), diff = as.vector(grid_data_sub_sm$mass_disp_sm - sim_data_cont3)),
                           grid_data_sub_sm[rep(seq_len(nrow(grid_data_sub_sm)), num_sim), ])

grid_data_sub_cont3 <- subset(sim_results_cont3, !is.na(mass_disp_sm)) %>%
  group_by(across(c(-sim, -sim_disp, -diff))) %>%
  summarize(mean_sim_disp = mean(sim_disp, na.rm = TRUE), sd_sim_disp = sd(sim_disp, na.rm = TRUE),
            mean_diff = mean(diff, na.rm = TRUE),
            var_diff = var(diff, na.rm = TRUE),
            p_two = wilcox.test(sim_disp, mu = unique(mass_disp_sm))$p.value,
            p_less = wilcox.test(sim_disp, mu = unique(mass_disp_sm), alternative = "less")$p.value,
            p_greater = wilcox.test(sim_disp, mu = unique(mass_disp_sm), alternative = "greater")$p.value,
            .groups = "drop")
grid_data_sub_cont3$p_two_adjust <- p.adjust(grid_data_sub_cont3$p_two, method = "fdr", n = nrow(grid_data_sub_cont3))
grid_data_sub_cont3$p_less_adjust <- p.adjust(grid_data_sub_cont3$p_less, method = "fdr", n = nrow(grid_data_sub_cont3))
grid_data_sub_cont3$p_greater_adjust <- p.adjust(grid_data_sub_cont3$p_greater, method = "fdr", n = nrow(grid_data_sub_cont3))

table(cut(grid_data_sub_cont3$p_less_adjust, c(0,0.05,1)), exclude = NULL)/nrow(grid_data_sub_cont3)

sm_summ_xl_xxl <- grid_data_sub_cont3 %>%
  group_by(n_xl_xxl) %>%
  summarize(higher = sum(mean_diff > 0 & p_less_adjust < 0.05)/n(),
            lower = sum(mean_diff < 0 & p_greater_adjust < 0.05)/n(), n = n())

sm_summ_xxl <- grid_data_sub_cont3 %>%
  group_by(n_xxl) %>%
  summarize(higher = sum(mean_diff > 0 & p_less_adjust < 0.05)/n(),
            lower = sum(mean_diff < 0 & p_greater_adjust < 0.05)/n(), n = n())

ggplot(grid_data_sub_cont3) +
  geom_violin(aes(x = as.factor(n_xl_xxl), y = mean_diff)) +
  geom_text(data = sm_summ_xl_xxl, aes(x = as.factor(n_xl_xxl), label = round(higher, 2)), y = 0.3) +
  geom_text(data = sm_summ_xl_xxl, aes(x = as.factor(n_xl_xxl), label = round(lower, 2)), y = -0.3) +
  scale_y_continuous("Mean Mass Dispersion Difference (Observed - Null)") +
  scale_x_discrete("# of Mammal Species > 100 kg") +
  coord_cartesian(ylim = c(-.3, .3))

ggplot(grid_data_sub_cont2) +
  geom_violin(aes(x = as.factor(n_xxl), y = mean_diff)) +
  geom_text(data = sm_summ_xxl, aes(x = as.factor(n_xxl), label = round(higher, 2)), y = 0.3) +
  geom_text(data = sm_summ_xxl, aes(x = as.factor(n_xxl), label = round(lower, 2)), y = -0.3) +
  scale_y_continuous("Mean Mass Dispersion Difference (Observed - Null)") +
  scale_x_discrete("# of Mammal Species > 1000 kg") +
  coord_cartesian(ylim = c(-.3, .3))

# plots
g1 <- ggplot(grid_data_sub_cont3) +
  geom_errorbar(aes(x = n, ymin = mean_sim_disp - 2 * sd_sim_disp,
                    ymax = mean_sim_disp + 2 * sd_sim_disp), width = 0, color = "grey70") +
  geom_point(aes(x = n, y = mean_sim_disp), shape = 21, fill = "grey100", show.legend = TRUE) +
  geom_point(aes(x = n, y = mass_disp_sm, fill = cut(plarge, quantile(plarge, c(0, .5, 1)), include.lowest = TRUE)), shape = 21) +
  coord_cartesian(ylim = c(.4, 2), xlim = c(0, max(grid_data_sub_cont$n) + 2)) +
  scale_x_continuous(name = "Species Richness", expand = c(0.001, 0.001)) +
  scale_y_continuous(name = expression("Mass Dispersion (log"[10]*"g)")) +
  scale_fill_manual(name = NULL, labels = c("< 18.5% Large Species", "> 18.5% Large Species", "Null Model"),
                    values = viridis_pal()(5)[c(2,4)]) +
  theme_classic(base_size = 30) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position = c(.73, .95), legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0), strip.background = element_blank(),
        legend.background = element_rect(fill = NA)) +
  guides(fill = guide_legend(override.aes = list("size" = 3)))

g2 <- ggplot(grid_data_sub_cont3) +
  geom_abline(slope = 1, intercept = 0) +
  geom_errorbar(aes(x = mass_disp_sm, ymin = mean_sim_disp - 2 * sd_sim_disp,
                    ymax = mean_sim_disp + 2 * sd_sim_disp), width = 0, color = "grey70") +
  geom_point(aes(x = mass_disp_sm, y = mean_sim_disp, fill = cut(plarge, quantile(plarge, c(0, .5, 1)), include.lowest = TRUE)), shape = 21, show.legend = TRUE) +
  #coord_cartesian(ylim = c(.5, 2.25), xlim = c(0, max(grid_data_sub_cont$n) + 2)) +
  scale_x_continuous(name = expression("Observed Mass Dispersion (log"[10]*"g)")) +
  scale_y_continuous(name = expression("Null Mass Dispersion (log"[10]*"g)")) +
  scale_fill_manual(name = NULL, labels = c("< 17.5% Large Species", "> 17.5% Large Species", "Null Model"),
                    values = viridis_pal()(5)[c(2,4)]) +
  theme_classic(base_size = 30) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position = c(.73, .95), legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0), strip.background = element_blank(),
        legend.background = element_rect(fill = NA)) +
  guides(fill = guide_legend(override.aes = list("size" = 3)))
gg <- ggarrange2(g1, g2, ncol = 2, widths = c(1,1),  draw = FALSE)
ggsave("Mammal Ranges Null 2 no large.pdf", gg, width = 20, height = 10)

#Null model 4 only small things####
num_sim <- 100
grid_data_sub_s <- grid_data_sub %>% filter(n_small >= 5)
sim_data_cont4 <- matrix(NA, nrow = nrow(grid_data_sub_s), ncol = num_sim)
pb <- txtProgressBar(min = 0, max = num_sim, style = 3)
for(i in 1:num_sim){
  setTxtProgressBar(pb, i)
  mammal_size_sim <- mammals_in_grids_cont %>% filter(size_cat == "small") %>% group_by(continent, size_cat, Trophic.Level) %>% mutate(mass = sample(BodyMass.log10)) %>% rename(binomial = IUCN_species)
  mammal_grids_sim <- merge(mammal_grids %>% dplyr::select(-BodyMass.Value, -size_cat, -Trophic.Level), mammal_size_sim, by = c("binomial", "continent"))
  grid_data_sim <- mammal_grids_sim %>% group_by(seqnum) %>% summarise(n = n(), n_mass = sum(!is.na(mass)), mean = mean(mass, na.rm = TRUE), var = var(mass, na.rm = TRUE), mass_disp = mean(dist(mass, method = "manhattan"), na.rm = TRUE))
  grid_data_sim_sub <- grid_data_sim %>% filter(n >= 5)
  sim_data_cont4[,i] <- grid_data_sim_sub$mass_disp
}
close(pb)

sim_results_cont4 <- cbind(data.frame(sim = rep(seq(1:num_sim), each = nrow(grid_data_sub_s)), sim_disp = as.vector(sim_data_cont4), diff = as.vector(grid_data_sub_s$mass_disp_s - sim_data_cont4)),
                           grid_data_sub_s[rep(seq_len(nrow(grid_data_sub_s)), num_sim), ])

grid_data_sub_cont4 <- subset(sim_results_cont4, !is.na(mass_disp_s)) %>%
  group_by(across(c(-sim, -sim_disp, -diff))) %>%
  summarize(mean_sim_disp = mean(sim_disp, na.rm = TRUE), sd_sim_disp = sd(sim_disp, na.rm = TRUE),
            mean_diff = mean(diff, na.rm = TRUE),
            var_diff = var(diff, na.rm = TRUE),
            p_two = wilcox.test(sim_disp, mu = unique(mass_disp_s))$p.value,
            p_less = wilcox.test(sim_disp, mu = unique(mass_disp_s), alternative = "less")$p.value,
            p_greater = wilcox.test(sim_disp, mu = unique(mass_disp_s), alternative = "greater")$p.value,
            .groups = "drop")
grid_data_sub_cont4$p_two_adjust <- p.adjust(grid_data_sub_cont4$p_two, method = "fdr", n = nrow(grid_data_sub_cont4))
grid_data_sub_cont4$p_less_adjust <- p.adjust(grid_data_sub_cont4$p_less, method = "fdr", n = nrow(grid_data_sub_cont4))
grid_data_sub_cont4$p_greater_adjust <- p.adjust(grid_data_sub_cont4$p_greater, method = "fdr", n = nrow(grid_data_sub_cont4))

table(cut(grid_data_sub_cont4$p_less_adjust, c(0,0.05,1)), exclude = NULL)/nrow(grid_data_sub_cont3)

small_summ_xl_xxl <- grid_data_sub_cont4 %>%
  group_by(n_xl_xxl) %>%
  summarize(higher = sum(mean_diff > 0 & p_less_adjust < 0.05)/n(),
            lower = sum(mean_diff < 0 & p_greater_adjust < 0.05)/n(), n = n())

small_summ_xxl <- grid_data_sub_cont4 %>%
  group_by(n_xxl) %>%
  summarize(higher = sum(mean_diff > 0 & p_less_adjust < 0.05)/n(),
            lower = sum(mean_diff < 0 & p_greater_adjust < 0.05)/n(), n = n())

ggplot(grid_data_sub_cont4) +
  geom_violin(aes(x = as.factor(n_xl_xxl), y = mean_diff)) +
  geom_text(data = small_summ_xl_xxl, aes(x = as.factor(n_xl_xxl), label = round(higher, 2)), y = 0.3) +
  geom_text(data = small_summ_xl_xxl, aes(x = as.factor(n_xl_xxl), label = round(lower, 2)), y = -0.3) +
  scale_y_continuous("Mean Mass Dispersion Difference (Observed - Null)") +
  scale_x_discrete("# of Mammal Species > 100 kg") +
  coord_cartesian(ylim = c(-.3, .3))

ggplot(grid_data_sub_cont4) +
  geom_violin(aes(x = as.factor(n_xxl), y = mean_diff)) +
  geom_text(data = small_summ_xxl, aes(x = as.factor(n_xxl), label = round(higher, 2)), y = 0.3) +
  geom_text(data = small_summ_xxl, aes(x = as.factor(n_xxl), label = round(lower, 2)), y = -0.3) +
  scale_y_continuous("Mean Mass Dispersion Difference (Observed - Null)") +
  scale_x_discrete("# of Mammal Species > 1000 kg") +
  coord_cartesian(ylim = c(-.3, .3))

# plots
g1 <- ggplot(grid_data_sub_cont4) +
  geom_errorbar(aes(x = n, ymin = mean_sim_disp - 2 * sd_sim_disp,
                    ymax = mean_sim_disp + 2 * sd_sim_disp), width = 0, color = "grey70") +
  geom_point(aes(x = n, y = mean_sim_disp), shape = 21, fill = "grey100", show.legend = TRUE) +
  geom_point(aes(x = n, y = mass_disp_s, fill = cut(plarge, quantile(plarge, c(0, .5, 1)), include.lowest = TRUE)), shape = 21) +
  coord_cartesian(ylim = c(0, 1.3), xlim = c(0, max(grid_data_sub_cont$n) + 2)) +
  scale_x_continuous(name = "Species Richness", expand = c(0.001, 0.001)) +
  scale_y_continuous(name = expression("Mass Dispersion (log"[10]*"g)")) +
  scale_fill_manual(name = NULL, labels = c("< 18.5% Large Species", "> 18.5% Large Species", "Null Model"),
                    values = viridis_pal()(5)[c(2,4)]) +
  theme_classic(base_size = 30) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position = c(.73, .95), legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0), strip.background = element_blank(),
        legend.background = element_rect(fill = NA)) +
  guides(fill = guide_legend(override.aes = list("size" = 3)))

g2 <- ggplot(grid_data_sub_cont4) +
  geom_abline(slope = 1, intercept = 0) +
  geom_errorbar(aes(x = mass_disp_s, ymin = mean_sim_disp - 2 * sd_sim_disp,
                    ymax = mean_sim_disp + 2 * sd_sim_disp), width = 0, color = "grey70") +
  geom_point(aes(x = mass_disp_s, y = mean_sim_disp, fill = cut(plarge, quantile(plarge, c(0, .5, 1)), include.lowest = TRUE)), shape = 21, show.legend = TRUE) +
  #coord_cartesian(ylim = c(.5, 2.25), xlim = c(0, max(grid_data_sub_cont$n) + 2)) +
  scale_x_continuous(name = expression("Observed Mass Dispersion (log"[10]*"g)")) +
  scale_y_continuous(name = expression("Null Mass Dispersion (log"[10]*"g)")) +
  scale_fill_manual(name = NULL, labels = c("< 18.5% Large Species", "> 18.5% Large Species", "Null Model"),
                    values = viridis_pal()(5)[c(2,4)]) +
  theme_classic(base_size = 30) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position = c(.73, .95), legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0), strip.background = element_blank(),
        legend.background = element_rect(fill = NA)) +
  guides(fill = guide_legend(override.aes = list("size" = 3)))
gg <- ggarrange2(g1, g2, ncol = 2, widths = c(1,1),  draw = FALSE)
ggsave("Mammal Ranges Null 2 small only.pdf", gg, width = 20, height = 10)

#Null model 5 only medium things####
num_sim <- 100
grid_data_sub_m <- grid_data_sub %>% filter(n_med >= 5)
sim_data_cont5 <- matrix(NA, nrow = nrow(grid_data_sub_m), ncol = num_sim)
pb <- txtProgressBar(min = 0, max = num_sim, style = 3)
for(i in 1:num_sim){
  setTxtProgressBar(pb, i)
  mammal_size_sim <- mammals_in_grids_cont %>% filter(size_cat == "medium") %>% group_by(continent, size_cat, Trophic.Level) %>% mutate(mass = sample(BodyMass.log10)) %>% rename(binomial = IUCN_species)
  mammal_grids_sim <- merge(mammal_grids %>% dplyr::select(-BodyMass.Value, -size_cat, -Trophic.Level), mammal_size_sim, by = c("binomial", "continent"))
  grid_data_sim <- mammal_grids_sim %>% group_by(seqnum) %>% summarise(n = n(), n_mass = sum(!is.na(mass)), mean = mean(mass, na.rm = TRUE), var = var(mass, na.rm = TRUE), mass_disp = mean(dist(mass, method = "manhattan"), na.rm = TRUE))
  grid_data_sim_sub <- grid_data_sim %>% filter(n >= 5)
  sim_data_cont5[,i] <- grid_data_sim_sub$mass_disp
}
close(pb)

sim_results_cont5 <- cbind(data.frame(sim = rep(seq(1:num_sim), each = nrow(grid_data_sub_m)), sim_disp = as.vector(sim_data_cont5), diff = as.vector(grid_data_sub_m$mass_disp_m - sim_data_cont5)),
                           grid_data_sub_m[rep(seq_len(nrow(grid_data_sub_m)), num_sim), ])

grid_data_sub_cont5 <- subset(sim_results_cont5, !is.na(mass_disp_m)) %>%
  group_by(across(c(-sim, -sim_disp, -diff))) %>%
  summarize(mean_sim_disp = mean(sim_disp, na.rm = TRUE), sd_sim_disp = sd(sim_disp, na.rm = TRUE),
            mean_diff = mean(diff, na.rm = TRUE),
            var_diff = var(diff, na.rm = TRUE),
            p_two = wilcox.test(sim_disp, mu = unique(mass_disp_m))$p.value,
            p_less = wilcox.test(sim_disp, mu = unique(mass_disp_m), alternative = "less")$p.value,
            p_greater = wilcox.test(sim_disp, mu = unique(mass_disp_m), alternative = "greater")$p.value,
            .groups = "drop")
grid_data_sub_cont5$p_two_adjust <- p.adjust(grid_data_sub_cont5$p_two, method = "fdr", n = nrow(grid_data_sub_cont5))
grid_data_sub_cont5$p_less_adjust <- p.adjust(grid_data_sub_cont5$p_less, method = "fdr", n = nrow(grid_data_sub_cont5))
grid_data_sub_cont5$p_greater_adjust <- p.adjust(grid_data_sub_cont5$p_greater, method = "fdr", n = nrow(grid_data_sub_cont5))

table(cut(grid_data_sub_cont5$p_less_adjust, c(0,0.05,1)), exclude = NULL)/nrow(grid_data_sub_cont3)

med_summ_xl_xxl <- grid_data_sub_cont5 %>%
  group_by(n_xl_xxl) %>%
  summarize(higher = sum(mean_diff > 0 & p_less_adjust < 0.05)/n(),
            lower = sum(mean_diff < 0 & p_greater_adjust < 0.05)/n(), n = n())

med_summ_xxl <- grid_data_sub_cont5 %>%
  group_by(n_xxl) %>%
  summarize(higher = sum(mean_diff > 0 & p_less_adjust < 0.05)/n(),
            lower = sum(mean_diff < 0 & p_greater_adjust < 0.05)/n(), n = n())

ggplot(grid_data_sub_cont5) +
  geom_violin(aes(x = as.factor(n_xl_xxl), y = mean_diff)) +
  geom_text(data = med_summ_xl_xxl, aes(x = as.factor(n_xl_xxl), label = round(higher, 2)), y = 0.3) +
  geom_text(data = med_summ_xl_xxl, aes(x = as.factor(n_xl_xxl), label = round(lower, 2)), y = -0.3) +
  scale_y_continuous("Mean Mass Dispersion Difference (Observed - Null)") +
  scale_x_discrete("# of Mammal Species > 100 kg") +
  coord_cartesian(ylim = c(-.3, .3))

ggplot(grid_data_sub_cont5) +
  geom_violin(aes(x = as.factor(n_xxl), y = mean_diff)) +
  geom_text(data = med_summ_xxl, aes(x = as.factor(n_xxl), label = round(higher, 2)), y = 0.3) +
  geom_text(data = med_summ_xxl, aes(x = as.factor(n_xxl), label = round(lower, 2)), y = -0.3) +
  scale_y_continuous("Mean Mass Dispersion Difference (Observed - Null)") +
  scale_x_discrete("# of Mammal Species > 1000 kg") +
  coord_cartesian(ylim = c(-.3, .3))

# Correlations####
# check correlations between variables
var_cors <- cor(grid_data_sub_clean[, c("n", "n_xl_xxl", "psmall", "plarge", "pcarn", "pherb",
                                        "plant_disp", "mean_plant", "mean_range", "range_disp",
                                        "footprint.2009", "deforestation",
                                        "bio1_mean", "bio4_mean", "bio5_mean", "bio6_mean",
                                        "bio12_mean", "bio15_mean", "elev_mean")])
colnames(var_cors) <- gsub("\n", " ", labels[colnames(var_cors)])
rownames(var_cors) <- gsub("\n", " ", labels[rownames(var_cors)])
var_cors_filt <- var_cors[upper.tri(var_cors)]
table(var_cors_filt > 0)
summary(var_cors_filt[var_cors_filt > 0])
summary(var_cors_filt[var_cors_filt < 0])
table(abs(var_cors_filt) > 0.5)
write.csv(var_cors, "correlations.csv")

# Relative Importance ####
# What are the best predictors of variance in general?
labels <- c("continent" = "Continent",
            "continentAsia" = "Asia", "continentOceania" = "Oceania", "continentEurope" = "Europe",
            "continentNorth America" = "North America", "continentSouth America" = "South America",
            "n" = "Species Richness", "n_xl_xxl" = "# Species > 100kg",
            "max_mass" = "Maximum Mass*", "min_mass" = "Minimum Mass*",
            "psmall" = "% Small-Sized*", "pmed" = "% Medium-Sized*",
            "plarge" = "% Large-Sized*", "mean_plant" = "Mean % Plant\nin Diet",
            "plant_disp" = "Dispersion of %\nPlant in Diet", "pherb" = "% Herbivore",
            "pomni" = "% Omnivore", "pcarn" = "% Carnivore",
            "mean_range" = "Mean Range Size", "range_disp" = "Dispersion of\nRange Size",
            "footprint.2009" = "Human Footprint", "deforestation" = "Deforestation",
            "elev_mean" = "Mean Elevation", "elev_var" = "Variance of\nElevation",
            "bio1_mean" = "Mean Annual Temp", "bio1_var" = "Variance of Mean\nAnnual Temp",
            "bio4_mean" = "Mean Temp\nSeasonality", "bio4_var" = "Variance of\nTemp Seasonality",
            "bio5_mean" = "Mean Max Temp", "bio5_var" = "Variance of\nMax Temp",
            "bio6_mean" = "Mean Min Temp", "bio6_var" = "Variance of\nMin Temp",
            "bio12_mean" = "Mean Annual\nPrecip", "bio12_var" = "Variance of\nAnnual Precip",
            "bio15_mean" = "Mean Precip\nSeasonality", "bio15_var" = "Variance of\nPrecip Seasonality")

## Observed dispersion ####
#(should avoid things that are relevant to calculating variance)
grid_data_sub_clean <- subset(grid_data_sub)
reg1 <- lm(mass_disp ~ n + n_xl_xxl + pcarn + pherb + plant_disp + mean_plant + mean_range + range_disp +
             continent + footprint.2009 + deforestation +
             bio1_mean + bio4_mean + bio5_mean + bio6_mean + bio12_mean + bio15_mean + elev_mean,
          #    s(Y, X, bs="sos", k=60),
           #   bio1_var + bio4_var + bio5_var + bio6_var +
           #   bio12_var + bio15_var + elev_var,
           data = grid_data_sub_clean, na.action = "na.omit")
relimp1 <- calc.relimp(reg1, type = "lmg")
sort(relimp1@lmg, decreasing = TRUE)

## Deviation from null ####
# remove pmed and pomni to prevent singularity issues
grid_data_sub_cont2_clean <- subset(grid_data_sub_cont2)
reg2 <- lm(mean_diff ~ n + n_xl_xxl + psmall + plarge + pcarn + pherb +
             plant_disp + mean_plant + mean_range + range_disp + continent + footprint.2009 + deforestation +
             bio1_mean + bio4_mean + bio5_mean + bio6_mean + bio12_mean + bio15_mean + elev_mean,
             #bio1_var + bio4_var + bio5_var + bio6_var + bio12_var + bio15_var + elev_var,
           data = grid_data_sub_cont2_clean, weights = 1/var_diff, na.action = "na.omit")
relimp2 <- calc.relimp(reg2, type = "lmg")
sort(relimp2@lmg, decreasing = TRUE)

## Figure 6 ####
relimp_df <- rbind(data.frame(var = relimp1@namen[2:length(relimp1@namen)],
                              data = "Raw Dispersion", lmg = relimp1@lmg),
                   data.frame(var = relimp2@namen[2:length(relimp2@namen)],
                              data = "Deviation From Null", lmg = relimp2@lmg)) %>%
  mutate(group = ifelse(grepl("bio.+_mean|elev_mean", var), "habitat",
                        ifelse(grepl("bio.+_var|elev_var", var), "variance", "community"))) %>%
  arrange(group, -lmg)

rects <- data.frame(xmin = c(0.5, 7.5), xmax = c(7.5, 20.5),
                    ymin = 0.2, ymax = 0.25)

texts <- data.frame(x = c(4, 14), y = 0.225,
                    label = c("Abiotic\nEnvironment", "Biotic Community"))

# without bootstrap limits
ggplot(relimp_df, aes(x = var, y = lmg, group = data, fill = data)) +
  geom_vline(xintercept = seq(1.5, 30.5, 1), color = 'grey90') +
  geom_hline(yintercept = 0.01, linetype = 'dashed') +
  geom_point(shape = 21, size = 2.5, position = position_dodge2(width = .6, preserve = "single", reverse = TRUE)) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey90", color = "black", inherit.aes = FALSE) +
  geom_text(data = texts, aes(x = x, y = y, label = label),
            angle = 270, size = 10, inherit.aes = FALSE) +
  coord_flip(ylim = c(0, 0.25)) +
  scale_x_discrete(name = NULL, labels = sub("\n", " ", labels), expand = expansion(add = 0.5), limits = rev(unique(relimp_df$var))) +
  scale_y_continuous(expression(paste("LMG (Average Partial ", R^2, ")")), expand = c(0,0), breaks = seq(0, 0.4, 0.1)) +
  scale_fill_manual(name = NULL, breaks = c("Raw Dispersion", "Deviation From Null"),
                    values = c("white", "black")) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.ticks.y = element_blank(),
        axis.line = element_blank(), axis.text = element_text(colour = "black"),
        axis.title.x = element_text(hjust = .3),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position="top", legend.margin = margin(0,120,0,0))
ggsave("Mammal Ranges Relative Importance.pdf", width = 10, height = 14)

#Model Averaging####
# Get average coefficient estimates for the useful parameters from above
## Observed dispersion ####
mass_disp_vars <- relimp_df %>%
  filter(data == "Raw Dispersion", lmg > 0.01) %>%
  arrange(-lmg) %>%
  pull(var)

reg1 <- lm(as.formula(paste0("mass_disp ~ ", paste(mass_disp_vars, collapse = " + "))),
           data = grid_data_sub_clean, na.action = "na.fail")
# standardize coefficients by partial sd to account for intercorrelation of variables
# https://search.r-project.org/CRAN/refmans/MuMIn/html/std.coef.html
drg1 <- dredge(reg1, beta = "partial.sd", trace = 2)
models1 <- model.avg(drg1)
coef1 <- coefTable(models1)
cis1 <- confint(models1)
mod_avg_coefs1 <- cbind.data.frame(term = rownames(coef1), coef1, lower = cis1[,"2.5 %"], upper = cis1[,"97.5 %"])

## Deviation from null ####
mean_diff_vars <- relimp_df %>%
  filter(data == "Deviation From Null", lmg > 0.01) %>%
  arrange(-lmg) %>%
  pull(var)

reg2 <- lm(as.formula(paste0("mean_diff ~ ", paste(mean_diff_vars, collapse = " + "))),
           data = grid_data_sub_cont2_clean, weights = 1/var_diff, na.action = "na.fail")
# standardize coefficients by partial sd to account for intercorrelation of variables
# https://search.r-project.org/CRAN/refmans/MuMIn/html/std.coef.html
drg2 <- dredge(reg2, beta = "partial.sd", trace = 2)
models2 <- model.avg(drg2)
coef2 <- coefTable(models2)
cis2 <- confint(models2)
mod_avg_coefs2 <- cbind.data.frame(term = rownames(coef2), coef2, lower = cis2[,"2.5 %"], upper = cis2[,"97.5 %"])

mod_avg_coefs_all <- rbind(cbind(mod_avg_coefs1, data = "Raw Dispersion"),
                           cbind(mod_avg_coefs2, data = "Deviation From Null"))
mod_avg_coefs_all$data <- factor(mod_avg_coefs_all$data, levels = c("Raw Dispersion", "Deviation From Null"))
mod_avg_coefs_all <- subset(mod_avg_coefs_all, term != "(Intercept)")

## Figure 7 ####
mod_avg_coefs_all$term <- factor(mod_avg_coefs_all$term,
                                 levels = rev(unique(c(grep("continent", mod_avg_coefs_all$term, value = TRUE), "n", "n_xl_xxl",
                                                       grep("mass", mod_avg_coefs_all$term, value = TRUE), "psmall", "pmed", "plarge",
                                                       grep("plant", mod_avg_coefs_all$term, value = TRUE), "pherb", "pomni", "pcarn",
                                                       grep("range", mod_avg_coefs_all$term, value = TRUE),"footprint.2009", "deforestation",
                                                       "elev_mean", str_sort(grep("bio.+_mean", mod_avg_coefs_all$term, value = TRUE), numeric = TRUE),
                                                       "elev_var", str_sort(grep("bio.+_var", mod_avg_coefs_all$term, value = TRUE), numeric = TRUE)
                                 ))))

rects <- data.frame(xmin = c(0.5, 7.5, 18.5), xmax = c(7.5, 18.5, 23.5),
                    ymin = 0.07, ymax = 0.11)

texts <- data.frame(x = c(4, 13, 21), y = 0.09,
                    label = c("Abiotic\nEnvironment", "Biotic Community", "Continent"))

ggplot(data = mod_avg_coefs_all, aes(x = term, group = data, fill = data)) +
  geom_vline(xintercept = seq(1.5, 22.5, 1), color = 'grey90') +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = lower < 0 & upper > 0), width = .6, linewidth = 1, position = position_dodge2(width = .8, preserve = "single", reverse = TRUE)) +
  geom_point(aes(y = Estimate), shape = 21, size = 2.5, position = position_dodge2(width = .6, preserve = "single", reverse = TRUE)) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey90", color = "black", inherit.aes = FALSE) +
  geom_text(data = texts, aes(x = x, y = y, label = label),
            angle = 270, size = 10, inherit.aes = FALSE) +
  coord_flip(ylim = c(-0.07, 0.11)) +
  scale_x_discrete(name = NULL, labels = sub("\n", " ", labels), expand = expansion(add = 0.5),
                   limits = rev(names(labels)[names(labels) %in% mod_avg_coefs_all$term])) +
  scale_y_continuous(name = "Model-Averaged Coefficient Estimate", expand = c(0,0), breaks = seq(-0.05, 0.05, 0.05)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.ticks.y = element_blank(),
        axis.line = element_blank(), axis.text = element_text(colour = "black"),
        axis.title.x = element_text(hjust = -.6),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position="top", legend.margin = margin(0,110,0,0)) +
  scale_fill_manual(name = NULL, breaks = c("Raw Dispersion", "Deviation From Null"),
                    values = c("white", "black")) +
  scale_color_manual(name = NULL, breaks = c("FALSE", "TRUE"), values = c("black", "grey60"), guide = "none")
ggsave("Mammal Ranges Dredge.pdf", width = 10, height = 14)

# Spatial autocorrelation ####
# adapted from https://strimas.com/ebp-workshop/subsampling.html
grid_data_nb <- poly2nb(grid_data_sub_clean$polygon, queen=TRUE)
grid_data_nblist <- nb2listw(grid_data_nb, zero.policy = TRUE)
lm.morantest(reg1, grid_data_nblist) # 0.64
lm.morantest(reg2, grid_data_nblist) # 0.58

dggs2 <- dgconstruct(res = 3)
grid_data_sub_clean$seqnum_big <- dgGEO_to_SEQNUM(dggs2, grid_data_sub_clean$X, grid_data_sub_clean$Y)$seqnum
grid2 <- st_wrap_dateline(dgearthgrid(dggs2))
grid2$centroid <- st_centroid(grid2)$geometry
grid2_keep <- grid2[sort(unique(grid_data_sub_clean$seqnum_big)), ]
table(table(grid_data_sub_clean$seqnum_big))
median(table(grid_data_sub_clean$seqnum_big))

# plot example of spatial subsampling process
ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = st_graticule(crs = "+proj=robin", lon = seq(-180, 180, 30),
                              lat = seq(-90, 90, 30)), color = "gray70") +
  geom_sf(data = robin_without, fill = "white", color = NA) +
  geom_sf(data = robin_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
  geom_sf(data = grid_data_sub_clean, aes(geometry = polygon_robin), fill = NA, color = "grey60") +
  geom_sf(data = grid_data_sub_clean %>% group_by(seqnum_big) %>% sample_n(size = 1),
          aes(geometry = polygon_robin), fill = "black", color = "grey60") +
  geom_sf(data = st_wrap_dateline(grid2_keep, options= c('WRAPDATELINE=YES', 'DATELINEOFFSET=80')),
          color = "grey20", fill = NA) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +
  theme_bw() +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea, color = "white", linewidth = 1),
        panel.grid.major = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.border = element_blank())
ggsave("Mammal Ranges Spatial Subsampling Example.pdf", width = 20, height = 10)

## Relative importance####
num_sim <- 100
reg1_keep <- list()
reg1_moran <- list()
pb <- txtProgressBar(min = 0, max = num_sim, style = 3)
for(i in 1:num_sim){
  setTxtProgressBar(pb, i)
  grid_data_sub_clean_keep <- grid_data_sub_clean %>%
    group_by(seqnum_big) %>%
    sample_n(size = 1) %>%
    ungroup()
  grid_data_sub_clean_keep$polygon_big <- grid2_keep$geometry
  grid_data_keep_nb <- poly2nb(grid_data_sub_clean_keep$polygon_big, queen=TRUE)
  grid_data_keep_nblist <- nb2listw(grid_data_keep_nb, zero.policy = TRUE)
  reg1_keep[[i]] <- lm(mass_disp ~ n + n_xl_xxl + pcarn + pherb + plant_disp + mean_plant + mean_range + range_disp +
                         continent + footprint.2009 + deforestation +
                         bio1_mean + bio4_mean + bio5_mean + bio6_mean + bio12_mean + bio15_mean + elev_mean,
                       #    s(Y, X, bs="sos", k=60),
                       #   bio1_var + bio4_var + bio5_var + bio6_var +
                       #   bio12_var + bio15_var + elev_var,
                       data = grid_data_sub_clean_keep, na.action = "na.omit")
  reg1_moran[[i]] <- lm.morantest(reg1_keep[[i]], grid_data_keep_nblist)
}
close(pb)

summary(sapply(reg1_moran, function(item) item$estimate[1]))

cl <- makeCluster(10)
capture.output(clusterEvalQ(cl, library(relaimpo)), file = nullfile())
spat_sub_sim1 <- parSapply(cl, reg1_keep, function(lst_item) calc.relimp(lst_item, type = "lmg")@lmg) %>%
  as.data.frame()
stopCluster(cl)
spat_sub_sim1$min <- apply(spat_sub_sim1[, 1:num_sim], 1, min)
spat_sub_sim1$max <- apply(spat_sub_sim1[, 1:num_sim], 1, max)
spat_sub_sim1$mean <- apply(spat_sub_sim1[, 1:num_sim], 1, mean)

grid_data_sub_cont2_clean$seqnum_big <- dgGEO_to_SEQNUM(dggs2, grid_data_sub_cont2_clean$X,
                                                        grid_data_sub_cont2_clean$Y)$seqnum
reg2_keep <- list()
reg2_moran <- list()
pb <- txtProgressBar(min = 0, max = num_sim, style = 3)
for(i in 1:num_sim){
  setTxtProgressBar(pb, i)
  grid_data_sub_cont2_clean_keep <- grid_data_sub_cont2_clean %>%
    group_by(seqnum_big) %>%
    sample_n(size = 1) %>%
    ungroup()
  grid_data_sub_cont2_clean_keep$polygon_big <- grid2_keep$geometry
  grid_data_keep_nb <- poly2nb(grid_data_sub_cont2_clean_keep$polygon_big, queen=TRUE)
  grid_data_keep_nblist <- nb2listw(grid_data_keep_nb, zero.policy = TRUE)
  reg2_keep[[i]] <- lm(mean_diff ~ n + n_xl_xxl + psmall + plarge + pcarn + pherb +
                         plant_disp + mean_plant + mean_range + range_disp + continent + footprint.2009 + deforestation +
                         bio1_mean + bio4_mean + bio5_mean + bio6_mean + bio12_mean + bio15_mean + elev_mean,
                       #bio1_var + bio4_var + bio5_var + bio6_var + bio12_var + bio15_var + elev_var,
                       data = grid_data_sub_cont2_clean_keep, weights = 1/var_diff, na.action = "na.omit")
  reg2_moran[[i]] <- lm.morantest(reg2_keep[[i]], grid_data_keep_nblist)
}
close(pb)

summary(sapply(reg2_moran, function(item) item$estimate[1]))

cl <- makeCluster(10)
capture.output(clusterEvalQ(cl, library(relaimpo)), file = nullfile())
spat_sub_sim2 <- parSapply(cl, reg2_keep, function(lst_item) calc.relimp(lst_item, type = "lmg")@lmg) %>%
  as.data.frame()
stopCluster(cl)
spat_sub_sim2$min <- apply(spat_sub_sim2[, 1:num_sim], 1, min)
spat_sub_sim2$max <- apply(spat_sub_sim2[, 1:num_sim], 1, max)
spat_sub_sim2$mean <- apply(spat_sub_sim2[, 1:num_sim], 1, mean)

### Figure 6 ####
relimp_df_sub <- rbind(spat_sub_sim1[, c("min", "mean", "max")] %>%
                     rownames_to_column("var") %>%
                     mutate(data = "Raw Dispersion"),
                   spat_sub_sim2[, c("min", "mean", "max")] %>%
                     rownames_to_column("var") %>%
                     mutate(data = "Deviation From Null")) %>%
  mutate(group = ifelse(grepl("bio.+_mean|elev_mean", var), "habitat",
                        ifelse(grepl("bio.+_var|elev_var", var), "variance", "community"))) %>%
  arrange(group, -mean)

rects <- data.frame(xmin = c(0.5, 7.5), xmax = c(7.5, 20.5),
                    ymin = 0.24, ymax = 0.30)

texts <- data.frame(x = c(4, 14), y = 0.27,
                    label = c("Abiotic\nEnvironment", "Biotic Community"))

ggplot(relimp_df_sub, aes(x = var, y = mean, group = data, fill = data)) +
  geom_vline(xintercept = seq(1.5, 30.5, 1), color = 'grey90') +
  geom_hline(yintercept = 0.01, linetype = 'dashed') +
  geom_errorbar(aes(ymin = min, ymax = max), width = .6, position = position_dodge(width = .8), linewidth = 1) +
  geom_point(shape = 21, size = 2.5, position = position_dodge2(width = .8, preserve = "single")) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey90", color = "black", inherit.aes = FALSE) +
  geom_text(data = texts, aes(x = x, y = y, label = label),
            angle = 270, size = 10, inherit.aes = FALSE) +
  coord_flip(ylim = c(0, 0.30)) +
  scale_x_discrete(name = NULL, labels = sub("\n", " ", labels), expand = expansion(add = 0.5), limits = rev(unique(relimp_df_sub$var))) +
  scale_y_continuous(expression(paste("LMG (Average Partial ", R^2, ")")), expand = c(0,0), breaks = seq(0, 0.4, 0.1)) +
  scale_fill_manual(name = NULL, breaks = c("Raw Dispersion", "Deviation From Null"),
                    values = c("white", "black")) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.ticks.y = element_blank(),
        axis.line = element_blank(), axis.text = element_text(colour = "black"),
        axis.title.x = element_text(hjust = .3),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position="top", legend.margin = margin(0,120,0,0))
ggsave("Mammal Ranges Relative Importance Subsampling.pdf", width = 10, height = 14)

## Model averaging####
### Observed dispersion ####
mass_disp_vars_sub <- relimp_df_sub %>%
  filter(data == "Raw Dispersion", mean > 0.01) %>%
  arrange(-mean) %>%
  pull(var)

cl <- makeCluster(10)
capture.output(clusterEvalQ(cl, {library(MuMIn);library(dplyr)}), file = nullfile())
clusterExport(cl, c("grid_data_sub_clean", "mass_disp_vars_sub"))
drg1_keep <- parLapply(cl, 1:num_sim, function(i) {
  grid_data_sub_clean_keep <- grid_data_sub_clean %>%
    select(-c(geometry, polygon, polygon_igh, polygon_robin)) %>%
    group_by(seqnum_big) %>%
    sample_n(size = 1) %>%
    ungroup()
  reg1_keep <- lm(as.formula(paste0("mass_disp ~ ", paste(mass_disp_vars_sub, collapse = " + "))),
                       data = grid_data_sub_clean_keep, na.action = "na.fail")
  # standardize coefficients by partial sd to account for intercorrelation of variables
  # https://search.r-project.org/CRAN/refmans/MuMIn/html/std.coef.html
  drg_tmp <- dredge(reg1_keep, beta = "partial.sd", trace = 2)
  models_tmp <- model.avg(drg_tmp)
  coef_tmp <- coefTable(models_tmp)
  cis_tmp <- confint(models_tmp)
  cbind.data.frame(term = rownames(coef_tmp), coef_tmp, lower = cis_tmp[,"2.5 %"], upper = cis_tmp[,"97.5 %"])
}) %>% bind_rows(.id = "sim")
stopCluster(cl)

drg1_keep_wtd <- drg1_keep %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  summarise(wtd_mean = Hmisc::wtd.mean(Estimate, w = 1/`Std. Error`^2, normwt=TRUE),
            wtd_sd = sqrt(Hmisc::wtd.var(Estimate, weights = 1/`Std. Error`^2, normwt=TRUE))) %>%
  mutate(lower = wtd_mean - 1.96 * wtd_sd, upper = wtd_mean + 1.96 * wtd_sd)

### Deviation from null ####
mean_diff_vars_sub <- relimp_df_sub %>%
  filter(data == "Deviation From Null", mean > 0.01) %>%
  arrange(-mean) %>%
  pull(var)

cl <- makeCluster(10)
capture.output(clusterEvalQ(cl, {library(MuMIn);library(dplyr)}), file = nullfile())
clusterExport(cl, c("grid_data_sub_cont2_clean", "mean_diff_vars_sub"))
drg2_keep <- parLapply(cl, 1:num_sim, function(i) {
  grid_data_sub_cont2_clean_keep <- grid_data_sub_cont2_clean %>%
    select(-c(geometry, polygon, polygon_igh, polygon_robin)) %>%
    group_by(seqnum_big) %>%
    sample_n(size = 1) %>%
    ungroup()
  reg2_keep <- lm(as.formula(paste0("mean_diff ~ ", paste(mean_diff_vars_sub, collapse = " + "))),
                  data = grid_data_sub_cont2_clean_keep, weights = 1/var_diff, na.action = "na.fail")
  # standardize coefficients by partial sd to account for intercorrelation of variables
  # https://search.r-project.org/CRAN/refmans/MuMIn/html/std.coef.html
  drg_tmp <- dredge(reg2_keep, beta = "partial.sd", trace = 2)
  models_tmp <- model.avg(drg_tmp)
  coef_tmp <- coefTable(models_tmp)
  cis_tmp <- confint(models_tmp)
  cbind.data.frame(term = rownames(coef_tmp), coef_tmp, lower = cis_tmp[,"2.5 %"], upper = cis_tmp[,"97.5 %"])
}) %>% bind_rows(.id = "sim")
stopCluster(cl)

drg2_keep_wtd <- drg2_keep %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  summarise(wtd_mean = Hmisc::wtd.mean(Estimate, w = 1/`Std. Error`^2, normwt=TRUE),
            wtd_sd = sqrt(Hmisc::wtd.var(Estimate, weights = 1/`Std. Error`^2, normwt=TRUE))) %>%
  mutate(lower = wtd_mean - 1.96 * wtd_sd, upper = wtd_mean + 1.96 * wtd_sd)

mod_avg_coefs_all_sub <- rbind(cbind(drg1_keep_wtd, data = "Raw Dispersion"),
                               cbind(drg2_keep_wtd, data = "Deviation From Null"))
mod_avg_coefs_all_sub$data <- factor(mod_avg_coefs_all_sub$data, levels = c("Raw Dispersion", "Deviation From Null"))
### Figure 7 ####
mod_avg_coefs_all_sub$term <- factor(mod_avg_coefs_all_sub$term,
                                 levels = rev(unique(c(grep("continent", mod_avg_coefs_all_sub$term, value = TRUE), "n", "n_xl_xxl",
                                                       grep("mass", mod_avg_coefs_all_sub$term, value = TRUE), "psmall", "pmed", "plarge",
                                                       grep("plant", mod_avg_coefs_all_sub$term, value = TRUE), "pherb", "pomni", "pcarn",
                                                       grep("range", mod_avg_coefs_all_sub$term, value = TRUE),"footprint.2009", "deforestation",
                                                       "elev_mean", str_sort(grep("bio.+_mean", mod_avg_coefs_all_sub$term, value = TRUE), numeric = TRUE),
                                                       "elev_var", str_sort(grep("bio.+_var", mod_avg_coefs_all_sub$term, value = TRUE), numeric = TRUE)
                                 ))))

rects <- data.frame(xmin = c(0.5, 7.5, 17.5), xmax = c(7.5, 17.5, 22.5),
                    ymin = 0.08, ymax = 0.12)

texts <- data.frame(x = c(4, 12.5, 20), y = 0.10,
                    label = c("Abiotic\nEnvironment", "Biotic Community", "Continent"))

ggplot(data = mod_avg_coefs_all_sub, aes(x = term, group = data, fill = data)) +
  geom_vline(xintercept = seq(1.5, 22.5, 1), color = 'grey90') +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = lower < 0 & upper > 0), width = .6, linewidth = 1, position = position_dodge2(width = .8, preserve = "single", reverse = TRUE)) +
  geom_point(aes(y = wtd_mean), shape = 21, size = 2.5, position = position_dodge2(width = .6, preserve = "single", reverse = TRUE)) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey90", color = "black", inherit.aes = FALSE) +
  geom_text(data = texts, aes(x = x, y = y, label = label),
            angle = 270, size = 10, inherit.aes = FALSE) +
  coord_flip(ylim = c(-0.07, 0.12)) +
  scale_x_discrete(name = NULL, labels = sub("\n", " ", labels), expand = expansion(add = 0.5),
                   limits = rev(names(labels)[names(labels) %in% mod_avg_coefs_all_sub$term])) +
  scale_y_continuous(name = "Model-Averaged Coefficient Estimate", expand = c(0,0), breaks = seq(-0.05, 0.05, 0.05)) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.ticks.y = element_blank(),
        axis.line = element_blank(), axis.text = element_text(colour = "black"),
        axis.title.x = element_text(hjust = -.6),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position="top", legend.margin = margin(0,110,0,0)) +
  scale_fill_manual(name = NULL, breaks = c("Raw Dispersion", "Deviation From Null"),
                    values = c("white", "black")) +
  scale_color_manual(name = NULL, breaks = c("FALSE", "TRUE"), values = c("black", "grey60"), guide = "none")
ggsave("Mammal Ranges Dredge Subsampling.pdf", width = 10, height = 14)

# ecosystem engineering? ####
## Figure 8 ####
#xl/xxl species vs null deviation
ggplot(grid_data_sub_cont2_clean, aes(x = n_xl_xxl, y = mean_diff)) +
  geom_violin(aes(group = n_xl_xxl)) +
  geom_quantile(quantiles = c(0.1, 0.5, 0.9), color = "black", linetype = "dashed", linewidth = 1) +
  scale_x_continuous(breaks = 0:20) +
  coord_cartesian(xlim = c(0,16)) +
  xlab("# of Species > 100 kg") +
  ylab(expression("Deviation from the Null (log"[10]*"g)")) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.ticks.y.right = element_blank(),
        axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position=c(.5,.95), legend.margin = margin(0,0,0,0),
        legend.direction = "horizontal")
ggsave("Mammal Ranges Ecosystem Engineering.pdf", width = 10, height = 10)

# plarge vs humans ####
## Figure 8 ####
#plarge vs footprint.2009
ggplot(grid_data_sub_cont2_clean, aes(x = footprint.2009, y = plarge, color = continent)) +
  geom_point(size = 2) +
  geom_quantile(quantiles = c(0.1, 0.5, 0.9), color = "black", linetype = "dashed", linewidth = 1) +
  scale_color_manual(NULL, values = unname(palette.colors()), guide = guide_legend(override.aes = list(size = 3))) +
  xlab("Mean Footprint Index 2009") +
  ylab("Proportion of Large Species") +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.ticks.y.right = element_blank(),
        axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position=c(.5,.95), legend.margin = margin(0,0,0,0),
        legend.direction = "horizontal")
  #guides(color = guide_colourbar(barwidth = unit(15, "lines")))
ggsave("Mammal Ranges Prop Large by Footprint.pdf", width = 10, height = 10)
