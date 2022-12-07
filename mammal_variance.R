library(sf)
library(terra)
library(grid)
library(rgbif)
library(maptools)
library(spData)
library(deeptime)
library(viridis)
library(ggtern)
library(stringr)
#devtools::install_github("r-barnes/dggridR", vignette=TRUE)
library(dggridR)
library(moments)
library(quantreg)
library(ggforce)
library(relaimpo)
library(MuMIn)
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

colsea = "#00509005"
colland = "#66666660"

#try different cut-offs of area of overlap of grid cells and ranges
#try different grid cell sizes
#remove cosmopolitan species? top 10%?

#IUCN ranges ####
#terrestrial mammals using IUCN ranges
#download from here: https://www.iucnredlist.org/resources/spatial-data-download
#then unzip all of the files into a folder named "MAMMALS_TERRESTRIAL_ONLY"
mammal_IUCN <- st_read(dsn = "MAMMALS_TERRESTRIAL_ONLY", seqnum = "MAMMALS_TERRESTRIAL_ONLY")

#a map of the world
data(wrld_simpl)
wrld_sf <- st_as_sf(wrld_simpl)

#get IUCN species, including synonyms
#both are downloaded from the IUCN website
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
mammal_traits$size_cat <- cut(mammal_traits$BodyMass.log10, c(0,2,4,10), labels = c("small", "medium", "large"), ordered_result = TRUE)

#merge in trait data
mammal_traits_IUCN <- merge(mammal_IUCN_synonyms, mammal_traits, by.x = "IUCN_species", by.y = "Scientific", all.x = TRUE)

#fill in traits based on synonyms if still missing
mammal_traits_IUCN[is.na(mammal_traits_IUCN$Trophic.Level), colnames(mammal_traits)[colnames(mammal_traits) %in% colnames(mammal_traits_IUCN)]] <- mammal_traits[match(mammal_traits_IUCN$species[is.na(mammal_traits_IUCN$Trophic.Level)], mammal_traits$Scientific), colnames(mammal_traits)[colnames(mammal_traits) %in% colnames(mammal_traits_IUCN)]]

#remove duplicates due to synonyms
mammal_traits_IUCN <- mammal_traits_IUCN %>% dplyr::select(-species) %>% unique()
mammal_traits_IUCN <- mammal_traits_IUCN[!is.na(mammal_traits_IUCN$Trophic.Level),]
mammal_traits_IUCN <- mammal_traits_IUCN[order(mammal_traits_IUCN$IUCN_species, mammal_traits_IUCN$Diet.Certainty, -(mammal_traits_IUCN$BodyMass.SpecLevel == 1)),]
mammal_traits_IUCN <- mammal_traits_IUCN[!duplicated(mammal_traits_IUCN$IUCN_species),]

mammal_IUCN_filt <- merge(mammal_IUCN, mammal_traits_IUCN, by.x = "binomial", by.y = "IUCN_species", all.x = TRUE)

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

#Get geographic data####
#human footprint index
#https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-2009-human-footprint
#https://www.nature.com/articles/sdata201667
#https://academic.oup.com/bioscience/article/52/10/891/354831
#~1km resolution
footprint.2009 <- app(rast("wildareas-v3-2009-human-footprint.tif"), fun=function(x) { x[x>50 | x<0] <- NA; return(x) })
footprint.1993 <- app(rast("wildareas-v3-1993-human-footprint.tif"), fun=function(x) { x[x>50 | x<0] <- NA; return(x) })

grid_filt_vect <- vect(grid_filt)
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
wrld_rast <- rasterize(wrld_sf, deforestation, fun=function(x, ...){0})
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
mammal_IUCN_clean <- st_make_valid(mammal_IUCN_filt)
mammal_IUCN_clean$SHAPE_Area <- st_area(mammal_IUCN_clean)
mammal_ranges <- mammal_IUCN_clean %>%
  st_drop_geometry() %>%
  group_by(binomial) %>%
  summarize(range_size = sum(SHAPE_Area), n_shapes = length(SHAPE_Area))
sf_use_s2(TRUE)
mammal_IUCN_clean <- left_join(mammal_IUCN_clean, mammal_ranges, by = "binomial")

# find which species intersect with which grid cells
# this can require a bit of time and RAM
intrsct <- terra::intersect(grid_filt_vect, vect(mammal_IUCN_clean))

# Extract areas from polygon objects then attach as attribute
# area in square meters
mammal_grids <- as.data.frame(intrsct)
mammal_grids$area <- expanse(intrsct)
mammal_grids$BodyMass.log10 <- log10(mammal_grids$BodyMass.Value)
mammal_grids$range_size.log10 <- log10(as.numeric(mammal_grids$range_size))

#determine continent of each seqnum
conts <- world %>%
  st_wrap_dateline() %>%
  st_make_valid() %>%
  select(continent) %>%
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
                        summarize(continent = continent[which.max(area)]),
                      by = "seqnum", all.x = TRUE)

#Order IUCN categories
mammal_grids$category <- factor(mammal_grids$category, ordered = TRUE,
                                levels = c("DD","LC","NT","VU","EN","CR","EW","EX"))

#Grid cell data####
#calculate number of species overlapping with grid cell, size stats, diet stats, range stats
grid_data <- mammal_grids %>% group_by(seqnum, continent) %>%
  summarize(n = n(),
            mean_mass = mean(BodyMass.log10, na.rm = TRUE), median_mass = median(BodyMass.log10, na.rm = TRUE),
            min_mass = min(BodyMass.log10, na.rm = TRUE), max_mass = max(BodyMass.log10, na.rm = TRUE),
            n_mass = sum(!is.na(BodyMass.Value)), var_mass = var(BodyMass.log10, na.rm = TRUE),
            mass_disp = mean(dist(BodyMass.log10, method = "manhattan"), na.rm = TRUE), #body mass dispersion https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12667
            skew_mass = skewness(BodyMass.log10, na.rm = TRUE), kurt_mass = kurtosis(BodyMass.log10, na.rm = TRUE),
            n_small = sum(BodyMass.log10 < 2, na.rm = TRUE), n_large = sum(BodyMass.log10 >= 4, na.rm = TRUE),
            n_med = sum(BodyMass.log10 >= 2 & BodyMass.log10 < 4, na.rm = TRUE),
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

#Figure 2####
lmsumm1 <- summary(lm(mass_disp~n, data = grid_data_sub))
lmsumm2 <- summary(lm(kurt_mass~n, data = grid_data_sub))
lmsumm3 <- summary(lm(skew_mass~n, data = grid_data_sub))

g1 <- ggplot(grid_data_sub, aes(x = n, y = mass_disp)) +
  geom_point(shape = 21, fill = "grey80") +
  geom_smooth(method = lm, se = FALSE, color = "black", linetype = "dashed") +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = "Mass Dispersion (log g)") +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA)) +
  annotate(geom = "text", x = 320, y = 2, hjust = 1, size = 10,
           label = deparse(bquote(p == .(round(lmsumm1$coefficients[2,4], 3))*","~ R^2 == .(round(lmsumm1$r.squared, 3)))), parse = T)
g2 <- ggplot(grid_data_sub, aes(x = n, y = kurt_mass)) +
  geom_point(shape = 21, fill = "grey80") +
  geom_smooth(method = lm, se = FALSE, color = "black", linetype = "dashed") +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = "Mass Kurtosis") +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA)) +
  annotate(geom = "text", x = 320, y = 5.7, hjust = 1, size = 10,
           label = deparse(bquote(p == .(round(lmsumm2$coefficients[2,4], 3))*","~ R^2 == .(round(lmsumm2$r.squared, 3)))), parse = T)
g3 <- ggplot(grid_data_sub, aes(x = n, y = skew_mass)) +
  geom_point(shape = 21, fill = "grey80") +
  geom_smooth(method = lm, se = FALSE, color = "black", linetype = "dashed") +
  scale_x_continuous(name = "Community Species Richness") +
  scale_y_continuous(name = "Mass Skewness") +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA)) +
  annotate(geom = "text", x = 320, y = 1.8, hjust = 1, size = 10,
           label = deparse(bquote(p == .(round(lmsumm3$coefficients[2,4], 3))*","~ R^2 == .(round(lmsumm3$r.squared, 3)))), parse = T)
gg <- ggarrange2(g1, g2, g3, ncol = 1, labels = c("(a)", "(b)", "(c)"),
                 label.args = list(gp = grid::gpar(font = 2, cex = 2)),  draw = FALSE)
ggsave("Mammal Ranges Stats by N.pdf", gg, width = 10, height = 15)

#Figure 1####
#mass moments and species richness
g1 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = mass_disp)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Mass Dispersion") +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = FALSE)
g1_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(mass_disp), fill = viridis(30, option = "plasma")) +
  labs(x = expression("Mass Dispersion (log"[10]*"g)"), y = "# of Communities") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
g2 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = kurt_mass)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Mass Kurtosis") +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = FALSE)
g2_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(kurt_mass), fill = viridis(30, option = "plasma")) +
  labs(x = "Mass Kurtosis", y = "# of Communities") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
g3 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = skew_mass)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Mass Skewness") +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = FALSE)
g3_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(skew_mass), fill = viridis(30, option = "plasma")) +
  labs(x = "Mass Skewness", y = "# of Communities") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
g4 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub, aes(geometry = geometry, color = n)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Species Richness") +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = FALSE)
g4_hist <- ggplot(data = grid_data_sub) +
  geom_histogram(aes(n), fill = viridis(30, option = "plasma")) +
  labs(x = "Species Richness", y = "# of Communities") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
gg <- ggarrange2(g4,g4_hist,g1,g1_hist,g2,g2_hist,g3,g3_hist, ncol = 2, widths = c(2,1),
                 labels = c("(a)", "", "(b)", "", "(c)", "", "(d)", ""),
                 label.args = list(gp = grid::gpar(font = 2, cex = 2.5)),  draw = FALSE)
ggsave("Mammal Ranges Mass Moments Maps.pdf", gg, width = 17.5, height = 25)

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
  scale_color_viridis(option = "plasma", guide = FALSE)
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
  scale_color_viridis(option = "plasma", guide = FALSE)
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
  scale_color_viridis(option = "plasma", guide = FALSE)
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
  scale_color_viridis(option = "plasma", guide = FALSE)
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
  scale_color_viridis(option = "plasma", guide = FALSE)
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
  scale_color_viridis(option = "plasma", guide = FALSE)
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
  scale_color_viridis(option = "plasma", guide = FALSE)
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
  scale_color_viridis(option = "plasma", guide = FALSE)
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
  scale_color_viridis(option = "plasma", guide = FALSE)
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
  scale_color_viridis(option = "plasma", guide = FALSE)
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
  scale_color_viridis(option = "plasma", guide = FALSE)
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
  scale_color_viridis(option = "plasma", guide = FALSE)
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
  geom_arc_bar(data = diet_gather, aes(x0 = x, y0 = y, r0 = 0, r = 1.1, amount = value,
                                       group = seqnum, fill = type), stat='pie') +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Dietary Composition of Communities") +
  theme(axis.text = element_text(colour = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_fill_viridis(name = "Diet", discrete = TRUE, guide = FALSE)
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
size_gather <- gather(grid_data_sub, "type", "value", c("n_small","n_med","n_large"))
size_gather$type <- factor(size_gather$type, levels = c("n_small","n_med","n_large"))
g2 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_arc_bar(data = size_gather, aes(x0 = X, y0 = Y, r0 = 0, r = 1.1, amount = value,
                                       group = seqnum, fill = type), stat='pie') +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Size Composition of Communities") +
  theme(axis.text = element_text(colour = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_fill_viridis(name = "Diet", discrete = TRUE, guide = FALSE, option = "plasma")
psize_gather <- gather(grid_data_sub, "type", "value", c("psmall","pmed","plarge"))
psize_gather$type <- factor(psize_gather$type, levels = c("psmall","pmed","plarge"))
g2_hist <- ggplot(psize_gather) +
  geom_histogram(aes(value * 100, fill = type), binwidth = 2.5, boundary = 25, show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, option = "plasma") +
  scale_x_continuous(name = NULL, limits = c(0, 100)) +
  scale_y_continuous(name = "# of Communities") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  facet_wrap(~type, ncol = 1, strip.position = "right", scales = "free_x",
             labeller = as_labeller(c("psmall" = "% <100g","pmed" = "% 100g - 10kg", "plarge" = "% >10kg")))
gg <- ggarrange2(g2, g2_hist, nrow = 1, widths = c(7,1), draw = FALSE)
ggsave("Mammal Ranges Size Pie Chart Map.pdf", gg, width = 30, height = 13.5)

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

sim_results_cont <- cbind(data.frame(sim = rep(seq(1:num_sim), each = nrow(grid_data_sub)), sim_disp = as.vector(sim_data_cont), diff = as.vector(grid_data_sub$mass_disp - sim_data_cont)),
                          grid_data_sub[rep(seq_len(nrow(grid_data_sub)), num_sim), ])

grid_data_sub_cont <- subset(sim_results_cont, !is.na(mass_disp)) %>%
  group_by(across(c(-sim, -sim_disp, -diff))) %>%
  summarize(mean_sim_disp = mean(sim_disp, na.rm = TRUE), sd_sim_disp = sd(sim_disp, na.rm = TRUE),
            mean_diff = mean(diff, na.rm = TRUE),
            var_diff = var(diff, na.rm = TRUE),
            p_less = wilcox.test(sim_disp, mu = unique(mass_disp), alternative = "less")$p.value,
            p_greater = wilcox.test(sim_disp, mu = unique(mass_disp), alternative = "greater")$p.value)
grid_data_sub_cont$p_less_adjust <- p.adjust(grid_data_sub_cont$p_less, method = "fdr", n = nrow(grid_data_sub_cont) * 2)
grid_data_sub_cont$p_greater_adjust <- p.adjust(grid_data_sub_cont$p_greater, method = "fdr", n = nrow(grid_data_sub_cont) * 2)

table(cut(grid_data_sub_cont$p_less_adjust, c(0,0.05,1)), exclude = NULL)/nrow(grid_data_sub_cont)

#Figure 5####
g1 <- ggplot(grid_data_sub_cont) +
  geom_errorbar(aes(x = n, ymin = mean_sim_disp - 2 * sd_sim_disp,
                    ymax = mean_sim_disp + 2 * sd_sim_disp), width = 0, color = "grey70") +
  geom_point(aes(x = n, y = mean_sim_disp), shape = 21, fill = "grey100", show.legend = TRUE) +
  geom_point(aes(x = n, y = mass_disp, fill = cut(plarge, quantile(plarge, c(0, .5, 1)), include.lowest = TRUE)), shape = 21) +
  coord_cartesian(ylim = c(.5, 2.25), xlim = c(0, max(grid_data_sub_cont$n) + 2)) +
  scale_x_continuous(name = "Community Species Richness", expand = c(0.001, 0.001)) +
  scale_y_continuous(name = "Mass Dispersion (log g)") +
  scale_fill_manual(name = NULL, labels = c("< 17.5% Large Species", "> 17.5% Large Species", "Null Model"),
                    values = viridis_pal()(5)[c(2,4)]) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position = c(.775, .95), legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0), strip.background = element_blank(),
        legend.background = element_rect(fill = NA)) +
  guides(fill = guide_legend(override.aes = list("size" = 3)))
lmsumm <- summary(lm(mean_diff~plarge, data = grid_data_sub_cont))
g2 <- ggplot(grid_data_sub_cont) +
  geom_point(aes(x = plarge, y = mean_diff, fill = cut(plarge, quantile(plarge, c(0, .5, 1)), include.lowest = TRUE)), shape = 21) +
  geom_smooth(aes(x = plarge, y = mean_diff), method = lm, se = FALSE, color = "black", linetype = "dashed") +
  annotate(geom = "text", x = .7, y = -.4, hjust = 1, size = 10,
           label = deparse(bquote(p == .(round(lmsumm$coefficients[2,4], 3))*","~ R^2 == .(round(lmsumm$r.squared, 3)))), parse = T) +
  scale_x_continuous(name = "Community Proportion of Large Species (>10kg)") +
  scale_y_continuous(name = "Mass Dispersion Deviation (log g)") +
  scale_fill_manual(name = NULL, labels = c("< 17.5% Large Species", "> 17.5% Large Species"),
                    values = viridis_pal()(5)[c(2,4)]) +
  theme_classic(base_size = 24) +
  theme(axis.ticks = element_line(color = "black"), axis.line = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"), panel.border = element_rect(fill = NA),
        legend.position = c(.225, .95), legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0), strip.background = element_blank(),
        legend.background = element_rect(fill = NA)) +
  guides(fill = guide_legend(override.aes = list("size" = 3)))
gg <- ggarrange2(g1, g2, ncol = 2, widths = c(1,1), labels = c("(a)", "(b)"),
                 label.args = list(gp = grid::gpar(font = 2, cex = 2.5)),  draw = FALSE)
ggsave("Mammal Ranges Simulation Variance - Ranges and Diet Maintained with plarge.pdf", gg, width = 20, height = 10)

#Figure 4####
#plot simulation value and simulation difference on map
g1 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub_cont, aes(geometry = geometry, color = mean_sim_disp)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Null Model Mean Mass Dispersion") +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = FALSE)
g1_hist <- ggplot(data = grid_data_sub_cont) +
  geom_histogram(aes(mean_sim_disp), fill = viridis(30, option = "plasma")) +
  labs(x = "Null Mass Dispersion (log g)", y = "# of Communities") +
  coord_cartesian(expand =FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))

g2 <- ggplot() +
  geom_sf(data = land, fill = colland) +
  geom_sf(data = grid_data_sub_cont, aes(geometry = geometry, color = mean_diff)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90)) +
  scale_x_continuous(expand = c(0,0), breaks = c(seq(-180,180,30)), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0,0), breaks = c(seq(-90,90,30)), minor_breaks = NULL) +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = NULL, title = "Observed Mass Dispersion - Null Model Mean Mass Dispersion") +
  theme(axis.text = element_text(color = "black")) +
  theme(panel.background = element_rect(fill = colsea), plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  scale_color_viridis(option = "plasma", guide = FALSE)
g2_hist <- ggplot(data = grid_data_sub_cont) +
  geom_histogram(aes(mean_diff), fill = viridis(30, option = "plasma")) +
  labs(x = "Deviation from the Null (log g)", y = "# of Communities") +
  coord_cartesian(expand =FALSE) +
  theme_classic(base_size = 25) +
  theme(axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(.5,.5,.5,.5), "cm"))
gg <- ggarrange2(g1, g1_hist, g2, g2_hist, ncol = 2, widths = c(2,1), labels = c("(a)", "", "(b)", ""),
                 label.args = list(gp = grid::gpar(font = 2, cex = 2.5)),  draw = FALSE)
ggsave("Mammal Ranges Simulation Values and Deviations - Ranges and Diet Maintained.pdf", gg, width = 20, height = 15)

# Relative Importance ####
# What are the best predictors of variance in general?
labels <- c("continent" = "Continent",
            "continentAsia" = "Asia", "continentOceania" = "Oceania", "continentEurope" = "Europe",
            "continentNorth America" = "North America", "continentSouth America" = "South America",
            "n" = "Species Richness",
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
grid_data_sub_clean <- subset(grid_data_sub, !is.na(mass_disp) & !is.na(elev_mean))
grid_data_sub_clean$continent <- factor(grid_data_sub_clean$continent)
reg1 <- lm(mass_disp ~ n + pcarn + pherb + plant_disp + mean_plant + mean_range + range_disp +
             continent + footprint.2009 + deforestation + bio1_mean + bio4_mean + bio5_mean +
             bio6_mean + bio12_mean + bio15_mean + bio1_var + bio4_var + bio5_var + bio6_var +
             bio12_var + bio15_var + elev_mean + elev_var,
           data = grid_data_sub_clean, na.action = "na.fail")
relimp1 <- calc.relimp(reg1, type = "lmg")
sort(relimp1@lmg, decreasing = TRUE)

## Deviation from null ####
# remove pmed and pomni to prevent singularity issues
grid_data_sub_cont_clean <- subset(grid_data_sub_cont, !is.na(mean_diff) & !is.na(elev_mean))
grid_data_sub_cont_clean$continent <- factor(grid_data_sub_cont_clean$continent)
reg2 <- lm(mean_diff ~ n + max_mass + min_mass + psmall + plarge + pcarn + pherb +
             plant_disp + mean_plant + mean_range + range_disp + continent + footprint.2009 + deforestation +
             bio1_mean + bio4_mean + bio5_mean + bio6_mean + bio12_mean + bio15_mean +
             bio1_var + bio4_var + bio5_var + bio6_var + bio12_var + bio15_var + elev_mean + elev_var,
           data = grid_data_sub_cont_clean, weights = 1/var_diff, na.action = "na.fail")
relimp2 <- calc.relimp(reg2, type = "lmg")
sort(relimp2@lmg, decreasing = TRUE)

## Figure 6 ####
relimp_df <- rbind(data.frame(var = relimp1@namen[2:length(relimp1@namen)],
                              data = "Raw Dispersion", lmg = relimp1@lmg),
                   data.frame(var = relimp2@namen[2:length(relimp2@namen)],
                              data = "Deviation From Null", lmg = relimp2@lmg)) %>%
  mutate(group = ifelse(grepl("bio.+_mean|elev_mean", var), "habitat", ifelse(grepl("bio.+_var|elev_var", var), "variance", "community"))) %>%
  arrange(data, group, -lmg)

rects <- data.frame(xmin = c(0.5, 7.5, 14.5), xmax = c(7.5, 14.5, 28.5),
                    ymin = .4, ymax = .7)

texts <- data.frame(x = c(4, 11, 21.5), y = .45,
                    label = c("Habitat\nVariance", "Habitat\nAverage", "Community"))

# without bootstrap limits
ggplot(relimp_df, aes(x = var, y = lmg, group = data, fill = data)) +
  geom_vline(xintercept = seq(1.5, 30.5, 1), color = 'grey90') +
  geom_hline(yintercept = 0.01, linetype = 'dashed') +
  geom_point(shape = 21, size = 2.5, position = position_dodge2(width = .6, preserve = "single", reverse = TRUE)) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey90", color = "black", inherit.aes = FALSE) +
  geom_text(data = texts, aes(x = x, y = y, label = label),
            angle = 270, size = 10, inherit.aes = FALSE) +
  coord_flip(ylim = c(0, 0.5)) +
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
ggsave("Mammal Ranges Relative Importance.pdf", width = 11, height = 14)

#Model Averaging####
# Get average coefficient estimates for the useful parameters from above
## Observed dispersion ####
mass_disp_vars <- relimp_df %>%
  filter(data == "Raw Dispersion", lmg > 0.01) %>%
  arrange(-lmg) %>%
  pull(var)

reg1 <- lm(as.formula(paste0("mass_disp ~ ", paste(mass_disp_vars, collapse = " + "))),
           data = grid_data_sub_clean, na.action = "na.fail")
drg1 <- dredge(reg1, beta = "sd", trace = 2)
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
           data = grid_data_sub_cont_clean, weights = 1/var_diff, na.action = "na.fail")
drg2 <- dredge(reg2, beta = "sd", trace = 2)
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
                                 levels = rev(unique(c(grep("continent", mod_avg_coefs_all$term, value = TRUE), "n",
                                                       grep("mass", mod_avg_coefs_all$term, value = TRUE), "psmall", "pmed", "plarge",
                                                       grep("plant", mod_avg_coefs_all$term, value = TRUE), "pherb", "pomni", "pcarn",
                                                       grep("range", mod_avg_coefs_all$term, value = TRUE),"footprint.2009", "deforestation",
                                                       "elev_mean", str_sort(grep("bio.+_mean", mod_avg_coefs_all$term, value = TRUE), numeric = TRUE),
                                                       "elev_var", str_sort(grep("bio.+_var", mod_avg_coefs_all$term, value = TRUE), numeric = TRUE)
                                 ))))

rects <- data.frame(xmin = c(0.5, 5.5, 14.5), xmax = c(5.5, 14.5, 19.5),
                    ymin = 1.45, ymax = 2.8)

texts <- data.frame(x = c(3, 10, 17), y = 1.85,
                    label = c("Habitat\nAverage", "Community", "Continent"))

ggplot(data = mod_avg_coefs_all, aes(x = term, group = data, fill = data)) +
  geom_vline(xintercept = seq(1.5, 22.5, 1), color = 'grey90') +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = lower < 0 & upper > 0), width = .6, linewidth = 1, position = position_dodge2(width = .8, preserve = "single", reverse = TRUE)) +
  geom_point(aes(y = Estimate), shape = 21, size = 2.5, position = position_dodge2(width = .6, preserve = "single", reverse = TRUE)) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey90", color = "black", inherit.aes = FALSE) +
  geom_text(data = texts, aes(x = x, y = y, label = label),
            angle = 270, size = 10, inherit.aes = FALSE) +
  coord_flip(ylim = c(-1.25,2)) +
  scale_x_discrete(name = NULL, labels = sub("\n", " ", labels), expand = expansion(add = 0.5),
                   limits = rev(names(labels)[names(labels) %in% mod_avg_coefs_all$term])) +
  scale_y_continuous(name = "Model-Averaged Coefficient Estimate", breaks = seq(-1, 1, 1)) +
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

# plarge vs humans ####
## Figure 8 ####
#plarge vs footprint.2009
ggplot(grid_data_sub_cont_clean, aes(x = footprint.2009, y = plarge, color = continent)) +
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
