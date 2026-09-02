library(terra)
library(ggplot2)
library(tidyterra)

Spat <- terra::rast("Predictors/gudenaa_landuse/Landuse_gudenaa_current.tif")

levels(Spat)


# 1) Put Exclude first, then Forest*, then Open*
ord <- c(
  "ForestDryPoor", "ForestWetPoor", "ForestDryRich", "ForestWetRich",
  "OpenDryPoor",  "OpenWetPoor",   "OpenDryRich",   "OpenWetRich", "Agriculture", "Urban", "ProductionForest","Other"
)



lev <- levels(Spat)[[1]]
lev$class <- factor(lev$class, levels = ord)
levels(Spat) <- lev

# 2) Use a *named* palette so colors map to labels reliably
pal <- c(
  ForestDryPoor = "#01665e",
  ForestWetPoor = "#35978f",
  ForestDryRich = "#80cdc1",
  ForestWetRich = "#c7eae5",
  OpenDryPoor   = "#f6e8c3",
  OpenWetPoor   = "#dfc27d",
  OpenDryRich   = "#bf812d",
  OpenWetRich   = "#8c510a",
  Agriculture = "#CC6677",
  Urban = "grey",
  ProductionForest = "#80739B",
  Other       = "white"
)


shp <- terra::vect("Datasets/Gudenåen_hovedløb/HovedlobMarts24.shp")
Spat <- terra::project(Spat, terra::crs(shp), "near")

ggplot() +
  geom_spatraster(data = Spat, maxcell = 4000000) +
   geom_spatvector(data = shp, color = "blue", linewidth = 1, fill = NA) +
  scale_fill_manual(
    name = "Current",
    values = pal,
    breaks = ord,          # legend order
    drop = FALSE,
    na.value = "#00000000"
  ) +
  theme_dark()


Spat <- terra::mask(Spat, shp)
Spat <- terra::crop(Spat, shp)

ggplot() +
  geom_spatraster(data = Spat, maxcell = 4000000) +
  scale_fill_manual(
    name = "Current",
    values = pal,
    breaks = ord,          # legend order
    drop = FALSE,
    na.value = "#00000000"
  ) +
  theme_dark()
