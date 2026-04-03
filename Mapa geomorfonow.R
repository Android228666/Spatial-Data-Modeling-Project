library(qgisprocess)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)

# ----------------------------------------------------------
# 1. Ustawienie QGIS 3.22 jako środowiska QGIS dla R
# ----------------------------------------------------------
Sys.setenv(R_QGISPROCESS_PATH = "C:/Program Files/QGIS 3.22.16/bin/qgis_process-qgis-ltr.bat")
qgis_configure(use_cached_data = FALSE)

qgis_enable_plugins(
  c("grassprovider", "processing_saga_nextgen"),
  quiet = TRUE
)

# ----------------------------------------------------------
# 2. Wczytanie danych wejściowych (kolejność: Y, X, Z, ID)
# ----------------------------------------------------------
tab <- read.csv("khomenko.txt", header = FALSE, sep = "")
colnames(tab) <- c("Y", "X", "Z", "ID")   # <-- POPRAWNE NAZWY

# ----------------------------------------------------------
# 3. Konwersja do sf z poprawnym CRS (EPSG:32628)
# ----------------------------------------------------------
tab_sf <- st_as_sf(
  tab,
  coords = c("X", "Y"), 
  crs = 32628       # <-- TWÓJ UKŁAD UTM
)

# ----------------------------------------------------------
# 4. Tworzenie rastra o resolution = 200 m
# ----------------------------------------------------------
raster_template <- rast(
  ext(tab_sf),
  resolution = 600,          # 200 m — OPTYMALNE
  crs = "EPSG:32628"
)

coords_mat <- st_coordinates(tab_sf)

df_raster1 <- rasterize(
  coords_mat,
  raster_template,
  tab$Z
)

plot(df_raster1, main = "DEM rasterized (200 m)")

# 1. Pobierz extents rastra
e <- ext(df_raster1)

xmin <- e$xmin
xmax <- e$xmax
ymin <- e$ymin
ymax <- e$ymax

# 2. Resolution (200m)
res <- 150

# 3. Zapis DEM
writeRaster(df_raster1, "dem_input.tif", overwrite = TRUE, filetype="GTiff")

# 4. Uruchomienie geomorfons Z USTAWIENIEM REGIONU GRASS
dem_geomorph <- qgis_run_algorithm(
  "grass7:r.geomorphon",
  elevation = "dem_input.tif",
  search = 20,
  flat = 0.2,
  forms = "geomorphon_output.tif",
  
  # --- KRYTYCZNE POPRAWKI ---
  GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", xmin, xmax, ymin, ymax),
  GRASS_REGION_CELLSIZE_PARAMETER = res
)

# 5. Wczytaj wynik
geom <- rast("geomorphon_output.tif")
plot(geom)

# ----------------------------------------------------------
# 5. Zapis DEM do pliku dla QGIS/GRASS
# ----------------------------------------------------------
writeRaster(df_raster1, "dem_input.tif", overwrite = TRUE, filetype = "GTiff")


# ----------------------------------------------------------
# 6. Uruchomienie algorytmu GRASS r.geomorphon
# ----------------------------------------------------------
dem_geomorph <- qgis_run_algorithm(
  "grass7:r.geomorphon",
  elevation = "dem_input.tif",
  search = 300,
  flat = 0.01,
  forms = "geomorphon_output.tif"
)

# ----------------------------------------------------------
# 7. Wczytanie wyniku do R
# ----------------------------------------------------------
geom <- rast("geomorphon_output.tif")

plot(geom, main = "Geomorphons (GRASS r.geomorphon)")

# ----------------------------------------------------------
# 8. Zapis wykresu PNG
# ----------------------------------------------------------
png(file = "geomorphons.png", width = 1200, height = 900, res = 150)
plot(geom, main = "Geomorphons (GRASS r.geomorphon)")
dev.off()
