library(here) # hilft beim organisieren der projektstruktur
library(terra) # zum bearbeiten räumlicher daten
library(tidyverse) # datenverarbeitung und visualisierung
library(rgbif) # herkunft von artenverbreitungsdaten
library(geodata) # herkunft von umweltdaten
library(ggspatial) # visualisierung räumlicher daten
Europe_ext <- c(-15,45,35,72)
Europe_map <- crop(x = world(resolution = 2,
                             level = 0),
                   y = Europe_ext)
plot(Europe_map, main = "Europe")
bioclim_data <- worldclim_global(var = "bio",
                                 res = 5)
plot(bioclim_data)
bioclim_data <- crop(bioclim_data, Europe_ext)
plot(bioclim_data)
gbif_ArnMon_download <- occ_data(scientificName = "Arnica montana",
                                 hasCoordinate = TRUE,
                                 limit = 1000,
                                 year="1999,2005",
                                 basisOfRecord = "HUMAN_OBSERVATION")
gbif_ArnMon <- gbif_ArnMon_download$data
head(gbif_ArnMon)
gbif_ArnMon <- gbif_ArnMon[,c("key",
                              "decimalLatitude",
                              "decimalLongitude",
                              "occurrenceStatus",
                              "coordinateUncertaintyInMeters")]
summary(gbif_ArnMon)
plot(Europe_map,
     axes = TRUE,
     col = "grey95",
     main = "Vorkommen in Europa")
points(x = gbif_ArnMon$decimalLongitude,
       y = gbif_ArnMon$decimalLatitude,
       col = "orange",
       pch = 16,
       cex = 0.5)
ggplot(gbif_ArnMon, aes(x = coordinateUncertaintyInMeters)) +
  geom_histogram(bins = 1e2) +
  theme_bw() +
  scale_y_continuous(trans = "log10")
gbif_ArnMon <- gbif_ArnMon[which(gbif_ArnMon$coordinateUncertaintyInMeters < 100), ]
background <- spatSample(x = bioclim_data,
                         size = nrow(gbif_ArnMon), # gleiche anzahl wie präsenzen
                         values = FALSE, # werte braucht es nicht
                         na.rm = TRUE, # nicht von na-zellen samplen
                         xy = TRUE) # mit koordinaten
head(background)
plot(Europe_map,
     axes = TRUE,
     col = "grey95",
     main = "Präsenzen und Absenzen in Europa")
points(background,
       col = "darkgrey",
       pch = 1,
       cex = 0.75)
points(x = gbif_ArnMon$decimalLongitude,
       y = gbif_ArnMon$decimalLatitude,
       col = "orange",
       pch = 16,
       cex = 0.5)
ArnMon_presences <- gbif_ArnMon[, c("decimalLongitude", "decimalLatitude")]
colnames(ArnMon_presences) <- c("longitude", "latitude")
ArnMon_presences$pa <- 1
ArnMon_absences <- as.data.frame(background)
colnames(ArnMon_absences) <- c("longitude", "latitude")
ArnMon_absences$pa <- 0
ArnMon_PA <- rbind(ArnMon_presences, ArnMon_absences)
bioclim_df <- terra::extract(x = bioclim_data,
                             y = ArnMon_PA[, c("longitude", "latitude")],
                             ID = FALSE)
ArnMon_complete_df <- cbind(ArnMon_PA, bioclim_df)
head(ArnMon_complete_df)
glm_ArnMon <- glm(pa ~ .,
                  data = ArnMon_complete_df[,-c(1,2)],
                  family = binomial())
predict_ArnMon <- predict(bioclim_data, glm_ArnMon, type = "response")
plot(predict_ArnMon)
forecast_data <- cmip6_world(model = "MPI-ESM1-2-HR",
                             ssp = "245",
                             time = "2061-2080",
                             var = "bioc",
                             res = 5,
                             path = here("data"))
names(forecast_data) <- names(bioclim_data)
forecast_data <- crop(x = forecast_data, y = Europe_ext)
forecast_presence <- predict(forecast_data, glm_ArnMon, type = "response")
plot(Europe_map,
     axes = TRUE,
     col = "grey95")
plot(forecast_presence, add = TRUE)
plot(Europe_map, add = TRUE, border = "grey5")
points(x = ArnMon_presences$longitude,
       y = ArnMon_presences$latitude,
       col = "black",
       pch = "+",
       cex = 0.75)
forecast_data_585 <- cmip6_world(
  model = "MPI-ESM1-2-HR",
  ssp = "585",
  time = "2061-2080",
  var = "bioc",
  res = 5,
  path = here("data")
)
names(forecast_data_585) <- names(bioclim_data)
forecast_data_585 <- crop(x = forecast_data_585, y = Europe_ext)
forecast_presence_585 <- predict(forecast_data_585, glm_ArnMon, type = "response")
plot(Europe_map,
     axes = TRUE,
     col = "grey95")
plot(forecast_presence_585, add = TRUE)
plot(Europe_map, add = TRUE, border = "grey5")
points(x = ArnMon_presences$longitude,
       y = ArnMon_presences$latitude,
       col = "black",
       pch = "+",
       cex = 0.75)
png("current_distribution.png", width = 1200, height = 900, res = 150)

predict_ArnMon_df <- as.data.frame(predict_ArnMon, xy = TRUE)

ggplot() +
  geom_raster(data = predict_ArnMon_df,
              aes(x = x, y = y, fill = lyr1)) +
  annotation_spatial(Europe_map, fill = NA, colour = "black") +
  scale_fill_viridis_c() +
  ggtitle("Potenzielle Verbreitung von Arnica montana (aktuell)") +
  theme_bw() +
  annotation_north_arrow(location = "tr", which_north = "true") +
  annotation_scale(location = "tl")

dev.off()
library(ggplot2)
library(ggspatial)
library(sf)

df_245 <- as.data.frame(forecast_presence, xy = TRUE)
colnames(df_245)[3] <- "value"

Europe_sf <- st_as_sf(Europe_map)

png("ssp245_distribution.png", width = 1200, height = 900, res = 150)

ggplot() +
  geom_raster(data = df_245,
              aes(x = x, y = y, fill = value)) +
  geom_sf(data = Europe_sf,
          fill = NA,
          color = "black",
          linewidth = 0.3) +
  scale_fill_viridis_c(name = "Eignung") +
  ggtitle("Potenzielle Verbreitung von Arnica montana (SSP245, 2061–2080)") +
  coord_sf() +
  annotation_north_arrow(location = "tr", which_north = "true") +
  annotation_scale(location = "tl") +
  theme_bw()

dev.off()
forecast_presence_585
plot(forecast_presence_585)
df_585 <- as.data.frame(forecast_presence_585, xy = TRUE)
colnames(df_585)[3] <- "value"
Europe_sf <- st_as_sf(Europe_map)

p585 <- ggplot() +
  geom_tile(data = df_585,
            aes(x = x, y = y, fill = value)) +
  geom_sf(data = Europe_sf,
          fill = NA,
          color = "black",
          linewidth = 0.3) +
  scale_fill_viridis_c(name = "Eignung") +
  ggtitle("Potenzielle Verbreitung von Arnica montana (SSP585, 2061–2080)") +
  coord_sf() +
  annotation_north_arrow(location = "tr", which_north = "true") +
  annotation_scale(location = "tl") +
  theme_bw()

print(p585)
ggsave("ssp585_distribution.png", plot = p585, width = 8, height = 6, dpi = 150)












