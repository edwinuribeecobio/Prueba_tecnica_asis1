# Al ejecutar por primera se debe activar la instalacion de paquetes 
#packages <- c("prioritizr", "prioritizrdata", "terra", "sf", "raster", #"fasterize", "ggplot2", "viridis", "ggspatial")

# Instalar los paquetes que no estan ya instalados
#install.packages(setdiff(packages, rownames(installed.packages())))

# install.packages("highs", repos = "https://cran.rstudio.com/") 

# Cargar librerias necesarias
library(prioritizr)
library(prioritizrdata)
library(terra)
library(sf)
library(raster)
library(fasterize)
library(ggplot2)
library(viridis)
library(ggspatial)

# Opcionales
# Quitar notacion cientifica
options(scipen=999)


# Establecer directorio de trabajo
setwd('D:/Convocatorias/Humboldt/Prueba_Tecnica_Asistente1')
# Insumos
costo <- raster('Capa_costos/RASTER/Beneficio_Neto_Total.tif')
model.distribucion = raster('Tremarctos ornatus/Tremarctos ornatus.tif')
hh <- raster('Huella Humana/IHEH_2018.tif')
paramos <- st_read('Paramos/Complejos de Paramos_Escala100k.shp')
runap <- st_read("D:/TNC/Biofisica_INPUTS/RUNAP.shp")



# Establecer area de estudio (ae) ----
ae <- model.distribucion
ae <- projectRaster(ae, crs = 32618, method = 'ngb')
ae.poly <- as.polygons(rast(ae), dissolve=TRUE)
ae.poly <- buffer(ae.poly, width=10000) 
ae <- crop(ae, extent(st_as_sf(ae.poly)))
# Crear plantilla
template <- raster(resolution = 4000, crs= crs(ae),ext = extent(ae))



raster_df <- as.data.frame(ae, xy = TRUE)
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = Tremarctos.ornatus)) +
  geom_sf(data = st_as_sf(ae.poly), color = "black", fill = NA, size = 0.0001,  alpha = 0.5) + 
  scale_fill_viridis_c(na.value = NA) +  
  theme_grey() +  guides(fill = "none") +
  labs(title = "Area de estudio")



# Homogenizar insumos ----

# Rasterizar insumos SHP
paramos <- fasterize(paramos, template)
runap <- st_crop(st_transform(runap, crs = st_crs(ae)), ae)
runap <- fasterize(runap, template)

# Iteraracion sobre insumos
#Lista
tiffs = list(costo = costo, model.distribucion = model.distribucion, hh = hh, paramos = paramos, runap = runap)
# Loop
for (i in 1:length(tiffs)){
  r <- tiffs[[i]]
  r1 <- projectRaster(r, crs = crs(ae), method = "ngb")
  r1 <- crop(r1, extent(ae))
  r1 <- resample(rast(r1), rast(template), method = "near")
  # Las capas resultantes se guardan con el sufijo ".h"
  assign(paste0(names(tiffs)[i], '.h'), r1)
  message(paste('Insumo preparado:', names(tiffs)[i]))
}



# Carpinteria específica por capas
# Costos
costo.h.t = costo.h/1000000000000
costo.h.t = mask(costo.h.t, ae.poly)
# Modelo de distribucion
# remplazar con ceros donde la especie no este presente, 
# dejar NA todo lo que esta por fuera del area de estudio.
model.distribucion.h[is.na(model.distribucion.h)] <- 0
model.distribucion.h = mask(model.distribucion.h, ae.poly)
# Paramos
paramos.h[is.na(paramos.h)] <- 0
paramos.h= mask(paramos.h, ae.poly)
# Huella Humana
# Invertir datos para priorizar los que tengan menor huella humana
hh.h.i <- 100 - hh.h
hh.h.i= mask(hh.h.i, ae.poly)
# runap
runap.h[is.na(runap.h)] <- 0
runap.h= mask(runap.h, ae.poly)

all_layers <- c(costo.h.t,model.distribucion.h,paramos.h, hh.h.i, runap.h)
raster_stack <- stack(all_layers)
names(raster_stack) <- c("Costos", "Modelo distribucion", "Paramos", 'Huella Humana', 'RUNAP')
plot(rast(raster_stack), col = viridis(100))


# Calcular presupuesto
budget <- terra::global(costo.h.t, "sum", na.rm = TRUE)[[1]] * 0.2
print(paste(budget, 'Billones de pesos'))


#Reunir objetos de conservacion
feature_layers = c(paramos.h, model.distribucion.h, hh.h.i)
# Cmabiar nombre de capas
names(feature_layers) <- c("Paramos", "Model_distribucion", "Huella_Humana")
# plot(feature_layers,col = viridis(100), axes = FALSE)

#Crear problema
p1 <-
  problem(costo.h.t, features = feature_layers) %>%
  add_min_shortfall_objective(budget) %>%
  add_relative_targets(0.10) %>%
  add_binary_decisions() %>%
  add_default_solver(gap = 0.1, verbose = FALSE)

# Calcular el número de unidades de planeacion
number_of_planning_units(p1)

s1 <- solve(p1)
# Visualizar resultados
s_transformed <- projectRaster(raster(s1), crs = 3857, method = 'ngb')
raster_df <- as.data.frame(s_transformed, xy = TRUE)

ggplot() +  annotation_map_tile(zoom = 8, type = "osm") + 
  geom_raster(data = na.omit(raster_df), aes(x = x, y = y, 
                                             fill = factor(Beneficio_Neto_Total))) +
  scale_fill_manual(values = c("grey", 'blue'), na.value = NA) +  
  theme_minimal() +
  labs(title = "Areas priorizadas S1", fill = "")

p2 <-
  problem(costo.h.t, features = feature_layers) %>%
  add_min_shortfall_objective(budget) %>%
  add_relative_targets(0.1) %>%
  add_binary_decisions() %>%
  add_default_solver(gap = 0.1, verbose = FALSE) %>% 
  add_locked_in_constraints(runap.h) %>% 
  add_boundary_penalties(penalty = 0.001, edge_factor = 0.05)
s2 <- solve(p2)


s_transformed <- projectRaster(raster(s2), crs = 3857, method = 'ngb')
raster_df <- as.data.frame(s_transformed, xy = TRUE)

ggplot() +  annotation_map_tile(zoom = 8, type = "osm") + 
  geom_raster(data = na.omit(raster_df), aes(x = x, y = y, 
                                             fill = factor(Beneficio_Neto_Total))) +
  scale_fill_manual(values = c("grey", 'blue'), na.value = NA) +  
  theme_minimal() +
  labs(title = "Areas priorizadas S2", fill = "")


# Calcular puntajes de importancia
rc <- p2 %>% eval_ferrier_importance(s2)

# Preparar para plot
s_transformed <- projectRaster(raster(rc[["total"]]), crs = 3857, method = 'ngb')
raster_df <- as.data.frame(s_transformed, xy = TRUE)

# Graficar los resultados
ggplot() +  annotation_map_tile(zoom = 8, type = "osm") + 
  geom_raster(data = na.omit(raster_df), aes(x = x, y = y, 
                                             fill = total))+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
  theme_minimal() +
  labs(title = "Areas priorizadas", fill = "")


# costo de la solucion 2
costo.total = eval_cost_summary(p2, s2)
# Resultado en Billones de pesos
print(costo.total$cost)

# Resumen del cumplimiento de las metas por cada objeto de conservacion 
p2_target_coverage <- eval_target_coverage_summary(p2, s2)
print(p2_target_coverage)
# Extraer el porcentaje promedio del cimplimiento de las metas de conservacion
print(mean(p2_target_coverage$met) * 100)


# costo de la solucion 1
costo.total = eval_cost_summary(p1, s1)
# Resultado en Billones de pesos
print(costo.total$cost)

# Resumen del cumplimiento de las metas por cada objeto de conservacion 
p1_target_coverage <- eval_target_coverage_summary(p1, s1)
print(p1_target_coverage)
# Extraer el porcentaje promedio del cumplimiento de las metas de conservacion
print(mean(p1_target_coverage$met) * 100)





