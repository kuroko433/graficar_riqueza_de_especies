##### script para graficar riqueza de especies en 2D #####

##### de paso para aprender a trabajar datos de PBDB

rm(list = ls())

#### establecer directorio de trabajo ####
setwd("direccion_de_tu_carpeta_de_trabajo")

###### librerias ########
library(dplyr) ##para manejar datos
library(ggplot2) ##para graficar
library(sf) ## para manejar datos espaciales
library(ggspatial) ##para graficar
library(hrbrthemes) ##para graficar
library(ggthemes) ##para graficar
library(gtools)
library(terra) ##manejar rasters
library(velociraptr) ##descargar datos fosiles desde paleobiology database
library(rgplates) ### reconstruir coordenadas paleogeograficas
library(purrr) ##### calcular riqueza de especies
library(smoothr) #### calcular riqueza de especies
library(vegan) ##manejar datos ecologicos
library(glue)
library(viridis)
##### riqueza mioceno ###########
## vamos a descargar datos de caballos fosiles (familia Equidae) y sus relativos del Mioceno
## a nivel de Genero (23.03 a 5.33 Ma atras)
caballo_miocene=downloadPBDB(Taxa="Equidae",StartInterval="Miocene",StopInterval="Miocene")

## quitamos ocurrencias que no llegan a genero
caballo_miocene=subset(caballo_miocene,!caballo_miocene$genus=="")
caballo_miocene<-caballo_miocene[,-c(16:18)]

### veamos cuantas especies tenemos
collections_miocene=tapply(caballo_miocene$collection_no,list(caballo_miocene$collection_no,caballo_miocene$genus),length)

collections_miocene[is.na(collections_miocene)]=0
collections_miocene[collections_miocene>1]=1
ngenus_miocene<-ncol(collections_miocene)

##tenemos 32 generos de caballos para el Mioceno
print(ngenus_miocene)

## veamos la completitud de nuestro ensamble (aproximacion tipo presencia-ausencia por coleccion)


specpool(collections_miocene) ##basado en CHAO 2 abrian 47 generos totales para el Mioceno

print(glue("en total tenemos {specpool(collections_miocene)[1]} generos para el Mioceno,
           y una completitud de ensambles del {round((specpool(collections_miocene)[1]/specpool(collections_miocene)[2])*100,2)}%"))


### vamos acalcular las paleolatitudes de nuestras pcurrencias ####
web<-getgws() ##llink a gplates

setgws(web) ##conectar a gplates

miocene_coordinates<-data.frame(lon=caballo_miocene$lng,lat=caballo_miocene$lat)

### ponemos la reconstruccion en 20 millones de años
### si decides cambiar el codigo a distintas edades del Mioceno
### puede que te de resultados ligeramente diferentes a los mios
### no es lo mas adecuado poner esta edad ya que el abanico de edades
### de las ocurrencias es muy grande, pero por fines de experimentacion
### y demostrativos plotearemos todo a los 20 millones de anos

miocene_coordinates<-as.data.frame(reconstruct(miocene_coordinates,age=20))

## juntamos datos
caballo_miocene=cbind(caballo_miocene[,c(8)],miocene_coordinates)

caballo_miocene=na.omit(caballo_miocene) ##quitamos ocurrencias que no se pudieron reconstruir

colnames(caballo_miocene)[1]<-c("Genus") ##arreglar nombre de columnas

####cargamos paleomapas
miocene_poly = st_read("paleo_map/20Ma_CM_v7.shp")


###mapa de los puntos####

miocene_points <- ggplot() + 
  geom_sf(data = miocene_poly, fill = "grey70", col = 'black', lwd = 0.2)+ 
  geom_point(data = caballo_miocene, aes(x = paleolong, y = paleolat),col="#8B2323")+
  labs(x = '', y = '') + ##Nombres de los ejes x e y
  ggtitle(label = 'Sampling Points') + ##titulo del mapa
  theme_ipsum_es() + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5), 
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom', 
        legend.key.width = unit(3.5, 'line'),
        panel.background = element_rect(fill = "white"),
        plot.subtitle = element_text(face = 'italic', hjust = 0.5)) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(8, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering())

###para guardar tu grafico
ggsave(plot = miocene_points, filename = glue('miocene_samples.jpg'), units = 'in', width = 10, height = 7, dpi = 300)

#####calculo de diversidad gamma ####

####pasamos nuestros puntos a geometrias
mio.caballo_sf <- st_as_sf(caballo_miocene, coords = c("paleolong", "paleolat"), crs = 4326)

### creamos un buffer
mio.buff <- sf::st_buffer(mio.caballo_sf, dist = 500000) |> 
  sf::st_union() |> 
  sf::st_sf() |>
  sf::st_transform(crs = terra::crs(mio.caballo_sf))

### creamos poligonos para cada genero
mio.caballo_poly <- mio.caballo_sf %>%
  group_by(Genus) %>%
  summarize(geometry = st_union(geometry)) %>%
  st_convex_hull()

### creamos una grilla espacial donde vamos a calcular la riqueza de generos por celda
### de acuerdo a la spuerposicion espacial de los generos de caballos
### en este caso es una grilla de 2 grado x 2 grado de latitud
### tambien puedes modificar esto e ir experimentando

mio.grid <- st_make_grid(mio.caballo_sf, cellsize = c(2,2))

###hacemos la interseccion entre los poligonos de nuestros generos, las celdas espaciales
mio.caballo_sr <- st_intersects(mio.grid, mio.caballo_poly)

### ahora contamos cuantos generos de caballo hay por celda espacial
mio.caballo_sr_count <- sapply(mio.caballo_sr, function(x) length(unique(mio.caballo_poly$Genus[x])))

### pasamos a dataframe
mio.richness_df <- data.frame(
  geometry = st_geometry(mio.grid),
  richness = mio.caballo_sr_count
)

### volvemos ahora pasamos a poligono (o spatial feature)
mio.richness_sf <- st_as_sf(mio.richness_df)

#### una vista rapida de como se ve nuestra grilla 

ggplot() +
  geom_sf(data = mio.richness_sf, aes(fill = richness)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Mapa de Diversidad Gamma")

### ahora recortaremos nuestra grilla de acuerdo a nuestro buffer espacial

ext(mio.richness_sf) ##extencion de nuestro spatial feature

# Crear un ráster vacío con la misma extensión y resolución que nuestro spatial feature
mio.raster_template <- rast(ext(mio.richness_sf), resolution = c(2,2), crs = crs(mio.richness_sf))

# Rasterizar el objeto sf
mio.richness_raster <- rasterize(mio.richness_sf, mio.raster_template, field = "richness")

## lo que hicimos fue rasterizar nuestra riqueza de especies

plot(mio.richness_raster)

##ahora vamos a recortar nuestro raster de acuerdo a los paleocontinentes
mio.richness_raster <-terra::crop(mio.richness_raster,miocene_poly)
mio.richness_raster <-terra::mask(mio.richness_raster,miocene_poly)

##vistazo rapido
plot(mio.richness_raster)

##ahora vamos a recortar nuestro raster de acuerdo a nuestro buffer espacial
mio.richness_raster <-terra::crop(mio.richness_raster,mio.buff)
mio.richness_raster <-terra::mask(mio.richness_raster,mio.buff)
mio.richness_raster[mio.richness_raster == 0] <- NA

##vistazo rapido
plot(mio.richness_raster)

##pasamos a dataframe

mio.gamma_sr_df<-terra::as.data.frame(mio.richness_raster, xy = T)

##### vamos a plotear la riqueza de especies #####

miocene_diversity <- ggplot() + 
  geom_sf(data = miocene_poly, fill = "grey70", col = 'black', lwd = 0.2)+ 
  geom_tile(data = mio.gamma_sr_df, aes(x = x, y = y, fill = richness))+
  scale_fill_viridis(option="turbo") +
  labs(x = '', y = '', fill="Genus Richness") + ##Nombres de los ejes x e y
  ggtitle(label = 'Miocene Equidae Gamma Richness') + ##titulo del mapa
  theme_ipsum_es() + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5), 
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom', 
        legend.key.width = unit(3.5, 'line'),
        panel.background = element_rect(fill = "white"),
        plot.subtitle = element_text(face = 'italic', hjust = 0.5)) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(8, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering())

###para guardar tu grafico
ggsave(plot = miocene_diversity, filename = glue('miocene_richness.jpg'), units = 'in', width = 10, height = 7, dpi = 300)


