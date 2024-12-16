#----
# working directory
setwd("/data/P-Prosjekter/41203800_oneimpact/04_tools/julia_to_rafael/")

#----
# libraries
# install.packages("package_name") # if the packages are not installed
library(terra) # for raster and vectors
library(sf) # for vectors, maybe we don't need it
library(tictoc) # for computing time

# if we want to use the functions from the oneimpact package to compute the ZOI/density
# remotes::install_github("NINAnor/oneimpact", ref = "HEAD")
library(oneimpact) # for computing ZOI
# more details on the package here: https://ninanor.github.io/oneimpact/
# and here: https://github.com/NINAnor/oneimpact

#----
# step 1: rasterization

# list vector layers
vects <- list.files("data/", pattern = ".gpkg", full.names = TRUE) |>
  grep(pattern = "cabins|roads", value = TRUE)

# load first layer
infrastructure <- terra::vect(vects[1])
# set a grid for rasterization based on the extent
# we could make the resolution as small as wanted if we want to test rasters with
# finer resolution
rr <- terra::rast(xmin = terra::ext(infrastructure)[1], resolution = 500,
                  extent = terra::ext(infrastructure), crs = terra::crs(infrastructure))
# load all vectors and rasterize
layers <- list()
for(i in 1:length(vects)) {
  layer_vect <- terra::vect(vects[i])
  if(grepl("cabins", vects[i]))
    layers[[i]] <- terra::rasterize(layer_vect, rr, fun = length)
  else
    layers[[i]] <- terra::rasterize(layer_vect, rr)
  layers[[i]] <- terra::ifel(is.na(layers[[i]]), 0, layers[[i]]) # set 0 where it is NA
}
layers_rast <- terra::rast(layers)
names(layers_rast) <- c("cabins", "roads_private", "roads_public")
plot(layers_rast)

#----
# step 2: compute the zone of influence - in R

# we do this in R
# but we could also compute that in GRASS to test it out and benchmark

# here we compute that for a single radius/moving window size, 5000m, for illustration
layers_zoi <- oneimpact::calc_zoi_cumulative(layers_rast, type = "bartlett" , radius = 5000)
plot(layers_zoi)

# but we could also compute that for many radii, to increase the number of layers to
# be annotated to the biological data
radii <- seq(1000, 10000, by=1000)
layers_zoi <- lapply(radii, function (x) oneimpact::calc_zoi_cumulative(layers_rast, type = "bartlett" , radius = x))
layers_zoi <- terra::rast(layers_zoi)
# check
plot(layers_zoi)

# the function above basically calls some helper functions to define the
# weight matrices for moving windows (we can check ?oneimpact::create_filter)
# and then uses terra::focal() to do the moving window process.
?terra::focal

#----
# step 3: extract raster values for biological data - points and steps
# (but maybe we could also try polygons?)

# create datasets - we can skip this
# data("reindeer")
# reindeer
# # points
# reindeer_pts <- terra::vect(as.data.frame(reindeer), geom = c("x", "y"), crs = "EPSG:25833", keepgeom = TRUE)
# plot(reindeer_pts)
# # export
# terra::writeVector(reindeer_pts, filename = "/data/scratch/tmp_bernardo/reindeer_points.gpkg", overwrite = TRUE)
# file.copy(from = "/data/scratch/tmp_bernardo/reindeer_points.gpkg",
#           to = "data/", overwrite = TRUE)
# terra::writeVector(reindeer_pts, "data/reindeer_points.shp")
#
# # steps
# linestrings <- list()
# reindeer_pts_sf <- sf::st_as_sf(reindeer_pts)
# inds <- unique(reindeer$animal_year_id)
# for(i in seq_along(inds)) {
#   print(i)
#   ind <- inds[i]
#   rr <- reindeer_pts_sf |> dplyr::filter(animal_year_id == ind)
#   linestring <- lapply(X = 1:(nrow(rr)-1), FUN = function(x) {
#     pair <- sf::st_combine(c(rr[x,]$geometry, rr[x + 1,]$geometry))
#     line <- dplyr::bind_cols(sf::st_drop_geometry(rr[x,]),
#                              geom = sf::st_cast(pair, "LINESTRING")) |>
#       sf::st_as_sf()
#
#     if(difftime(rr$t[x+1], rr$t[x], units = "hours") > 5) line <- NULL
#     return(line)
#   })
#   linestrings[[i]] <- dplyr::bind_rows(linestring)
# }
# reindeer_lines <- dplyr::bind_rows(linestrings)
# plot(terra::vect(reindeer_lines))
# # export
# sf::st_write(reindeer_lines, dsn = "/data/scratch/tmp_bernardo/reindeer_lines.gpkg")
# file.copy(from = "/data/scratch/tmp_bernardo/reindeer_lines.gpkg",
#           to = "data/", overwrite = TRUE)
# sf::st_write(reindeer_lines, dsn = "data/reindeer_lines.shp", delete_dsn = TRUE)

#----
# step 3.1: extract covariates for points

# read vector
reindeer_pts <- terra::vect("data/reindeer_points.gpkg")
plot(reindeer_pts)

# extract for points
tic()
reindeer_covar_pts <- terra::extract(layers_zoi, reindeer_pts)
toc()

#----
# step 3.2: extract covariates for lines

# read vector
reindeer_lines <- terra::vect("data/reindeer_lines.gpkg")
plot(reindeer_lines)

# extract for points
tic()
reindeer_covar_lin <- terra::extract(layers_zoi, reindeer_lines)
# then we need to summarize values by ID - for instance, compute the mean
reindeer_covar_lin_split <- split(reindeer_covar_lin, reindeer_covar_lin$ID)
str(reindeer_covar_lin_split, max.level = 1)
reindeer_covar_lin_summary <- lapply(reindeer_covar_lin_split,
                                     function (x) apply(x, MARGIN = 2, mean, na.rm = TRUE)) |>
  do.call(what = "rbind") |>
  data.frame()
reindeer_covar_lin_summary
toc()


#----
# step 4 - spatial prediction - we can develop that later, we can add here a simple GLM-type model and then predict it in space.
