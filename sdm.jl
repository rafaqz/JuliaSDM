#----
# libraries
nothing
using Rasters
using ArchGDAL
import GeometryOps as GO
using GeoDataFrames
using GLMakie
using ImageFiltering
using DataFrames
using Statistics

# if we want to use the functions from the oneimpact package to compute the ZOI/density
# remotes::install_github("NINAnor/oneimpact", ref = "HEAD")

#----
# step 1: rasterization

# list vector layers
filepaths = readdir("data/"; join=true)
gpkgs = filter(x -> occursin(".gpkg", x), filepaths)
gpkg_labels = Symbol.(first.(splitext.(basename.(gpkgs))))
geoms = NamedTuple(gpkg_labels .=> gpkgs)

# load first layer
geoms.reindeer_area
template = GeoDataFrames.read(geoms.reindeer_area)
# set a grid for rasterization based on the extent
# we could make the resolution as small as wanted if we want to test rasters with
# finer resolution
# load all vectors and rasterize
res = 50
st = map(geoms) do path
    layer_vect = GeoDataFrames.read(path)
    rasterize(count, layer_vect; 
        to=template, res, missingval=0
    )
end |> RasterStack

# Plot the input data
Rasters.rplot(st) # rplot is for multi-plot stacks, Makie.plot will do this one day

#----
# step 2: compute the zone of influence - in R

# we do this in R
# but we could also compute that in GRASS to test it out and benchmark

# here we compute that for a single radius/moving window size, 5000m, for illustration

radus = 1000
r = radus ÷ res
d = 2r + 1
@time layers_zoi = maplayers(st[(:reindeer_cabins, :reindeer_roads_public, :reindeer_roads_private)]) do rast
    data = mapwindow(sum, rast, (d, d))
    Raster(data, dims(rast); missingval=0)
end

# Plot zone of influence results
Rasters.rplot(layers_zoi)

# but we could also compute that for many radii, to increase the number of layers to
# be annotated to the biological data

# radii = [100, 500, 2000]
# layers_zois = map(radii) do radius
#     r = radius ÷ res
#     d = 2r + 1
#     @show d
#     @time maplayers(st) do rast
#         mapwindow(sum, rast, (d, d))
#     end
# end
# check
# Rasters.rplot(layers_zois[end])


#----
# step 3: extract raster values for biological data - points and steps
# (but maybe we could also try polygons?)

#----
# step 3.1: extract covariates for points

# read vector
reindeer_pts = GeoDataFrames.read("data/reindeer_points.gpkg")
Makie.plot(reindeer_pts.geom)

# extract points
@time reindeer_covar_pts = extract(layers_zoi, reindeer_pts)

#----
# step 3.2: extract covariates for lines

# read vector
reindeer_lines = GeoDataFrames.read("data/reindeer_lines.gpkg")
Makie.plot(reindeer_lines.geom)

# Take zonal statistics for lines
@time means = zonal(mean, layers_zoi; of=reindeer_lines, threaded=true);

#----
# step 4 - spatial prediction - we can develop that later, we can add here a simple GLM-type model and then predict it in space.

using Pkg

import SpeciesDistributionModels as SDM

