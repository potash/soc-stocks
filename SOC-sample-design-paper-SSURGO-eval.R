
# need latest from GH
library(SoilTaxonomy)

library(aqp)
library(soilDB)
library(sf)

# need latest sharpshootR from GH
library(sharpshootR)

library(terra)
library(rasterVis)
library(viridisLite)


## copy / paste viewport bounding-box from SoilWeb
## click somewhere on the map
## press 'b', BBOX is copied to the clipboard


# https://casoilresource.lawr.ucdavis.edu/gmap/?loc=40.00715,-88.28981,z17
bb <- '-88.2983 40.0030,-88.2983 40.0109,-88.2819 40.0109,-88.2819 40.0030,-88.2983 40.0030'



## assemble AOI polygon into WKT
wkt <- sprintf('POLYGON((%s))', bb)

## init sf polygon
x <- st_as_sfc(wkt)

# set CRS as GCS WGS84
st_crs(x) <- 4326

## get overlapping map unit keys
# could also use SDA_query() with more elaborate SQL
m <- SDA_spatialQuery(x, what = 'mukey')

## compose SQL to return component details for these map unit keys
sql <- sprintf("SELECT mukey, cokey, compname, comppct_r FROM component WHERE mukey IN %s AND majcompflag = 'Yes'", format_SQL_in_statement(as.integer(m$mukey)))

## send to SDA, result is a data.frame
s <- SDA_query(sql)

## get OSD morphology + extended summaries 
osd <- fetchOSD(unique(s$compname), extended = TRUE)

## check out results
str(osd, 1)




# get gSSURGO grid here
# result is a terra SpatRaster
mu <- mukey.wcs(aoi = x, db = 'gssurgo')

# extract mukeys / RAT for thematic mapping
rat <- cats(mu)[[1]]

# TODO: annoying
# set mukey to integer
rat$mukey <- as.integer(rat$mukey)

# check:
plot(mu, col = viridis(100))

# get thematic data from SDA
# dominant component
# depth-weighted average
# sand, silt, clay (RV)
p <-  get_SDA_property(property = c("sandtotal_r","silttotal_r","claytotal_r"),
                       method = "Dominant Component (Numeric)", 
                       mukeys = rat$mukey,
                       top_depth = 0,
                       bottom_depth = 25)

head(p)
str(p)

cn <-  get_SDA_property(property = c("compname"),
                       method = "Dominant Component (Category)", 
                       mukeys = as.integer(rat$mukey))


# re-create raster attribute table with aggregate soil properties
rat <- merge(rat, p, by.x = 'mukey', by.y = 'mukey', sort = FALSE, all.x = TRUE)

# add compname
rat <- merge(rat, cn[, c('mukey', 'compname')], by.x = 'mukey', by.y = 'mukey', sort = FALSE, all.x = TRUE)

# TODO: annoying that base::merge() re-orders
#       mukey is a character inherited from rat

# re-pack RAT
levels(mu) <- rat

# convert raster + RAT --> stack of values
# note this includes extra, ID rasters
mu.data <- catalyze(mu)


## categorical data
activeCat(mu) <- 8

col.set <- RColorBrewer::brewer.pal(9, 'Set1')
cols <- colorRampPalette(col.set, space = 'Lab', interpolate = 'spline')(50)
plot(mu, col = cols)


# extract specific soil properties
ssc <- mu.data[[c("sandtotal_r","silttotal_r","claytotal_r")]]

# graphical check
# note implicit simplification via maxpixels
rasterVis::levelplot(
  ssc, 
  main = 'Sand, Silt, Clay (RV) 25-50cm\nDominant Component',
  margin = FALSE, 
  scales = list(draw = FALSE), 
  col.regions = viridis,
  maxpixels = 1e4
)

## soil texture class of fine earth fraction

# copy grid, will replace cell values
ssc.class <- ssc[[1]]

# classify sand, clay fractions
# retain all possible texture classes
values(ssc.class) <- ssc_to_texcl(
  sand = values(ssc[['sandtotal_r']]), 
  clay = values(ssc[['claytotal_r']]), 
  droplevels = FALSE
)

# name for RAT
names(ssc.class) <- 'soil.texture'

# check
plot(ssc.class, col = cols)

# note that all possible texture classes are included in the RAT
cats(ssc.class)[[1]]


## setup factor levels of ST
data("ST_unique_list")

osd$SPC$soilorder <- droplevels(factor(osd$SPC$soilorder, levels = ST_unique_list$order, ordered = TRUE))
osd$SPC$suborder <- droplevels(factor(osd$SPC$suborder, levels = ST_unique_list$suborder, ordered = TRUE))
osd$SPC$greatgroup <- droplevels(factor(osd$SPC$greatgroup, levels = ST_unique_list$greatgroup, ordered = TRUE))
osd$SPC$subgroup <- droplevels(factor(osd$SPC$subgroup, levels = ST_unique_list$subgroup, ordered = TRUE))

## convert horizon boundary distinctness -> vertical distance
# see manual page
osd$SPC$hzd <- hzDistinctnessCodeToOffset(
  osd$SPC$distinctness, 
  codes = c('very abrupt', 'abrubt', 'clear', 'gradual', 'diffuse')
)

## arrange sketches according to soil classification


# suggest ordering of dendrogram leaves
# catenary position
taxa.order <- c('typic argiudolls', 'oxyaquic argiudolls', 'aquic argiudolls', 'typic endoaquolls')
m <- match(osd$SPC$subgroup, taxa.order)

o <- order(rank(m, ties.method = 'first'))

profile_id(osd$SPC)[o]

plotSPC(osd$SPC, plot.order = o)

SoilTaxonomyDendrogram(
  spc = osd$SPC, 
  y.offset = 0.4, 
  scaling.factor = 0.016, 
  cex.taxon.labels = 1,
  cex.id = 0.66,
  cex.names = 0.66,
  width = 0.3, 
  name.style = 'center-center', 
  depth.axis = list(line = -3.5),
  hz.distinctness.offset = 'hzd'
)

SoilTaxonomyDendrogram(
  spc = osd$SPC, 
  rotationOrder = profile_id(osd$SPC)[o],
  y.offset = 0.4, 
  scaling.factor = 0.016, 
  cex.taxon.labels = 1,
  cex.id = 0.85,
  cex.names = 0.75,
  width = 0.3, 
  name.style = 'center-center', 
  depth.axis = list(line = -3.5),
  hz.distinctness.offset = 'hzd'
)

SoilTaxonomyDendrogram(
  spc = osd$SPC, 
  rotationOrder = profile_id(osd$SPC)[order(osd$SPC$subgroup)],
  y.offset = 0.4, 
  scaling.factor = 0.016, 
  cex.taxon.labels = 1,
  cex.id = 0.85,
  cex.names = 0.75,
  width = 0.3, 
  name.style = 'center-center', 
  depth.axis = list(line = -3.5),
  hz.distinctness.offset = 'hzd'
)



## 3D geomorphic summary

# there may be records missing from SPC / geomorphic component
nm <- intersect(profile_id(osd$SPC), osd$geomcomp$series)

# keep only those series that exist in both
sub <- subset(osd$SPC, profile_id(osd$SPC) %in% nm)

## inverse problem: extra records in geomcomp summaries
# subset geomcopm
geomcomp.sub <- subset(osd$geomcomp, subset = series %in% profile_id(sub))

# viz geomorphic proportion summary, results contain clustering object
res <- vizGeomorphicComponent(geomcomp.sub)
print(res$fig)


# arrange according to clustering of geomorphic component
par(mar = c(0, 0, 0, 0))
plotProfileDendrogram(
  sub,
  clust = res$clust,
  dend.y.scale = 3, y.offset = 0.2,
  scaling.factor = 0.01,
  width = 0.3,
  name.style = 'center-center',
  depth.axis = FALSE,
  hz.depths = TRUE,
  hz.distinctness.offset = 'hzd',
  cex.names = 0.6,
  cex.id = 0.6
)


## 2D geomorphic summary
# there may be records missing from SPC / hill slope position
nm <- intersect(profile_id(osd$SPC), osd$hillpos$series)

# keep only those series that exist in both
sub <- subset(osd$SPC, profile_id(osd$SPC) %in% nm)

## inverse problem: extra records in hill slope summaries
# subset hillpos
hillpos.sub <- subset(osd$hillpos, subset = series %in% profile_id(sub))

# viz hill slope proportion summary, results contain clustering object
res <- vizHillslopePosition(hillpos.sub)
print(res$fig)


# arrange according to clustering of hillslope position
par(mar = c(0, 0, 0, 0))
plotProfileDendrogram(
  sub, 
  clust = res$clust, 
  dend.y.scale = 3, y.offset = 0.2,
  scaling.factor = 0.01, 
  width = 0.3, 
  name.style = 'center-center', 
  depth.axis = FALSE, 
  hz.depths = TRUE, 
  hz.distinctness.offset = 'hzd', 
  cex.names = 0.6, 
  cex.id = 0.6
)


## borrowing ideas from this tutorial:
## https://ncss-tech.github.io/AQP/soilDB/exploring-geomorph-summary.html
##
hp.cols <- RColorBrewer::brewer.pal(n = 5, name = 'Set1')[c(2, 3, 4, 5, 1)]

# re-order hillslope proportions according to clustering
hp <- hillpos.sub[res$order, ]
nm <- names(hp[, 2:6])

par(mar = c(0.5, 0, 0, 2))
layout(matrix(c(1,2)), widths = c(1,1), heights = c(2,1))
plotProfileDendrogram(sub, res$clust, dend.y.scale = 3, scaling.factor = 0.012, y.offset = 0.2, width = 0.32, name.style = 'center-center', cex.names = 0.7, shrink = TRUE, cex.id = 0.55)

## TODO: encode Shannon entropy: values are computed row-wise, data plotted as columns
matplot(y = hp[, 2:6], type = 'b', lty = 1, pch = 16, axes = FALSE, col = hp.cols, xlab = '', ylab = '', xlim = c(0.5, length(sub) + 1))
# grid(nx = 0, ny = NULL)
axis(side = 4, line = -1, las = 1, cex.axis = 0.7)
# axis(side = 2, line = -3, las = 1, cex.axis = 0.7)
legend('top', legend = rev(nm), col = rev(hp.cols), pch = 16, bty = 'n', cex = 0.8, pt.cex = 2, horiz = TRUE, inset = c(0.01, 0.01))
mtext('Probability', side = 2, line = -2, font = 2)


## TODO: encode Shannon entropy
par(mar = c(0.5, 0, 0, 2))
layout(matrix(c(1,2)), widths = c(1,1), heights = c(2,1))
plotProfileDendrogram(sub, res$clust, dend.y.scale = 3, scaling.factor = 0.012, y.offset = 0.2, width = 0.32, name.style = 'center-center', cex.names = 0.7, shrink = TRUE, cex.id = 0.55)

sp <- c(1.5, rep(1, times = length(sub) - 1))
barplot(height = t(as.matrix(hp[, 2:6])), beside = FALSE, width = 0.5, space = sp, col = hp.cols,  axes = FALSE, xlab = '', ylab = '', xlim = c(0.5, length(sub) + 1), ylim = c(0, 1.2))

legend(x = 0.75, y = 1.2, legend = rev(nm), col = rev(hp.cols), pch = 15, bty = 'n', cex = 0.8, pt.cex = 1.25, horiz = TRUE)
mtext('Probability', side = 2, line = -2, font = 2)



