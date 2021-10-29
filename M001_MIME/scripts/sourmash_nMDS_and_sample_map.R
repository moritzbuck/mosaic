library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("gghighlight")
library(data.table)
library(vegan)

world <- ne_countries(scale = "medium", returnclass = "sf")
setwd("/home/moritz/kadath/data/M001_MIME/dbs")
context_md = fread("context_md.csv")
setkey(context_md, "V1")
context_md['ERR3589578']$lat = 74.8283
context_md['ERR3589578']$lon = 76.2096
oden_md = data.frame(fread("~/kadath/projects/mosaic/M001_MIME/config_nd_metadata/oden_general_table_samples.csv"))
mosaic_md = data.frame(fread("~/kadath/projects/mosaic/M001_MIME/config_nd_metadata/mosaic_general_table_sample.csv"))
mosaic_md$sample.nr.in.metadata = gsub("?","", mosaic_md$sample.nr.in.metadata, fixed = TRUE)
mosaic_md$sample.nr.in.metadata = gsub(" ","", mosaic_md$sample.nr.in.metadata, fixed = TRUE)
mosaic_md[mosaic_md$sample_label == "PS122_CN_DNA_703", ]$sample.nr.in.metadata = "703"
mosaic_md = mosaic_md[mosaic_md$sample.nr.in.metadata != "",]
mosaic_md = mosaic_md[mosaic_md$sample.nr.in.metadata != "nosample",]
row.names(mosaic_md) = gsub("^", "MOSAIC-MIME-DNA-", mosaic_md$sample.nr.in.metadata)
mosaic_md = mosaic_md[grep("MIME",row.names(sig_simis), val=TRUE) ,]
for(col in names(context_md)) set( context_md, i=which(context_md[[col]]=='not collected'), j=col, value=NA)
oden_md =  oden_md[grepl("Seawater", oden_md$sample_name), ]
row.names(oden_md) = gsub("Seawater-", "ODIN-2013-DNA-", oden_md$sample_name)

numeric_cols = c("chlorophyll", "salinity", "tot_nitro")
for(col in numeric_cols) context_md[, col := as.numeric(gsub(",",".",context_md[,get(col)]))]

context_md$temp = as.numeric(gsub(",",".",context_md$temp))

month2temp = function(x) abs((x +6 -2) %% 12 -6)

mosaic_gps = read.csv("~/kadath/projects/mosaic/M001_MIME/config_nd_metadata/mosaic_gps.tsv", row.names=1, sep="\t")

sig_simis = read.csv("samples_nd_context_all_simis_w_abundance.csv", row.names=1, as.is = TRUE)

selection = !grepl("ODafIN", row.names(sig_simis))

#mds = data.table(monoMDS(t(as.dist(1-sig_simis[selection,selection])), model = "global", stress = 2)$points)
mds = data.table(predict(prcomp(1-sig_simis[selection,selection])))[,1:5]

mds[, samples := row.names(sig_simis)[selection]]

mds[, set := 'context']
mds[grepl("MOSAIC",mds$samples), set := 'mosaic']
mds[grepl("ODIN",mds$samples), set := 'odin-2013']
mds[, context := context_md[samples]$context]
mds[, lon := context_md[samples]$lon]
mds[, lat := context_md[samples]$lat]
mds[grepl("MOSAIC",mds$samples), lat := mosaic_md[grep("MOSAIC",mds$samples, val=TRUE), "approx_lat"]]
mds[grepl("ODIN",mds$samples), lat := oden_md[grep("ODIN",mds$samples, val=TRUE), "latitude"]]
mds[, ocean := "atlantic"]
mds[lon < -89.2 | lon > 125, ocean := "pacific"]
mds[, month := sapply(strsplit(context_md[samples]$collection_date, "-"), function(x) month2temp(as.numeric(x[2])) )]
mds[, tot_nitro := context_md[samples]$tot_nitro]
mds[, salinity := context_md[samples]$salinity]


#mds[, temp := abs(context_md[samples]$temp)]

#ggplot(mds, aes(x=PC1, y=PC2, col=set))+geom_point(size=20)+theme_bw()+theme(text = element_text(size=50))
ggplot(mds, aes(x=PC1, y=PC2, shape=set, col=abs(lat)))+geom_point(size=20)+theme_bw()+theme(text = element_text(size=50)) + scale_color_viridis_c(option = "magma")

sub = context_md[lon != "NA", ]
context_gps =  context_md[lon != "NA", c("V1","lat","lon", "context")]
oden_gps = oden_md[oden_md$Contents == "Seawater",c("Station", "latitude","longitude")]
context_sf = st_as_sf(x = context_gps, coords = c("lon", "lat"), crs=  "+proj=longlat")
oden_sf = st_as_sf(x = oden_gps, coords = c("longitude", "latitude"), crs=  "+proj=longlat")
mosaic_sf = st_as_sf(x = mosaic_md, coords = c("approx_lon", "approx_lat"), crs=  "+proj=longlat")

ggplot(data = world) + geom_sf(fill="floralwhite", alpha=0.5, col="lightgray")+geom_sf(data=context_sf, shape=18, size=3, mapping=aes(col = context))+
geom_sf(data=oden_sf, shape=17, size=6, col = "red")+geom_sf(data=mosaic_sf, shape=17, size=6, col = "green")

ggsave("~/kadath/projects/mosaic/M001_MIME/figures/World.pdf", width=17, height = 9)


proj="+proj=laea +lat_0=60 +lon_0=0 +ellps=WGS84 +units=m"
zoom=1
world_pretty = st_transform(world, crs = proj)
context_pretty = st_transform(st_as_sf(x = context_gps[context_gps$lat > 0,], coords = c("lon", "lat"), crs=  "+proj=longlat"), crs = proj)
oden_pretty = st_transform(oden_sf, crs = proj)
mosaic_pretty = st_transform(mosaic_sf, crs = proj)

xwind = c(min(st_coordinates(context_pretty)[,"X"]), max(st_coordinates(context_pretty)[,"X"]))*zoom
ywind = c(min(st_coordinates(context_pretty)[,"Y"]), max(st_coordinates(context_pretty)[,"Y"]))*zoom

ggplot(data = world_pretty) + geom_sf(fill="floralwhite", alpha=0.5, col="lightgray")+
  geom_sf(data=context_pretty, mapping = aes(col=context), shape=18, size=3)+geom_sf(data=mosaic_pretty, shape=17, size=6, col = "orange")+geom_sf(data=oden_pretty, shape=17, size=6, col = "purple")+
  coord_sf(xlim=xwind, ylim=ywind, expand=TRUE)

ggsave("~/kadath/projects/mosaic/M001_MIME/figures/Northern.pdf", width=10, height = 9)


proj="+proj=laea +lat_0=-60 +lon_0=-30 +ellps=WGS84 +units=m"
zoom=1.5
world_pretty = st_transform(world, crs = proj)
context_pretty = st_transform(st_as_sf(x = context_gps[context_gps$lat < 0,], coords = c("lon", "lat"), crs=  "+proj=longlat"), crs = proj)

xwind = c(min(st_coordinates(context_pretty)[,"X"]), max(st_coordinates(context_pretty)[,"X"]))*zoom
ywind = c(min(st_coordinates(context_pretty)[,"Y"]), max(st_coordinates(context_pretty)[,"Y"]))*zoom

ggplot(data = world_pretty) + geom_sf(fill="floralwhite", alpha=0.5, col="lightgray")+
  geom_sf(data=context_pretty, mapping = aes(col=context), shape=18, size=3)+coord_sf(xlim=xwind, ylim=ywind, expand=TRUE)

ggsave("~/kadath/projects/mosaic/M001_MIME/figures/Southern.pdf", width=7, height = 9)
