library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("gghighlight")
library(data.table)
world <- ne_countries(scale = "medium", returnclass = "sf")

context_md = fread("context_md_with_manual.csv")
setkey(context_md, "short_name")

sig_simis = read.csv("samples_nd_context_all_simis_w_abundance.csv", row.names=1, as.is = TRUE)

selection = !grepl("ODIadsN", row.names(sig_simis))

mds = data.table(monoMDS(t(as.dist(1-sig_simis[selection,selection])), model = "global")$points, stress = 2)
mds = data.table(predict(prcomp(1-sig_simis[selection,selection])))[,1:5]

mds[, samples := row.names(sig_simis)[selection]]

mds[, set := 'context']
mds[grepl("MOSAIC",mds$samples), set := 'mosaic']
mds[grepl("ODIN",mds$samples), set := 'odin-2013']
mds[, context := context_md[samples]$context]

ggplot(mds, aes(x=PC3, y=PC2, col=set))+geom_point(size=20)+theme_bw()+theme(text = element_text(size=50))
