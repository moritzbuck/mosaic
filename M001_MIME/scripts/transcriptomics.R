library(data.table)
library(ggplot2)
library(pheatmap)
library(vegan)
library(ggrepel)

#mapping_folder = "/home/moritz/data/M001_MIME/mappings/MOSAIC-MIME-RNA"
mapping_folder = "/home/moritz/temp/mappings/MOSAIC-MIME-RNA"
binset_folder = "/home/moritz/data/M001_MIME/binsets/COMPLETE-SET-AND-BULK"
md_folder = "/home/moritz/projects/mosaic/M001_MIME/config_nd_metadata/"

#rpkm = fread(file.path(mapping_folder, "tpm.csv"))
#genes = rpkm$derep_gene
#rpkm[,derep_gene := NULL]

#rpkm = fread(file.path(mapping_folder, "KO_relative_abundance.csv"))
#rpkm = rpkm[!is.na(KO)]
#genes = rpkm$KO
#rpkm[,KO := NULL]
#genes = genes[rowSums(rpkm) > 0]
#rpkm = rpkm[rowSums(rpkm) > 0]

rpkm = fread(file.path(mapping_folder, "COG_category_relative_abundance.csv"))
rpkm = rpkm[!is.na(COG_category)]
genes = rpkm$COG_category
rpkm[,COG_category := NULL]
genes = genes[rowSums(rpkm) > 0]
rpkm = rpkm[rowSums(rpkm) > 0]


gene_md = fread(file.path(binset_folder, "arctic_gene_atlas_basics.csv"))
setkey(gene_md, "cluster_id")
md = fread(file.path(md_folder, "exp_data.csv"))
setkey(md, "V1")
md[, V2 := NULL]

sample_pca = function(samps, x="PC1", y="PC2", labs = "treatment", col = "light", title ="" ){
exp22 = rpkm[,..samps]
genes22 = genes[rowSums(exp22) > 0]
exp22 = exp22[rowSums(exp22) > 0,]
not_hypos = gene_md[genes22]$hypo_func != "hypothetical protein"

pca = predict(prcomp(vegdist(t(exp22))))[,1:4]
pca = data.table(pca)
pca[, sample := names(exp22)]
pca = cbind(pca,md[pca$sample])

ggplot(pca, aes_string(x=x, y=y, col=col, label=labs)) +geom_point()+geom_text_repel(max.overlaps=30)+theme_minimal()+ggtitle(title )
}

sample_pca(grep("Sample22",names(rpkm), value=TRUE), title = "Sample 22")
sample_pca(grep("Sample30",names(rpkm), value=TRUE), title = "Sample 30")
sample_pca(grep("Sample30",names(rpkm), value=TRUE), title = "Sample 30")
sample_pca(grep("[DL]", names(rpkm), value=TRUE, invert = TRUE), title = "time 0s", labs = "experiment")
sample_pca(names(rpkm), labs = "experiment" , title = "all RNA")

temp = copy(rpkm)
temp[, genes := genes]

melted = melt(temp, id.vars = "genes", variable.name = "sample", value.name="rpkm")
rm(temp)

gene_md[, gtdbtk_classif := NULL]
gene_md[, representative_nucls := NULL]
gene_md[, representative_aas := NULL]
gene_md[, mOTU := NULL]
gene_md[, members := NULL]

melted = cbind(melted, gene_md[melted$genes])

by_symbol = melted[,.(rpkm = sum(rpkm)),by=c('sample','symbol')]
by_symbol = cbind(by_symbol,md[as.vector(by_symbol$sample)])

casty = dcast(by_symbol, formula=symbol~V1, value.var='rpkm')
casty = casty[!symbol %in% c(NA,"-") ]
symbs = casty$symbol

by_ko = melted[,.(rpkm = sum(rpkm)),by=c('sample','KO')]
by_ko = cbind(by_ko,md[as.vector(by_ko$sample)])


setkey(by_symbol, "symbol")
setkey(by_ko, "ko")


by_symbol[!is.na(symbol) & symbol !=  "-"]
