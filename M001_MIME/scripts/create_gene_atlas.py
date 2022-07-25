import os
from tqdm import tqdm
from os.path import join as pjoin
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas
import json
from subprocess import call
import pandas

seqid = 0.95
cov = 0.8
covmode = 1
temp_folder = "/home/moritz/temp"
threads = 24
path = "/home/moritz/temp/binsets/COMPLETE-SET-AND-BULK/"
fpath = "/home/moritz/data/M001_MIME/binsets/COMPLETE-SET-AND-BULK/"

if not os.path.exists(pjoin(path,  "arctic_gene_atlas.ffn")):
    for file in tqdm(os.listdir(pjoin(path, "clean_bins"))):
        call("cat {file} >> {ass}".format(file = pjoin(path, "clean_bins", file, file +".ffn"), ass = pjoin(path,  "arctic_gene_atlas.ffn")), shell=True)

bindat = pandas.read_csv(pjoin(path,  "COMPLETE-SET-AND-BULK_basics.csv" ), index_col=0)

call(f"mmseqs easy-cluster --min-seq-id {seqid} --cov-mode {covmode} -c {cov} --threads {threads} {path}/arctic_gene_atlas.ffn {path}/arctic_gene_atlas {temp_folder}/mmseqs_temp > {logfile}", shell=True)

clusters = []
record = None
count = 1
put_new_record = False
member_buffer = []
with open(f"{path}/arctic_gene_atlas_all_seqs.fasta") as handle:
    for l in tqdm(handle):
        if l.startswith(">") and not put_new_record :
            cluster_name = l[1:-1]
            put_new_record = True
            member_buffer += [l[1:-1].split(" ")[0]]
        elif l.startswith(">"):
            if put_new_record:
                if record  :
                    record['members'] = member_buffer[:-1]
                    clusters += [record]
                    member_buffer = [l[1:-1].split(" ")[0]]
                record = { "cluster_id" : f"arctic_gene_atlas-{count:06d}", "hypo_func" : " ".join(l[1:].split(" ")[1:]).strip(), "representative" : cluster_name}
                count += 1
            else :
                put_new_record = False
        else :
            if put_new_record:
                record["representative_nucls"] = l[:-1]
                record["representative_aas"] = str(Seq(l[:-1]).translate())
                put_new_record = False

with open(f"{path}/arctic_gene_atlas.json", "w") as handle:
    json.dump(clusters, handle, indent = 2, sort_keys=True)

with open(f"{path}/arctic_gene_atlas.json") as handle:
    clusters = json.load(handle)


with open(f"{path}/arctic_gene_atlas.ffn", "w") as handle:
    SeqIO.write([SeqRecord(id = f['cluster_id'], seq = Seq(f['representative_nucls']), description = f['hypo_func'].strip()) for f in tqdm(clusters)], handle, "fasta")

with open(f"{path}/arctic_gene_atlas.faa", "w") as handle:
    SeqIO.write([SeqRecord(id = f['cluster_id'], seq = Seq(f['representative_nucls']).translate(), description = f['hypo_func'].strip()) for f in tqdm(clusters)], handle, "fasta")

cds2bin = lambda cds : "_".join(cds.split("_")[:-1])

for clust in tqdm(clusters):
    clust['mOTU'] = set(bindat.loc[[cds2bin(cds) for cds in clust['members']],"mOTU"])

for clust in tqdm(clusters):
    clust['gtdbtk_classif'] = set(bindat.loc[[cds2bin(cds) for cds in clust['members']]].gtdbtk_classif)


"ln -s arctic_gene_atlas.ffn arctic_gene_atlas.fna"
"emapper.py --cpu 24 -m mmseqs --output arctic_gene_atlas -d bact -i arctic_gene_atlas.faa"

emapper_annots = {}
with open(f"{path}/arctic_gene_atlas.emapper.annotations") as handle:
    header = None
    counter = 0
    while not header and counter < 15:
        head = handle.readline().rstrip().split("\t")
        if head[0] == "#query" :
            header = head
        counter += 1
    for l in tqdm(handle):
        fields = l.rstrip().split("\t")
        emapper_annots[fields[0]] = {a : b for a,b in zip(header[1:], fields[1:])}


for clust in tqdm(clusters):
    if clust['cluster_id'] in emapper_annots:
        annot = emapper_annots[clust['cluster_id']]
        clust['root_eggNOG'] = annot['eggNOG_OGs'].split(",")[0]
        assert clust['root_eggNOG'].endswith("@1|root")
        clust['root_eggNOG'] = clust['root_eggNOG'].replace("@1|root","")
        clust['COG_category'] = annot['COG_category']
        clust['symbol'] = annot['Preferred_name']
        clust['KO'] = annot['KEGG_ko']
        clust['Description'] = annot['Description']
    else:
        clust['root_eggNOG'] = "NA"
        clust['COG_category'] = "NA"
        clust['symbol'] = "NA"
        clust['KO'] = "NA"
        clust['Description'] = "NA"

bindat = pandas.read_csv(pjoin(fpath,  binset_name + "_basics.csv" ), index_col=0)

classif = csv2dict(f"{temp_folder}/gtdbtk_r207/gtdbtk.ar53.summary.tsv", sep = "\t")
classif.update(csv2dict(f"{temp_folder}/gtdbtk_r207/gtdbtk.bac120.summary.tsv", sep = "\t"))

binset_stats = csv2dict(pjoin(root_folder, "binsets", binset_name, binset_name + "_basics2.csv"))
cleanz = {k : {'gtdbtk_classif_r207' : v['classification'], 'gtdbtk_notes_r207' : ";".join([ field + "=" + v[field].replace(" ","_") for field in ['note','classification_method','warnings'] if v[field] != "N/A"]), 'translation_table_r202' : v["translation_table"]} for k,v in classif.items()}

for k,v in cleanz.items():
    binset_stats[k].update(v)

#with open(pjoin(path,  "gene_clusters/gene2gene_clusters.json" )) as handle:
#    gene2gc = json.load(handle)

#clusters = {c['cluster_id'] : c for c in clusters}


#for clust in tqdm(clusters.values()):
#    if clust['representative'] in gene2gc:
#        clust['gene_cluster'] = gene2gc[clust['representative']]

#    else:
#        clust['root_eggNOG'] = "NA"
#        clust['COG_category'] = "NA"
#        clust['symbol'] = "NA"

#tt = [{gene2gc[vv] for vv in v['members'] if vv in gene2gc} for v in clusters.values()]

gene_md = pandas.DataFrame.from_records(clusters)
gene_md.index = gene_md.cluster_id
del gene_md['cluster_id']
gene_md[gene_md == "-"] = "NA"
gene_md.members = [";".join(l) for l in tqdm(gene_md.members)]
gene_md.to_csv(f"{path}/arctic_gene_atlas_basics.csv")
gene_md.to_csv(f"{fpath}/arctic_gene_atlas_basics.csv")


 nifHs = {"ERR594289.26_00811" , "SRR10912806.7_00626" , "MOSAIC-MIME-BINNING-ClusterAss-2_bin-406_03565" , "bulkbins-MOSAIC-MIME-BINNING-ClusterAss-1_bin-689_01844" , "bulkbins-MOSAIC-MIME-BINNING-2203_bin-08_00805" , "MOSAIC-MIME-BINNING-ClusterAss-1_bin-945_01309" , "bulkbins-MOSAIC-MIME-BINNING-2203_bin-60_30822" , "bulkbins-MOSAIC-MIME-BINNING-2203_bin-08_00823" , "MOSAIC-MIME-BINNING-ClusterAss-3a_bin-200_02370" , "bulkbins-MOSAIC-MIME-BINNING-2903_bin-087_01006"}


bin2gc = {"_".join(g.split("_")[:-1]) : set() for k, c in tqdm(gene_md.iterrows()) for g in c['members'].split(";")}

for k, c in tqdm(gene_md.iterrows()) :
    for g in c['members'].split(";"):
        bin2gc["_".join(g.split("_")[:-1])].add(c.name)

bin2egg = {"_".join(g.split("_")[:-1]) : set() for k, c in tqdm(gene_md.iterrows()) for g in c['members'].split(";")}

for c in tqdm(gene_md.iterrows()):
    for g in c[1]['members'].split(";"):
        if c[1]['root_eggNOG'] != "NA":
            bin2egg["_".join(g.split("_")[:-1])].add(c[1]['root_eggNOG'])

bin2ko ={"_".join(g.split("_")[:-1]) : set() for k, c in tqdm(gene_md.iterrows()) for g in c['members'].split(";")}

for c in tqdm(gene_md.iterrows()):
    for g in c[1]['members'].split(";"):
        if str(c[1]['KO']).startswith("ko:"):
            for ko in c[1]['KO'].split(","):
                bin2ko["_".join(g.split("_")[:-1])].add(ko)

for k,v in motupan_dat.items():
    bins = v.get('MAGs', "").split(";") + [] if v.get('SUBs', "") != v.get('SUBs', "") else v.get('SUBs', "").split(";")
    for vv in bins:
        if vv != "":
            bindat.loc[vv,'mOTU'] = k
            bindat.loc[vv,'representative'] = v['representative']

all_motus = {v : [] for v in bindat.mOTU if v != "NoMOTU"}
for k,v in zip(bindat.index, bindat.mOTU):
    if v  != "NoMOTU":
        all_motus[v] += [k]


os.makedirs(f"{path}/mOTUs", exist_ok = True)
for m,bins_ in tqdm(all_motus.items()):
    os.makedirs(f"{path}/mOTUs/{m}", exist_ok = True)
    os.makedirs(f"{path}/mOTUs/{m}/fnas", exist_ok = True)
    for b in bins_:
        shutil.copy(f"{path}/bins/{b}.fna", f"{path}/mOTUs/{m}/fnas")

title2log(f"regenerate checkm file", logfile)
for m,bins_ in tqdm(all_motus.items()):
    with open(f"{path}/mOTUs/{m}/checkm.txt", "w") as handle :
        handle.writelines(["Bin Id\tCompleteness\tContamination\n"] + [ f"{bin_}\t{bindat.loc[bin_,'percent_completion']}\t{bindat.loc[bin_,'percent_redundancy']}\n" for bin_ in bins_])

#threads = int(threads)
boots = 10
min_genomes = 10
remove_singleton_gcs = True
module_cutoff = 0.75
#temp_folder = pjoin(config_file['temp_folder'], "binsets", binset_name)
gc_folder = pjoin(temp_folder, "gene_clusters")
kegg_folder = pjoin(temp_folder, "KEGGs")
eggnog_folder = pjoin(temp_folder, "EGGNOGs")


title2log(f"Running mOTUpan on mOTUs >{min_genomes} genomes for eggNOG_OGs", logfile)
freetxt_line(f" and parsing mOTUpan-outs, removing singleton GCs : {remove_singleton_gcs}", logfile)

egg_pangenome_stats = {}
os.makedirs(f"{eggnog_folder}/mOTUs", exist_ok = True)
bindat['mOTUpan_completeness'] = None
for motu, bins_ in tqdm(all_motus.items()):
    bins_ = set(bins_)
    if len(bins_) > min_genomes:
        if motu not in egg_pangenome_stats:
#            with open(pjoin(temp_folder, "bin2gc.json") , "w") as handle:
#                json.dump({k : list(v) for k,v in bin2egg.items() if k in bins_}, handle, indent = 2, sort_keys = True)
#            call(f"mOTUpan.py --name {motu} --gene_clusters_file {temp_folder}/bin2gc.json  --length_  -o {eggnog_folder}/mOTUs/{motu}.tsv --boots {boots}  >>{logfile} 2>&1 ", shell = True)
            with open(f"{eggnog_folder}/mOTUs/{motu}.tsv") as handle:
                header = {}
                main = []
                for l in handle:
                    l = l.strip()
                    if l != "#" and not l.startswith("#mOTUlizer"):
                        if l.startswith("#"):
                            header[l.split("=")[0][1:]] = "=".join(l.split("=")[1:])
                        else :
                            main += [l]
                cols = main[0].split()
                gcs = [{k : ll for k,ll in zip(cols, l.split())} for l in main[1:] if l != '']
                if remove_singleton_gcs:
                    gcs = [d for d in gcs if int(d['genome_occurences']) > 1]
                core = []
                accessory = []
                for d in gcs:
                    if d['type'] == 'core':
                        core.append(d['trait_name'])
                    else :
                        accessory.append(d['trait_name'])
                        for dat in header['genomes'].split(";"):
                            dat = dat.split(":")
                            bin_ = dat[0]
                            new_comp = [dd for dd in dat[1:] if dd.startswith("posterior_complete=")][0]
                            bindat.loc[bin_,'mOTUpan_completeness'] = float(new_comp.replace("posterior_complete=",""))
                            egg_pangenome_stats[motu] = {
                            'genome_count' : int(header['genome_count']),
                            'core_length': int(header['core_length']),
                            'mean_prior_completeness' : float(header['mean_prior_completeness']),
                            'mean_posterior_completeness' : float(header['mean_posterior_completeness']),
                            "nb_bootstraps" : 10 if "bootstrapped_mean_false_positive_rate" not in header else header['bootstrapped_nb_reps'].count("<"),
                            "fpr" : 1 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_false_positive_rate'].split(";")[0]),
                            "sd_fpr" : 0 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_false_positive_rate'].split(";")[1].split("=")[1]),
                            "lowest_false" : 1 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_lowest_false_positive'].split(";")[0]),
                            "sd_lowest_false" : 0 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_lowest_false_positive'].split(";")[1].split("=")[1]),
                            'core' : ";".join(core),
                            'accessory' : ";".join(accessory)
                            }

kegg_pangenome_stats = {}
os.makedirs(f"{kegg_folder}/mOTUs", exist_ok = True)
for motu, bins_ in tqdm(all_motus.items()):
    bins_ = set(bins_)
    if len(bins_) > min_genomes:
        if motu not in kegg_pangenome_stats:
#            with open(pjoin(temp_folder, "bin2gc.json") , "w") as handle:
#                json.dump({k : list(v) for k,v in bin2ko.items() if k in bins_}, handle, indent = 2, sort_keys = True)
#            call(f"mOTUpan.py --name {motu} --gene_clusters_file {temp_folder}/bin2gc.json  --length_  -o {kegg_folder}/mOTUs/{motu}.tsv --boots {boots}  >>{logfile} 2>&1 ", shell = True)
            with open(f"{kegg_folder}/mOTUs/{motu}.tsv") as handle:
                header = {}
                main = []
                for l in handle:
                    l = l.strip()
                    if l != "#" and not l.startswith("#mOTUlizer"):
                        if l.startswith("#"):
                            header[l.split("=")[0][1:]] = "=".join(l.split("=")[1:])
                        else :
                            main += [l]
                cols = main[0].split()
                gcs = [{k : ll for k,ll in zip(cols, l.split())} for l in main[1:] if l != '']
                if remove_singleton_gcs:
                    gcs = [d for d in gcs if int(d['genome_occurences']) > 1]
                core = []
                accessory = []
                for d in gcs:
                    if d['type'] == 'core':
                        core.append(d['trait_name'])
                    else :
                        accessory.append(d['trait_name'])
                        for dat in header['genomes'].split(";"):
                            dat = dat.split(":")
                            bin_ = dat[0]
                            new_comp = [dd for dd in dat[1:] if dd.startswith("posterior_complete=")][0]
                            kegg_pangenome_stats[motu] = {
                            'genome_count' : int(header['genome_count']),
                            'core_length': int(header['core_length']),
                            'mean_prior_completeness' : float(header['mean_prior_completeness']),
                            'mean_posterior_completeness' : float(header['mean_posterior_completeness']),
                            "nb_bootstraps" : 10 if "bootstrapped_mean_false_positive_rate" not in header else header['bootstrapped_nb_reps'].count("<"),
                            "fpr" : 1 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_false_positive_rate'].split(";")[0]),
                            "sd_fpr" : 0 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_false_positive_rate'].split(";")[1].split("=")[1]),
                            "lowest_false" : 1 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_lowest_false_positive'].split(";")[0]),
                            "sd_lowest_false" : 0 if "bootstrapped_mean_false_positive_rate" not in header else float(header['bootstrapped_mean_lowest_false_positive'].split(";")[1].split("=")[1]),
                            'core' : ";".join(core),
                            'accessory' : ";".join(accessory)
                            }



module_path = pjoin(os.path.dirname(anvio.__file__), "data/misc/KEGG/modules/")
#title2log("Loading kegg_paths",logfile)
module2name = {}
module2def = {}
module2ko = {}

def s(keggs, def_line, andline = True):
    def_line = def_line.split() if andline else def_line.split(",")
    i = 0
    blocks = []
    fwd = None
    rev = None
    for v in def_line:
        blocks += [i]
        if v.startswith("("):
            fwd = v.count("(")
            rev = v.count(")")
        elif fwd:
            fwd += v.count("(")
            rev += v.count(")")
        if fwd == rev:
            i += 1
            fwd = None
            rev = None
    operation = " " if andline else ","
    big_operation = lambda x: sum(x)/len(x) if andline else max(x)
    and_block = [operation.join([b for i,b in enumerate(def_line) if blocks[i] == block]) for block in set(blocks)]
    and_block = [b[1:-1] if b.startswith("(") else b for b in and_block]
    return big_operation([s(keggs, b, not andline) if ("," in b or " " in b) else int((b in keggs) if not b.startswith("-") else True) for b in and_block])


for m in os.listdir(module_path):
    with open(pjoin(module_path, m)) as handle:
        lines = handle.readlines()
        for l in lines:
            if l.startswith("DEFINITION"):
                module2def[m] = " ".join(l.split()[1:])
            if l.startswith("NAME"):
                module2name[m] = " ".join(l.split()[1:])

motu2modules = {}
for motu, dd in tqdm(kegg_pangenome_stats.items()):
    kos = set(dd['core'].split(";"))
    completes = { module : s([k[3:] for k in kos], defline)for module, defline in module2def.items()}
    completes = { k : v for k, v in completes.items() if v >= module_cutoff}
    motu2modules[motu] = completes

motu_stats = {v['mOTU'] : v.to_dict() for k,v in  bindat.iterrows() if k == v['representative']}
for k in motu_stats:
    if k in egg_pangenome_stats:
        keys = ['core_length','nb_bootstraps', 'fpr', 'sd_fpr', 'lowest_false', 'sd_lowest_false', 'core', 'accessory']
        egg_stats = { "gc_" + kk : v for kk, v in egg_pangenome_stats[k].items() if kk in keys}
        kegg_stats = { "kegg_" + kk : v for kk, v in kegg_pangenome_stats[k].items() if kk in keys}
        motu_stats[k].update(egg_stats)
        motu_stats[k].update(kegg_stats)
        motu_stats[k]['ko_modules'] = ";".join(motu2modules[k])
