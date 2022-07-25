script , binset_name, config_file, root_folder, out_folder, threads = "" , "COMPLETE-SET-AND-BULK", "/home/moritz/projects/mosaic/M001_MIME/config_nd_metadata/M001_config.json", "/home/moritz/data/M001_MIME/", "/home/moritz/data/M001_MIME/binsets/COMPLETE-SET-AND-BULK/", 23
mapping_name = "MOSAIC-MIME-RNA"

import sys, os
from tqdm import tqdm
from os.path import join as pjoin
sys.path.append(os.getcwd())

from workflow.scripts.utils import generate_config, title2log, freetxt_line, dict2file, gff2anvio, csv2dict
import shutil
from subprocess import call
from os.path import basename
import json
from Bio import SeqIO
import re
#from anvio.summarizer import ContigSummarizer
from tempfile import NamedTemporaryFile
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
#from anvio.dbops import ContigsSuperclass
#from anvio.utils import export_sequences_from_contigs_db
import pandas
import json

config_file = generate_config(config_file)

mapping_path = f"/home/moritz/data/M001_MIME/mappings/{mapping_name}/"
temp_folder = pjoin(config_file['temp_folder'], "binsets", binset_name)

bindat = pandas.read_csv(pjoin(out_folder,  binset_name + "_basics.csv" ), index_col=0, sep =",", low_memory = False)
gene_md = pandas.read_csv(f"{out_folder}/arctic_gene_atlas_basics.csv", index_col = 0)
motu_md = pandas.read_csv(f"{out_folder}/motustats.csv", index_col = 0)
tpm = pandas.read_csv(f"{mapping_path}/tpm.csv", index_col = 0)

id2gene_cluster = {vv : k for k,v in tqdm(gene_md.members.to_dict().items()) for vv in v.split(";")}

####### making function to bin mappings #######
weird =  []
bin2ko = { bin_ : set() for bin_ in bindat.index}
for genes, kos in tqdm(zip(gene_md.members, gene_md.KO)):
    for gene in genes.split(";"):
        if "_".join(gene.split("_")[:-1]) in bin2ko:
            if kos == kos:
                for ko in kos.split(","):
                    bin2ko["_".join(gene.split("_")[:-1])].add(ko)
        else :
            weird += ["_".join(gene.split("_")[:-1])]

ko2bin = {kk : set() for k in bin2ko.values() for kk in k}
for bin_, kos in tqdm(bin2ko.items()):
    for ko in kos:
        ko2bin[ko].add(bin_)

with open(f"{out_folder}/COMPLETE-SET-AND-BULK.gff") as handle:
    prokka_annots = [l for l in tqdm(handle) if 'CDS' in l]

prokka_annots = [ dict([( pp.split("=")[0], pp.split("=")[1]) for pp in p.split("\t")[8].split(";")])  for p in tqdm(prokka_annots)]
id2name = {v['ID'] : v['Name'].split("_")[0] for v in tqdm(prokka_annots) if "Name" in v}
name2id = {v : [] for v in id2name.values()}
for k,v in tqdm(id2name.items()):
    name2id[v] += [k]

name2gene_cluster = { k : {id2gene_cluster[vv] for vv in v} for k,v in tqdm(name2id.items())}

####### checking out nifHs / nifKs ######

# writing candidate genes
put_nifK = set(gene_md.index[[False if x != x else "ko:K02591" in x for x in gene_md.KO]])
put_nifH = set(gene_md.index[[False if x != x else "ko:K02588" in x for x in gene_md.KO]])
put_nifD = set(gene_md.index[[False if x != x else "ko:K02586" in x for x in gene_md.KO]])
os.makedirs(f"{temp_folder}/genes/nifH/", exist_ok = True)

with open(f"{temp_folder}/genes/nifH/put_nifH.faa", "w") as handle:
    handle.writelines([f">eggnog__{k}\n{v['representative_aas']}\n" for k, v in gene_md.loc[set(put_nifH).difference(name2gene_cluster['nifH'])].iterrows()])
    handle.writelines([f">prokka_{k}\n{v['representative_aas']}\n" for k, v in gene_md.loc[name2gene_cluster['nifH'].difference(put_nifH)].iterrows()])
    handle.writelines([f">both_{k}\n{v['representative_aas']}\n" for k, v in gene_md.loc[name2gene_cluster['nifH'].intersection(put_nifH)].iterrows()])

with open(f"{temp_folder}/put_nifK.faa", "w") as handle:
    handle.writelines([f">{k}\n{v['representative_aas']}\n" for k, v in gene_md.loc[put_nifK].iterrows()])

with open(f"{temp_folder}/put_nifD.faa", "w") as handle:
    handle.writelines([f">{k}\n{v['representative_aas']}\n" for k, v in gene_md.loc[put_nifD].iterrows()])

# loading db and writing them with candidates into a bmff
nifH_ali_seqs = []
for s in SeqIO.parse("/home/moritz/dbs/proteins/nifH/nifH_ingroup.faa", "fasta"):
    org_name = [t for t in s.description.split("=")]
    org_name = "_".join([t for i,t in enumerate(org_name[1:]) if org_name[i].endswith("OS")][0].split(" ")[:-1])
    org_name = org_name.split("_(")[0]
    id_ = s.id.split("|")[1]
    s.id = f"nifH_{org_name}_{id_}"
    s.description = ""
    nifH_ali_seqs += [s]

for s in SeqIO.parse("/home/moritz/dbs/proteins/nifH/nifH_outgroup.faa", "fasta"):
    org_name = [t for t in s.description.split("=")]
    org_name = "_".join([t for i,t in enumerate(org_name[1:]) if org_name[i].endswith("OS")][0].split(" ")[:-1])
    org_name = org_name.split("_(")[0]
    id_ = s.id.split("|")[1]
    s.id = f"outgroup_{org_name}_{id_}"
    s.description = ""
    nifH_ali_seqs += [s]

for s in SeqIO.parse(f"{temp_folder}/genes/nifH/put_nifH.faa", "fasta"):
    nifH_ali_seqs += [s]

with open(f"{temp_folder}/genes/nifH/nifH_w_db.faa", "w") as handle:
    SeqIO.write(nifH_ali_seqs, handle, "fasta")

call(f"kalign {temp_folder}/genes/nifH/nifH_w_db.faa > {temp_folder}/genes/nifH/nifH_w_db.ali.faa", shell = True)
call(f"fasttree < {temp_folder}/genes/nifH/nifH_w_db.ali.faa > {temp_folder}/genes/nifH/nifH_w_db.tree", shell = True)

# remove highly probable noin-nifH based on tree

not_nifHs = {"arctic_gene_atlas-5688836", "arctic_gene_atlas-745799", "arctic_gene_atlas-3769438"}
put_nifH = put_nifH.difference(not_nifHs)

# matching up nifHs to possible nifKs
got_nifK = {"_".join(vv.split("_")[:-1]) for v in gene_md.loc[put_nifK].members for vv in v.split(";")}
got_nifD = {"_".join(vv.split("_")[:-1]) for v in gene_md.loc[put_nifD].members for vv in v.split(";")}

put_nifH = list(name2gene_cluster['nifH'].union(put_nifH))
matches = { k : sum([ "_".join(vv.split("_")[:-1]) in got_nifK for vv in v.split(";")])/len(v.split(";")) for k,v in zip(put_nifH, gene_md.loc[put_nifH].members)}
matches = [k for k,v in matches.items() if v > 0 and k != "arctic_gene_atlas-6292501"]
nif_taxas = { k : { bindat.loc["_".join(vv.split("_")[:-1]),"gtdbtk_classif_r207"] for vv in v.split(";") if bindat.loc["_".join(vv.split("_")[:-1]),"gtdbtk_classif_r207"] == bindat.loc["_".join(vv.split("_")[:-1]),"gtdbtk_classif_r207"] } for k,v in zip(matches, gene_md.loc[matches].members)}
nif_exps = tpm.loc[nif_taxas.keys()]
{ bindat.loc["_".join(vv.split("_")[:-1]),"gtdbtk_classif_r207"] for vv in v.split(";") if bindat.loc["_".join(vv.split("_")[:-1]),"gtdbtk_classif_r207"] == bindat.loc["_".join(vv.split("_")[:-1]),"gtdbtk_classif_r207"] }
#the only candidate nifH that has no nifK but seams clearly a nifH, weird-ish
hmmmm = "arctic_gene_atlas-2561075"

key_fixers_ = "d__Bacteria;p__Myxococcota;c__UBA9042;o__JABWCM01;f__JABWCM01"


# Other nitrogen processes

nitrifiers = ko2bin['ko:K10944'].intersection(ko2bin['ko:K10945']).intersection(ko2bin['ko:K10946'])

nars = ko2bin['ko:K00370'].intersection(ko2bin['ko:K00371']).intersection(ko2bin['ko:K00374'])
naps = ko2bin['ko:K02567'].intersection(ko2bin['ko:K02568'])
nirKSs = ko2bin['ko:K00368'].union(ko2bin['ko:K15864'])
norBCs = ko2bin['ko:K04561'].intersection(ko2bin['ko:K02305'])
nosZs = ko2bin['ko:K00376']
denitrifiers = (nars.union(naps)).intersection(nirKSs).intersection(norBCs).intersection(nosZs)

# nosZ

egg_nosZ = set(gene_md.index[[False if x != x else "ko:K00376" in x for x in gene_md.KO]])
prok_nosZ = name2gene_cluster['nosZ']


os.makedirs(f"{temp_folder}/genes/nosZ/", exist_ok = True)

with open(f"{temp_folder}/genes/nosZ/put_nosZ.faa", "w") as handle:
    handle.writelines([f">eggnog__{k}\n{v['representative_aas']}\n" for k, v in gene_md.loc[set(egg_nosZ).difference(prok_nosZ)].iterrows()])
    handle.writelines([f">prokka_{k}\n{v['representative_aas']}\n" for k, v in gene_md.loc[prok_nosZ.difference(egg_nosZ)].iterrows()])
    handle.writelines([f">both_{k}\n{v['representative_aas']}\n" for k, v in gene_md.loc[prok_nosZ.intersection(egg_nosZ)].iterrows()])

nosZ_ali_seqs = []
for s in SeqIO.parse("/home/moritz/dbs/proteins/nosZ/curated_nosZ.faa", "fasta"):
    if "OS=" in  s.description:
        org_name = [t for t in s.description.split("=")]
        org_name = "_".join([t for i,t in enumerate(org_name[1:]) if org_name[i].endswith("OS")][0].split(" ")[:-1])
        org_name = org_name.split("_(")[0]
    else :
        org_name = "NA"
    id_ = s.id.split("|")[1] if "|" in s.id else s.id
    s.id = f"nosZ_{org_name}_{id_}"
    s.description = ""
    nosZ_ali_seqs += [s]

with open("/home/moritz/dbs/proteins/nosZ/nosZ_clade_II.txt") as handle:
    cladeII = {l.strip() for l in handle}
for s in SeqIO.parse("/home/moritz/dbs/proteins/nosZ/outgroup_nosZ.faa", "fasta"):
    if "OS=" in  s.description:
        org_name = [t for t in s.description.split("=")]
        org_name = "_".join([t for i,t in enumerate(org_name[1:]) if org_name[i].endswith("OS")][0].split(" ")[:-1])
        org_name = org_name.split("_(")[0]
    else :
        org_name = "NA"
    id_ = s.id.split("|")[1] if "|" in s.id else s.id
    if id_ in cladeII:
        s.id = f"nosZ_cladeII_{org_name}_{id_}"
    else :
        s.id = f"outgroup_{org_name}_{id_}"

    s.description = ""
    nosZ_ali_seqs += [s]

for s in SeqIO.parse(f"{temp_folder}/genes/nosZ/put_nosZ.faa", "fasta"):
    nosZ_ali_seqs += [s]

with open(f"{temp_folder}/genes/nosZ/nosZ_w_db.faa", "w") as handle:
    SeqIO.write(nosZ_ali_seqs, handle, "fasta")

call(f"kalign {temp_folder}/genes/nosZ/nosZ_w_db.faa > {temp_folder}/genes/nosZ/nosZ_w_db.ali.faa", shell = True)
call(f"fasttree < {temp_folder}/genes/nosZ/nosZ_w_db.ali.faa > {temp_folder}/genes/nosZ/nosZ_w_db.tree", shell = True)


# nrfA

egg_nrfA = set(gene_md.index[[False if x != x else "ko:K03385" in x for x in gene_md.KO]])
prok_nrfA = name2gene_cluster['nrfA']


os.makedirs(f"{temp_folder}/genes/nrfA/", exist_ok = True)

with open(f"{temp_folder}/genes/nrfA/put_nrfA.faa", "w") as handle:
    handle.writelines([f">eggnog__{k}\n{v['representative_aas']}\n" for k, v in gene_md.loc[set(egg_nrfA).difference(prok_nrfA)].iterrows()])
    handle.writelines([f">prokka_{k}\n{v['representative_aas']}\n" for k, v in gene_md.loc[prok_nrfA.difference(egg_nrfA)].iterrows()])
    handle.writelines([f">both_{k}\n{v['representative_aas']}\n" for k, v in gene_md.loc[prok_nrfA.intersection(egg_nrfA)].iterrows()])

nrfA_ali_seqs = []
for s in SeqIO.parse("/home/moritz/dbs/proteins/nrfA/curated_nrfAs.faa", "fasta"):
    if "OS=" in  s.description:
        org_name = [t for t in s.description.split("=")]
        org_name = "_".join([t for i,t in enumerate(org_name[1:]) if org_name[i].endswith("OS")][0].split(" ")[:-1])
        org_name = org_name.split("_(")[0]
    else :
        org_name = "NA"
    id_ = s.id.replace("|", "_")
    s.id = f"nrfA_{s.id}"
    s.description = ""
    nrfA_ali_seqs += [s]

for s in SeqIO.parse("/home/moritz/dbs/proteins/nrfA/outgroup_nrfAs.faa", "fasta"):
    if "OS=" in  s.description:
        org_name = [t for t in s.description.split("=")]
        org_name = "_".join([t for i,t in enumerate(org_name[1:]) if org_name[i].endswith("OS")][0].split(" ")[:-1])
        org_name = org_name.split("_(")[0]
    else :
        org_name = "NA"
    id_ = s.id.replace("|", "_")
    s.id = f"outgroup_{org_name}_{s.id}"

    s.description = ""
    nrfA_ali_seqs += [s]

for s in SeqIO.parse(f"{temp_folder}/genes/nrfA/put_nrfA.faa", "fasta"):
    nrfA_ali_seqs += [s]

with open(f"{temp_folder}/genes/nrfA/nrfA_w_db.faa", "w") as handle:
    SeqIO.write(nrfA_ali_seqs, handle, "fasta")

call(f"kalign {temp_folder}/genes/nrfA/nrfA_w_db.faa > {temp_folder}/genes/nrfA/nrfA_w_db.ali.faa", shell = True)
call(f"fasttree < {temp_folder}/genes/nrfA/nrfA_w_db.ali.faa > {temp_folder}/genes/nrfA/nrfA_w_db.tree", shell = True)
{"arctic_gene_atlas-542595", "arctic_gene_atlas-2357413", "arctic_gene_atlas-745799", "arctic_gene_atlas-3038205", "arctic_gene_atlas-3602825", "arctic_gene_atlas-3073417", "arctic_gene_atlas-2201712", "arctic_gene_atlas-4252661", "arctic_gene_atlas-4103301", "arctic_gene_atlas-5367056", "arctic_gene_atlas-4559447", "arctic_gene_atlas-1486204", "arctic_gene_atlas-5688836", "arctic_gene_atlas-2193172", "arctic_gene_atlas-362020", "arctic_gene_atlas-3268297", "arctic_gene_atlas-4483513", "arctic_gene_atlas-2797851", "arctic_gene_atlas-5388019", "arctic_gene_atlas-1833927", "arctic_gene_atlas-682943", "arctic_gene_atlas-979869", "arctic_gene_atlas-3044705", "arctic_gene_atlas-5792431", "arctic_gene_atlas-2561075", "arctic_gene_atlas-5821544", "arctic_gene_atlas-4342488", "arctic_gene_atlas-333240", "arctic_gene_atlas-1711838", "arctic_gene_atlas-2034977", "arctic_gene_atlas-5170627", "arctic_gene_atlas-3475739", "arctic_gene_atlas-3769438", "arctic_gene_atlas-4994843}
