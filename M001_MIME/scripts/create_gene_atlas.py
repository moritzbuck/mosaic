import os
from tqdm import tqdm
from os.path import join as pjoin

seqid = 0.95
cov = 0.8
covmode = 1
temp_folder = "/home/moritz/temp"
threads = 24
path = "/home/moritz/data/M001_MIME/binsets/TEMP-COMPLETE-SET/"
if not os.path.exists(pjoin(path,  "arctic_gene_atlas.ffn"))):
    for file in tqdm(os.listdir(pjoin(path, "bins"))):
        call("cat {file} >> {ass}".format(file = pjoin(path, "bins", file, file +".ffn"), ass = pjoin(path,  "arctic_gene_atlas.ffn")), shell=True)


call(f"mmseqs easy-cluster --min-seq-id {seqid} --cov-mode {covmode} -c {cov} --threads {threads} {path}/arctic_gene_atlas.ffn {path}/arctic_gene_atlas {temp_folder}/mmseqs_temp > {logfile}", shell=True)
