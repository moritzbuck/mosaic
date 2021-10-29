import json
from os.path import join as pjoin
from subprocess import call
import os
from tqdm import tqdm
from sourmash import MinHash
from sourmash.signature import SourmashSignature
from sourmash.signature import load_signatures
from multiprocessing import Pool
import gzip
import shutil
import re
import pandas

def load_sig(f):
    with gzip.open(f) as handle:
        ll = "".join( [l.decode() for l in handle.readlines()]).replace('\n','')
    return list(load_signatures(ll)) [0]

threash = 0.001
nb_matches = 10000

# our  minhash sigs
all_sub_sigs = [f"/home/moritz/data/M001_MIME/libraries/{s}/{s}.sig.gz" for s in os.listdir('/home/moritz/data/M001_MIME/libraries') if os.path.exists(f"/home/moritz/data/M001_MIME/libraries/{s}/{s}.sig.gz")]

# all SRAnsack sigs
all_ref_sigs = [f"data/libraries/{s}/{s}.sig.gz" for s in os.listdir('data/libraries') if os.path.exists(f"data/libraries/{s}/{s}.sig.gz")]

block_size = 50

# blockwise get minhash simi

sig_blocks = [all_ref_sigs[i:(i+block_size)] for i in list(range(0,len(all_ref_sigs), block_size))]
sub_sigs = {v.split("/")[-1][:-7] : load_sig(v) for v in tqdm(all_sub_sigs)}

all_2013_sigs = {s[:-7] : load_sig(f"/home/moritz/data/M001_MIME/other_data/sigs/{s}") for s in os.listdir("/home/moritz/data/M001_MIME/other_data/sigs/")}

def run_bloc2(i):
    if os.path.exists("blocks_odin/block_{}.csv".format(i)):
        return False
    ref_sigs = {v.split("/")[-1][:-7] : load_sig(v) for v in sig_blocks[i]}
    dists = {(k,l) : v.similarity(w, ignore_abundance=True) for k,v in all_2013_sigs.items() for l,w in ref_sigs.items()}
    dists = {k : v for k,v in dists.items() if v > threash and  k[0] != k[1]}
    print("Done bloc:", i)
    with open("blocks_odin/block_{}.csv".format(i), "w") as handle:
        handle.writelines(["{}\t{}\t{}\n".format(k[0],k[1],v) for k,v in dists.items()] )
    return True



pool2 = Pool(processes=24)
p2 = pool2.map_async(run_bloc2, list(range(len(sig_blocks))))

def run_bloc(i):
    if os.path.exists("blocks_mime/block_{}.csv".format(i)):
        return False
    ref_sigs = {v.split("/")[-1][:-7] : load_sig(v) for v in sig_blocks[i]}
    dists = {(k,l) : v.similarity(w, ignore_abundance=True) for k,v in sub_sigs.items() for l,w in ref_sigs.items()}
    dists = {k : v for k,v in dists.items() if v > threash and  k[0] != k[1]}
    print("Done bloc:", i)
    with open("blocks_mime/block_{}.csv".format(i), "w") as handle:
        handle.writelines(["{}\t{}\t{}\n".format(k[0],k[1],v) for k,v in dists.items()] )
    return True

pool = Pool(processes=24)
p = pool.map_async(run_bloc, list(range(len(sig_blocks))))



lines = []
for f in os.listdir("blocks_mime/"):
    with open("blocks_mime/" + f) as handle:
        lines += handle.readlines()
lines = [l.replace(".si","") for l in lines]
# load the simis from the saved file, because
lines = [l.split("\t") for l in lines]
pairs = {(l[0], l[1]) : float(l[2]) for l in lines}

lines2 = []
for f in os.listdir("blocks_odin/"):
    with open("blocks_odin/" + f) as handle:
        lines2 += handle.readlines()

# load the simis from the saved file, because
lines2 = [l.split("\t") for l in lines2]
pairs2 = {(l[0], l[1]) : float(l[2]) for l in lines2}



# mainly cleaning up GPS coords for the context

def preproc_sra_md(sra_md):
    for k,v in sra_md.items():
        if v.get('sample_attributes'):
            if v['sample_attributes'] == v['sample_attributes']:
                v.update(v['sample_attributes'])
    sra_md = pandas.DataFrame.from_dict({v['SRA_ID'] : v for k,v in sra_md.items()}, orient = "index")

    del sra_md['TARGETED_LOCI']
    del sra_md['LIBRARY_CONSTRUCTION_PROTOCOL']
    del sra_md['study_id']
    del sra_md['LIBRARY_SOURCE']
    del sra_md['LIBRARY_STRATEGY']
    del sra_md['Platform']
    del sra_md['SRA_ID']
    del sra_md['Statistics']
    del sra_md['sample_attributes']
    del sra_md['sample_id']
    del sra_md['sample_name']
    sra_md.loc['ERR3589581', 'Latitude End'] = "71.0704"
    sra_md.loc['ERR3589581', 'Longitude End'] = "174.9916"

    sra_md['coord'] = ["NA" if v in ['missing', 'not collected', "Missing", "N/A", 'not applicable','Not Applicable', 'Not applicable', "NULL", "Unknown"] else v  for v in sra_md.lat_lon]

    def tolat_lon(lat, lon):
        if lat >= 0 :
            lat = str(lat) + " N "
        else :
            lat = str(-lat) + " S "
        if lon >= 0 :
            lon = str(lon) + " E"
        else :
            lon = str(-lon) + " W"
        return lat + lon


    is_nan = lambda l : l != l
    for k in sra_md.index:
        if not is_nan(sra_md.loc[k]['Latitude Start']) and is_nan(sra_md.loc[k]['Latitude Start']):
            sra_md['coord'][k] = tolat_lon(float(isnan(sra_md.loc[k]['Latitude Start'])), float(isnan(sra_md.loc[k]['Longitude Start'])))

    for k in sra_md.index:
        if not is_nan(sra_md.loc[k]['Latitude Start']) and not is_nan(sra_md.loc[k]['Latitude End']):
            lat = (float(sra_md.loc[k]['Latitude Start']) + float(sra_md.loc[k]['Latitude End']))/2
            lon = (float(sra_md.loc[k]['Longitude Start']) + float(sra_md.loc[k]['Longitude End']))/2
            sra_md['coord'][k] = tolat_lon(lat, lon)

    clean_char = set("0123456789,.-")
    to_dec = lambda v : v[0] + v[1]/60 + v[2]/3600

    for k in sra_md.index:
        if  not is_nan(sra_md.loc[k]['geographic location (latitude)']) and  not is_nan(sra_md.loc[k]['geographic location (longitude)']):
            lat =  sra_md.loc[k]['geographic location (latitude)'].replace(" DD", "")
            lon =  sra_md.loc[k]['geographic location (longitude)'].replace(" DD", "")
            if all([c in clean_char for c in lat + lon]):
                lat = float(lat.replace(",", "."))
                lon = float(lon.replace(",", "."))
                sra_md['coord'][k] = tolat_lon(lat, lon)
            else :
                lat_tri = [float(l) for l in re.split("[^0-9.,]",lat) if l != ""]
                lon_tri = [float(l) for l in re.split("[^0-9.,]",lon) if l != ""]
                if len(lat_tri) == 3 and len(lon_tri) == 3:
                    if ("N" in lat or "S" in lat) and ("E" in lon or "W" in lon):
                        lat_str = str(to_dec(lat_tri)) + (" N " if "N" in lat else " S ")
                        lon_str = str(to_dec(lon_tri)) + (" E " if "E" in lat else " W ")
                        sra_md['coord'][k] = lat_str + lon_str
                elif k.startswith("ERR133318"):
                    lon = lon + " E"
                    sra_md['coord'][k] = lat + " " + lon
                elif k == 'ERR1299101':
                    lon = lon.replace("'", "")
                    sra_md['coord'][k] = tolat_lon(float(lat), float(lon))


    sra_md = sra_md.fillna("NA")

    for id, lat_lon in sra_md.coord[sra_md.coord.apply(lambda x : len(str(x).split()) != 4 and x != 'NA')].items():
        lat_lon_tri = [float(l) for l in re.split("[^0-9.]",lat_lon) if l != ""]
        print("none decimal e.g. 50032.74 50008.79 coords")
        if len(lat_lon_tri) ==6:
            lat_tri = lat_lon_tri[0:3]
            lon_tri = lat_lon_tri[3:]
        else :
            lat_tri = lat_lon_tri[0:2] + [0]
            lon_tri = lat_lon_tri[2:] + [0]
        if ("N" in lat_lon or "S" in lat_lon) and ("E" in lat_lon or "W" in lat_lon):
            lat_str = str(to_dec(lat_tri)) + (" N " if "N" in lat_lon else " S ")
            lon_str = str(to_dec(lon_tri)) + (" E " if "E" in lat_lon else " W ")
            sra_md['coord'][id] = lat_str + lon_str
        else :
            sra_md['coord'][id] = tolat_lon(to_dec(lat_tri), to_dec(lon_tri))



    lat_long2lat = lambda x: (float(x.split()[0]) * (-1 if "S" in x else 1) ) if x != "NA" else x
    lat_long2lon = lambda x: (float(x.split()[2]) * (-1 if "W" in x else 1) ) if x != "NA" else x

    sra_md['lat'] = sra_md.coord.apply(lat_long2lat)
    sra_md['lon'] = sra_md.coord.apply(lat_long2lon)

    id2tax = {v[1] : v[0]  for k,v in sra_md[['taxon', 'taxon_id']].iterrows() if v[0] != ""}
    sra_md['taxon'] = [id2tax[v] for v in sra_md.taxon_id]

    bad_study = sra_md.loc['ERR4193663']['study']

    lats = sra_md.loc[sra_md.study == bad_study, 'lat']
    lons = sra_md.loc[sra_md.study == bad_study, 'lon']
    sra_md.loc[sra_md.study == bad_study, 'lat'] = lons
    sra_md.loc[sra_md.study == bad_study, 'lon'] = lats
    sra_md.loc[sra_md.study == bad_study, 'coord'] =  [tolat_lon(float(b),float(a)) for a, b in zip(lats, lons)]
    return sra_md


# Load SRAnsack metadata

with open("/home/moritz/data/SRAprov/data/dbs/sra_data.r1.json") as handle :
     sra_md = json.load(handle)

with open("/home/moritz/data/SRAprov/data/dbs/sra_data.json") as handle :
     sra_md.update(json.load(handle))
sra_md = {v['SRA_ID'] : v for k,v in sra_md.items()}

# cleanup
sra_md = preproc_sra_md(sra_md)

# extract  metagenomic contexts
cntxt = {k[1] for k,v in pairs.items() if v > 0.1}
broad_cntxt = {k[1] for k,v in pairs.items() if v > 0.05}

odin_cntxt = {k[1] for k,v in pairs2.items() if v > 0.1}
odin_broad_cntxt = {k[1] for k,v in pairs2.items() if v > 0.05}


sra_md = sra_md.loc[odin_broad_cntxt.union(broad_cntxt)]
sra_md['context'] = ["close_context" if i in cntxt or i in odin_cntxt else "broad_context" for i in sra_md.index]
#sra_md = sra_md[[c  for c in sra_md.columns if sum(sra_md[c] == "NA") < 208]]

# some more cleanup
uniq_cols = {c.lower() : [] for c in sra_md.columns}
for c in sra_md.columns:
    uniq_cols[c.lower()] += [c]

uniq_cols["collection_date"] = ["collection_date", "collection date"]
uniq_cols["depth"] = uniq_cols.get("depth", []) + ["mean composite depth", "Depth (m)", "water depth (m)", "sampling depth", "sample_depth", "depth_sample", "water_depth"]
uniq_cols["temp"] = uniq_cols.get("temp", []) + ["temperature", 'temp (oC)', 'temp, degrees C', 'Temperature (C?)', 'water_temp', 'temperature_C', 'Temperature_degree_C', 'temperature_deg_c', 'temperature (celsius)']

to_vacate = ["NA" , "Illumina HiSeq 4000 paired end sequencing", "Illumina HiSeq 3000 paired end sequencing", 'Illumina HiSeq 2000 paired end sequencing; marine metagenome_Illumina_F']
for k,v in uniq_cols.items():
    if len(v) > 1:
        for i in sra_md.index:
            gg = list(set([vv for vv in sra_md.loc[i,v] if vv not in to_vacate]))
            assert len(gg) < 2, gg
            if len(gg) == 1:
                sra_md.loc[i,k] = gg[0]
        for vv in v:
            if vv != k:
                del sra_md[vv]
sra_md.columns = [c.lower() for c in sra_md.columns]

sra_md['collection_year'] = [";".join({l for l in re.split(r"[-/ ]", d) if len(l) == 4}) for d in sra_md['collection_date']]


#sra_md[sorted(sra_md.columns)].to_csv("/home/moritz/data/M001_MIME/dbs/context_md_raw.csv")
#keep_cols  = pandas.read_csv("/home/moritz/data/M001_MIME/dbs/context_md_pruned.csv", index_col=0).columns
keep_cols = ['ammonium', 'chlorophyll', 'chlorophyll sensor', 'collection_date',
       'collection_year', 'context', 'coord', 'depth', 'diss_oxygen',
       'environment (feature)', 'environment (feature) further details',
       'environment (material)', 'fluor', 'geo_loc_name', 'lat_lon',
       'nb_bases', 'nb_reads', 'nitrate', 'nitrate sensor', 'nitrite',
       'oxygen sensor', 'phosphate', 'project name', 'protocol description',
       'salinity', 'salinity sensor', 'sampling platform', 'sampling site',
       'sampling station', 'silicate', 'size fraction lower threshold',
       'size fraction upper threshold', 'study', 'temp',
       'title', 'tot_nitro', 'lat', 'lon']

sra_md[keep_cols].to_csv("/home/moritz/data/M001_MIME/dbs/context_md.csv")


# load manualy fixed md
sra_md = pandas.read_csv("/home/moritz/data/M001_MIME/dbs/context_md_with_manual.csv", index_col=0)

all_context_sigs = {s : f"data/libraries/{s}/{s}.sig.gz" for s in cntxt.union(broad_cntxt)}

#for s_id , sig in tqdm(sub_sigs.items()) :
#    if s_id.startswith("MOSAIC"):
#        sub_sigs[s_id] = load_sig(sig)

sub_sigs.update(all_2013_sigs)

# updateing sigs collection with the context
for s_id , sig in tqdm(all_context_sigs.items()) :
    sub_sigs[s_id] = load_sig(sig)

# computing distance matrix (with abundances)
dists = {l : {k : v.similarity(w) for k,v in sub_sigs.items()} for l,w in tqdm(sub_sigs.items())}

for k,v in tqdm(sub_sigs.items()):
    for l, w in sub_sigs.items():
        if l not in dists:
            dists[l] = dict()
        if k not in dists[l]:
            dists[l][k] = v.similarity(w)

dists = pandas.DataFrame.from_dict(dists).fillna(0)
#dists.index = [i if i not in sra_md.index else sra_md['short_name'][i] for i in dists.index]
dists.index = ["ODIN-2013-" + i if i.startswith("Brine") or i.startswith("Seawater") else i for i in dists.index]

# dump data into files
dists.to_csv("/home/moritz/data/M001_MIME/dbs/samples_nd_context_all_simis_w_abundance.csv", index_label = 'short_name')

with open("/home/moritz/data/M001_MIME/dbs/SRAnsack_context.txt", "w") as handle:
    handle.writelines(["\n".join(cntxt)])
with open("/home/moritz/data/M001_MIME/dbs/SRAnsack_broad_context.txt", "w") as handle:
    handle.writelines(["\n".join(broad_cntxt)])
context_folder = "/home/moritz/data/M001_MIME/dbs/context_bins/"
sransack_loc = "/home/moritz/data/SRAprov/data/libraries/"
for f in cntxt.union(broad_cntxt):
    for bin_ in os.listdir(pjoin(sransack_loc, f, "bins")):
        shutil.copy(pjoin(sransack_loc, f, "bins", bin_), context_folder)
call(f"unpigz {context_folder}/*.gz", shell=True)
for f in os.listdir(context_folder):
    shutil.move(pjoin(context_folder, f), pjoin(context_folder, f.replace(".fa", ".fna")))

with open("/home/moritz/projects/mosaic/M001_MIME/config_nd_metadata/M001_sample_data.json") as handle:
    mosaic_md = json.load(handle)
