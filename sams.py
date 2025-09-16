import samap
import scanpy as sc
import samalg
import anndata as ad
import pandas as pd

# -------------------------------------------------------- #
# TO DO: edit ID names in human h5ad to be ensembl ids
import mygene
#fn2 = '/Users/dr/cr/pr/class_II_project/in/human_raw_wlabels.h5ad'
#fn2 = ad.read_h5ad(fn2)
# print("var columns:", list(fn2.var.columns))
# print("first var_names:", fn2.var_names[:10].tolist())

symbols = pd.Index(fn2.var_names.astype(str))
mg = mygene.MyGeneInfo()
res = mg.querymany(symbols.tolist(),
                   scopes="symbol",
                   fields="ensembl.gene",
                   species="human",
                   as_dataframe=True).reset_index()

# get genes from mygene output
def flat(x):
    if isinstance(x, list): x = x[0]
    if isinstance(x, dict): return x.get("gene")
    return x

res["ENSG"] = res["ensembl.gene"].map(flat)
mapper = (res.dropna(subset=["ENSG"])
            .drop_duplicates("query")
            .set_index("query")["ENSG"])
new_ids = pd.Series(symbols.map(mapper), index=symbols, dtype="string")
keep = new_ids.notna().to_numpy()
fn2 = fn2[:, keep].copy()
fn2.var["ensembl_id"] = new_ids[keep].to_numpy()
fn2.var_names = fn2.var["ensembl_id"]
fn2.var_names_make_unique()

if fn2.var.index.name in fn2.var.columns:
    fn2.var.rename(columns={fn2.var.index.name: fn2.var.index.name + "_col"}, inplace=True)
fn2.write('/Users/dr/cr/pr/class_II_project/in/human_raw_wlabels_new.h5ad')

# -------------------------------------------------------- #
# load datasets
fn1 = '/Users/dr/cr/pr/class_II_project/in/mouse_raw_wlabels.h5ad'
fn2 = '/Users/dr/cr/pr/class_II_project/in/human_raw_wlabels_new.h5ad'

filenames = {'mm':fn1,'hs':fn2}

# -------------------------------------------------------- #
from samalg import SAM
import sys

# make sam object
sam1=SAM()
sam1.load_data(fn1)

sam2=SAM()
sam2.load_data(fn2)

sams = {'mm':sam1,'hs':sam2}

# -------------------------------------------------------- #
from samap.mapping import SAMAP
# TO DO: move maps into ~/in/maps

sam1.run()  # optionally sam1.run(k=20, npcs=50, standardize=True)
sam2.run()

# -------------------------------------------------------- #
# IMPORTANT CHECK: make sure maps are in subfolders called "mmhs" and "hsmm" or will get error
# CHECK 2: make sure whatever gene naming convention you use in the map is the same as var names
sm = SAMAP(
        sams,
        f_maps = '/Users/dr/cr/pr/class_II_project/in/maps/',
    )

# -------------------------------------------------------- #
sm.run(neigh_from_keys=["cell_type"], umap = True)
