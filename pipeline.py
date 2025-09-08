# -------------------------------------------------------- #
# loading packages
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# -------------------------------------------------------- #
# loading files
spleen1 = sc.read_10x_mtx(
    "dr/cr/pr/class_II_project/in/human/", # directory
    var_names="gene_symbols", # use symbol genes for homologene and also cell typist
    cache=True,
)
spleen2 = sc.read_10x_mtx(
    "dr/cr/pr/class_II_project/in/mouse/",
    var_names="gene_symbols",
    cache=True,
)
#spleen = spleen1
#spleen = spleen2
spleen.var_names_make_unique()


# -------------------------------------------------------- #
# quality control
sc.pp.filter_cells(spleen, min_genes=250)
sc.pp.filter_genes(spleen, min_cells=10)


# mitochondrial genes, "MT-" for human, "Mt-" for mouse
spleen.var["mt"] = spleen.var_names.str.startswith("MT-")
#spleen.var["mt"] = spleen.var_names.str.startswith("Mt-")
# ribosomal genes
spleen.var["ribo"] = spleen.var_names.str.startswith(("RPS", "RPL"))
#spleen.var["ribo"] = spleen.var_names.str.startswith(("Rps", "Rpl"))
# hemoglobin genes
spleen.var["hb"] = spleen.var_names.str.contains("^HB[^(P)]")
#spleen.var["hb"] = spleen.var_names.str.contains("^Hb[^(p)]")

sc.pp.calculate_qc_metrics(
    spleen, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

spleen.obs["log10_genes_per_umi"] = (
    np.log10(spleen.obs["n_genes_by_counts"] + eps) / np.log10(spleen.obs["total_counts"] + eps)
)

#spleen.obs["log10_genes_per_umi"] = np.log10(spleen.obs["n_genes_by_counts"].clip(lower=1) / spleen.obs["total_counts"].clip(lower=1))

keep_cells = (
    (spleen.obs["total_counts"] >= 250) &
    (spleen.obs["n_genes_by_counts"] >= 200) &
    (spleen.obs["pct_counts_mt"] <= 10) &
    spleen.obs["log10_genes_per_umi"] >=0.8)
    )
)

spleen = spleen[keep_cells].copy()

# -------------------------------------------------------- #
# doublet detection and removal
vc = spleen.obs["donor_id"].value_counts() # TO DO: double check name of donor id variable in metadata
keep = vc[vc >= 100].index
spleen_big = spleen[spleen.obs["donor_id"].isin(keep)].copy()

# choose a safe PC count globally
max_pcs = int(min(spleen_big.n_obs, spleen_big.n_vars) - 1)
npcs = max(1, min(15, max_pcs))

import scanpy as sc
sc.pp.scrublet(
    spleen_big,
    batch_key="donor_id",
    n_prin_comps=npcs,
    threshold=0.25,
)

spleen.obs["doublet_score"] = np.nan
spleen.obs["predicted_doublet"] = np.nan
spleen.obs.loc[spleen_big.obs_names, "doublet_score"] = spleen_big.obs["doublet_score"].values
spleen.obs.loc[spleen_big.obs_names, "predicted_doublet"] = spleen_big.obs["predicted_doublet"].values

spleen.obs["predicted_doublet"] = spleen.obs["predicted_doublet"].astype(str)
spleen.obs["doublet_score"] = spleen.obs["doublet_score"].astype(float)

mask_keep = ~(spleen.obs["is_doublet"].fillna(False))
spleen = spleen[mask_keep].copy()

# ---------- MOUSE ONLY BELOW ---------- #
# ---------------------------------------------- #
# -------------------------------------------------------- #
# for mouse only -> convert genes w/ homologene using homologene R script.

map = pd.read_csv("map.csv").dropna(subset=["human"]).set_index("mouse")["human"]

actual_genenames = spleen.var_names.isin(map.index)
spleen_mapped = spleen[:, actual_genenames].copy()
spleen_mapped.var_names = map.loc[spleen_mapped.var_names].values
spleen_mapped = spleen_mapped[:, ~spleen_mapped.var_names.duplicated()].copy()
# TO DO: double check it looks okay
# spleen = spleen_mapped

# -------------------------------------------------------- #
# ---------------------------------------------- #
# ---------- MOUSE ONLY ABOVE ---------- #


# -------------------------------------------------------- #
# varibale feature analysis
sc.pp.highly_variable_genes(spleen, n_top_genes=2000, batch_key="sample")
sc.pl.highly_variable_genes(spleen)

# -------------------------------------------------------- #
# normalization, logarithimization, etc.
sc.pp.normalize_total(spleen, target_sum=1e4)
sc.pp.log1p(spleen)
sc.pp.highly_variable_genes(
    spleen,
    layer="counts",
    n_top_genes=2000,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    flavor="seurat_v3",
)

# -------------------------------------------------------- #
# scaling and packages
sc.tl.pca(spleen)
sc.pl.pca_variance_ratio(spleen, n_pcs=50, log=True)

# -------------------------------------------------------- #
# harmony implementation
#import scanpy as sc
import scanpy.external as sce
import harmonypy

sce.pp.harmony_integrate(
    spleen,
    key='biosample_id',          # batch column in spleen.obs
    basis='X_pca',               # input embedding
    adjusted_basis='X_pca_harmony',  # output will be stored here
    max_iter_harmony=20,         # iterations
    theta=2.0                    # batch diversity penalty
)

sc.pp.neighbors(spleen, use_rep='X_pca_harmony', n_pcs=30)
sc.tl.umap(spleen)
sc.tl.leiden(spleen, key_added="leiden_harmony", resolution=1.0)
sc.pl.umap(spleen, color=["biosample_id", "leiden_harmony"])
#sc.tl.louvain(spleen, key_added="louvain", flavor="igraph")
#sc.pl.umap(spleen, color=["biosample_id", "louvain"])


# -------------------------------------------------------- #
# cell annotation w/ cell typist



# -------------------------------------------------------- #
# making files
