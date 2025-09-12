# run samap
# TO DO: ensure that gene names in murine dataset are MURINE symbol!
# # TO DO: transfer cell annotations to file with murine symbols (using cell barcode)

# -------------------------------------------------------- #
# load datasets
fn1 = '/dr/cr/pr/class_II_project/in/spleen1.h5ad'
fn2 = '/dr/cr/pr/class_II_project/in/spleen2.h5ad'

filenames = {'mm':fn1,'hs':fn2}

# -------------------------------------------------------- #
from samalg import SAM
# make sam object
sam1=SAM()
sam1.load_data(fn1)

sam2=SAM()
sam2.load_data(fn2)

sams = {'pl':sam1,'sc':sam2}

# -------------------------------------------------------- #
# run samap
from samap.mapping import SAMAP
# TO DO: move maps into ~/in/maps
sm = SAMAP(
        sams,
        f_maps = '/dr/cr/pr/class_II_project/in/maps/',
    )

sm.run(pairwise=True, neigh_from_keys=["cell_type"], umap = True)
