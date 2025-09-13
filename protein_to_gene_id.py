import pandas as pd, mygene

def stripver(x): return str(x).split(".", 1)[0]

def build_prot2gene(id_list):
    mg = mygene.MyGeneInfo()
    ids = [stripver(x) for x in id_list]
    ids = list(dict.fromkeys(ids))
    human = [x for x in ids if x.startswith("ENSP")]
    mouse = [x for x in ids if x.startswith("ENSMUSP")]
    m = {}

    def q(batch, species):
        if not batch: return
        res = mg.querymany(batch, scopes="ensembl.protein",
                           fields="ensembl.gene", species=species,
                           as_dataframe=False, verbose=False)
        for r in res:
            if r.get("notfound"): continue
            ens = r.get("ensembl")
            gene = None
            if isinstance(ens, dict): gene = ens.get("gene")
            elif isinstance(ens, list):
                for it in ens:
                    if isinstance(it, dict) and "gene" in it:
                        gene = it["gene"]; break
            if r.get("query") and gene:
                m[stripver(r["query"])] = stripver(gene)

    q(human, "human"); q(mouse, "mouse")
    return m

df = pd.read_csv("/Users/dr/cr/pr/class_II_project/in/maps/hsmm/mm_to_hs.txt", sep="\t", header=None)

all_prots = pd.unique(pd.concat([df[0], df[1]]).astype(str).map(stripver))

prot2gene = build_prot2gene(all_prots)

df[0] = df[0].astype(str).map(lambda x: prot2gene.get(stripver(x), x))
df[1] = df[1].astype(str).map(lambda x: prot2gene.get(stripver(x), x))

df.to_csv("/Users/dr/cr/pr/class_II_project/in/maps/out_mm_to_hs.txt", sep="\t", header=False, index=False)
