# implementing homologene to convert names from murine to human

# TO DO: check if vector has header or not.
genenames <- read.delim("myfile.tsv", header = TRUE, stringsAsFactors = FALSE)
#is.data.frame(genenames)
#genenames <- genenames[-1, ]
#is.vector(genenames)
#genenames <- genenames[-1]

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("homologene")
BiocManager::install("homologene.data")

library(homologene)
library(homologene.data)

# get homologs
new_genenames <- homologene(mouse_genes,
    inTax    = 10090,
    outTax   = 9606,
    db       = homologene.data::homologeneData)

# make map
map <- data.frame(mouse = genenames) |>
  merge(new_genenames[, c("input_gene", "homolog_gene")],
        by.x = "mouse", by.y = "input_gene", all.x = TRUE)
names(map) <- c("mouse", "human")

head(map)
