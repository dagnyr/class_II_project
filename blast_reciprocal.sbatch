#!/bin/bash
#SBATCH --job-name=blast_hs_mm
#SBATCH --output=blast_%j.out
#SBATCH --error=blast_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal

# load ncbi blast
module load biology ncbi-blast+/2.16.0

# move to home directory
cd $HOME

# load fasta (protein)
HS=Homo_sapiens.GRCh38.pep.all.fa
MM=Mus_musculus.GRCm39.pep.all.fa

# make BLAST databases
makeblastdb -in $HS -dbtype prot -out hs_db
makeblastdb -in $MM -dbtype prot -out mm_db

# run reciprocal blast
blastp -query $HS -db mm_db -out hs_to_mm.txt -evalue 1e-6 -outfmt 6 -num_threads 8
blastp -query $MM -db hs_db -out mm_to_hs.txt -evalue 1e-6 -outfmt 6 -num_threads 8
