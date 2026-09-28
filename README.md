# Barrnapster

**Barrnapster** is an end-to-end pipeline designed to streamline the extraction, classification, and visualization of 18S rRNA eukaryotic sequences from long-read metagenomic assemblies. 

Originally built to parse through noisy environmental sequences to find protist hosts, Barrnapster automates the tedious transitions between standard bioinformatics tools. It predicts rRNA genes, extracts the 18S sequences, classifies them against the PR2 database, filters out bacterial and mitochondrial misassignments, generates publication-ready diversity visualizations, and retrieves the full-length original contigs for your target eukaryotes.

##  How It Works

Barrnapster executes a 6-step pipeline in a single command:

1. **Prediction:** Runs `barrnap` to identify eukaryotic rRNA regions within the input assembly.
2. **Filtering:** Parses the resulting GFF file to isolate exclusively 18S rRNA predictions.
3. **Extraction:** Uses `bedtools getfasta` to slice the 18S nucleotide sequences from the original assembly.
4. **Classification:** Uses `vsearch` with the SINTAX algorithm to classify the 18S sequences against the PR2 (Protist Ribosomal Reference) database.
5. **Visualization:** Cleans the taxonomic output, removes non-target hits (Bacteria, Mitochondria), and leverages `matplotlib`/`pandas` to generate horizontal bar charts, composition donuts, and summary statistics across multiple taxonomic ranks.
6. **Contig Retrieval:** Maps the filtered eukaryotic 18S sequences back to the raw assembly and extracts the full parent contigs into a clean FASTA file for downstream genomic analysis.

---

##  Installation & Dependencies

Barrnapster relies on a mix of command-line bioinformatics utilities and Python data science libraries. The easiest and most reproducible way to install these dependencies is via Conda/Mamba.

### 1. Create a Conda Environment

Create a new environment with all required tools (specifying `bioconda` and `conda-forge` channels):

```bash
conda create -n barrnapster -c conda-forge -c bioconda \
    python=3.11 \
    barrnap \
    bedtools \
    vsearch \
    biopython \
    pandas \
    numpy \
    matplotlib \
    -y
```
After creating the environment, activate it.
```bash
conda activate barrnapster
```
### 2. Setting up the PR2 database
You will need to install the PR2 database from their website [here](https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_UTAX.fasta.gz).
After installing, unzip the database.

### 3. Running the pipeline
The input for the tool is a metagenome assembly (a FASTA file with contigs). 
```bash
python barrnapster.py -a [assembly.fasta] -d [path/to/PR2_database] -o [output_folder]
```



