<!--[![DOI](https://zenodo.org/badge/829822589.svg)](https://zenodo.org/doi/10.5281/zenodo.13347495)-->
# Overview of the scMnT method
![scMnT preview](https://github.com/18parkky/scMnT/blob/main/Overview.jpg)

# Table of contents
1. [Installation](#installation)<br/>
2. [Commands & Parameters](#commands--parameters)<br/>
3. [Output files](#output-files)<br/>
4. [Using scMnT](#using-scmnt)<br/> 
5. [Tutorial - Identifying MSI status in single-cell resolution](#tutorial---identifying-msi-status-in-single-cell-resolution-using-cancer-cell-data-from-a-study-by-kinker-et-al)<br/>
6. [Citation](#citation)<br/>

# scMnT
scMNT is a specialized extension of NanoMnT, developed for a distinct purpose: quantifying microsatellite instability (MSI) intensity at single-cell resolution. It provides a suite of Python scripts that analyze reads mapped to mononucleotide-repeat microsatellite (MS) regions.

## Installation
To run scMnT, minimap2 and the following Python packages are required:
  
  - Matplotlib >= 3.7.1
  - Numpy >= 1.20.3
  - Pysam >= 0.20.0
  - Pandas >= 2.0.0
  - Scipy >= 1.7.1
  - Seaborn >= 0.13.0

After installing these packages, use the following commands to install scMnT.
*You may consider creating a new Conda environment for scMnT.
```
conda create -n scmnt python=3.12.5 # mamba create -n scmnt python=3.12.5
git clone https://github.com/18parkky/scMnT.git
cd scMnT
pip install -e .
```

scMnT requires a tab-delimited file to specify MS loci to analyze that includes the following columns: `sequence`, `standard`, `motif`, `type`, `repeat`, `start`, `end`, `length`.
Files of this format can be easily be created using [Krait](https://github.com/lmdu/krait) (v1.3.3≥, not to be confused with Krait2) `Search for SSR` function, which needs to be installed separately.
MS loci files for mononucleotide-repeat regions (≥10 bp) in the human genome (GRCh38 and T2T-CHM13) and the mouse genome (GRCm38) are available in scMnT GitHub repository under:`~/scMnT/ref/`. 
- `~/scMnT/ref/total` : Contains all MS regions. We recommend using these files when the sequencing data is not high-coverage.
- `~/scMnT/ref/filtered` : Contains only MS regions flanked by sequences with high sequencing complexity. Although these files are smaller (representing only about a fraction of the size of those in ~/scMnT/ref/total) we recommend using them whenever possible, as MS allele estimation tends to be more accurate in these regions.

## Commands & Parameters
scMnT has three steps, all of which can be run with a single command: `scMnT`.
Each step can also be called separately: (1) `getAlleleTable`, (2) `scMnT-score`, (3) `scMnT-find`.

### 1) `scMnT`
The `scMnT` command runs the entire scMnT workflow, from (1) collecting reads that map to MS loci (getAlleleTable), (2) labeling MS information for each cell onto the user-provided Scanpy object (scMnT-score), then (3) identifying MSI cell types by comparing MSI score distributions (scMnT-find). `scMnT` requires the following files as inputs. **In most cases, running this single command is all that’s needed.**
- Scanpy object where cell type (or any other categorical cell annotation, e.g., Leiden clusters) information is available at .obs['CellType'] and cell barcodes (e.g., TAGTTGGCAGGTCCAC-1) must be available at .obs.index. Please ensure that each cell type has enough number of cells (>100).
- If only using a single sample as input: a BAM where each read is tagged by its CB (CB tag) and UMI (UB tag)
- If using multiple samples as input: a BAM where each read is tagged by its CB (CB tag) and UMI (UB tag)
- a list of MS loci (which can be generated using Krait)
- the reference genome in FASTA format

```
usage: scMnT [-h] -a PATH_SCANPY_OBJ -b PATH_BAM -s PATH_STR_TSV -r PATH_REFERENCE_GENOME -n FILENAME [-f FLANKING] [-m MINIMUM_LOCI] [-min_l MINIMUM_MS_LENGTH] [-max_l MAXIMUM_MS_LENGTH] [-t THREADS] [-out DIR_OUT]
```

**Required parameters**:
- `-a`: Path to Scanpy object (.h5ad) where cell type information is available at .obs['CellType'] and cell barcodes (e.g., TAGTTGGCAGGTCCAC-1) must be available at .obs.index. Tumor cells must be labeled 'Tumor' (e.g., Tumor 1, Tumor_2, Tumor III)
- `-s`: PATH to STR list file generated using Krait
- `-r`: PATH to reference genome (must be same as the one used for generating BAM)
- `-n`: Name of the output files

**Optional parameters**:
- `-b`: PATH to input BAM file (must be sorted and indexed, required if sample sheet is not provided)
- `-ss`: Tab-delimited file with two columns: (1) sample name and (2) path to its BAM file. Don't include column name. Required when running scMnT on multiple samples
- `-f`: Length of flanking sequences to collect for each read (default: 6)
- `-m`: Minimum number of MS loci required per cell to be included in the analysis (default: 10). Cells with fewer loci will be excluded from analysis.
- `-min_l`: Minimum length of the MS loci to assess (default: 10)
- `-max_l`: Maximum length of the MS loci to assess (default: 24)
- `--resume`: Resume from the previous run (default: False)
- `-t`: Number of threads to use for multiprocessing (default: 4)
- `-out`: Directory to write output files (default: current directory)

<details>
<summary><h3>2) <code>getAlleleTable</code> Click to see description</h3></summary>
The `getAlleleTable` command requires (1) a BAM where each read is tagged by its CB (CB tag) and UMI (UB tag) and (2) a list of MS loci (generated using Krait), and produces a single TSV file, namely, the allele table. The allele table contains information about reads that aligned to any of the given MS loci. 

In a typical scenario, you will be using the BAM file generated by Cell Ranger (e.g., `possorted_genome_bam.bam`) and the Krait file under `~/scMnT/ref/` as the inputs.

```
usage: getAlleleTable [-h] -b PATH_BAM -s PATH_STR_TSV -r PATH_REFERENCE_GENOME -n FILENAME [-f FLANKING] [-t THREADS] [-out DIR_OUT]
```

**Required parameters**:
- `-b`: PATH to input BAM file (must be sorted and indexed)
- `-s`: PATH to STR list file generated using either Krait
- `-r`: PATH to reference genome (must be same as the one used for generating BAM)
- `-n`: Name of the output files

**Optional parameters**:
- `-f`: Length of flanking sequences to collect for each read (default: 6)
- `-t`: Number of threads to use for multiprocessing (default: 4)
- `-out`: Directory to write output files (default: current directory)
</details>

<details>
<summary><h3>3) <code>scMnT-score</code> Click to see description</h3></summary>
The `scMnT-score` command requires (1) an Scanpy object file (h5ad) where cell type information must be available at `adata.obs['CellType']` and (2) the allele table generated using `getAlleleTable` as inputs. `scMnT-score` then calculates the MSI score for each cell (saving it to `adata.obs['MSI_score']`) and writes the resulting Scanpy object file to disk.

```
usage: scMnT-score [-h] -a PATH_SCANPY_OBJ -at PATH_ALLELE_TABLE -n FILENAME [-min_l MINIMUM_MS_LENGTH] [-max_l MAXIMUM_MS_LENGTH] [-out DIR_OUT]
```
**Required parameters**:
- `-a`: Path to Scanpy object (.h5ad) where cell type information is available at .obs['CellType'] and cell barcodes (e.g., TAGTTGGCAGGTCCAC-1) must be available at .obs.index. Tumor cells must be labeled 'Tumor' (e.g., Tumor 1, Tumor_2, Tumor III)
- `-at`: PATH to the Allele Table (generated using getAlleleTable.py)
- `-n`: Name of the output files

**Optional parameters**:
- `-min_l`: Minimum length of the MS loci to assess (default: 10)
- `-max_l`: Maximum length of the MS loci to assess (default: 24)
- `-out`: Directory to write output file to (default: current directory)
</details>

<details>
<summary><h3>4) <code>scMnT-find</code> Click to see description</h3></summary>
The `scMnT-find` command finds cell types that are highly likely to be MSI by comparing MSI score distributions. As inputs, `scMnT-find` requires the Scanpy object generated using `scMnT-score`. `scMnT-find` will then produce a TSV file that tells the MSI status of each cell type. 

```
usage: scMnT-find [-h] -a PATH_SCANPY_OBJ -n FILENAME [-m MINIMUM_LOCI] [-out DIR_OUT]
```

**Required parameters**:
- Scanpy object 

- `-a`: Path to Scanpy object (.h5ad) where cell type information is available at .obs['CellType'] and cell barcodes (e.g., TAGTTGGCAGGTCCAC-1) must be available at .obs.index. Tumor cells must be labeled 'Tumor' (e.g., Tumor 1, Tumor_2, Tumor III)
- `-n`: a comma separated list of normal cell types. (e.g., monocyte,T,fibroblast)

**Optional parameters**:
- `-m`: Minimum number of MS loci required per cell to be included in the analysis (default: 10). Cells with fewer loci will be excluded from analysis.
- `-out`: Directory to write output files (default: current directory)
</details>

## Output files

### 1) `scMnT` → (1) allele table, (2) Scanpy object labeled with MSI score, and (3) MSI information of each cell type.
See below for detailed explanation of each output file.

### 2) `getAlleleTable` → allele table
As mentioned above, `getAlleleTable` produces the 'allele table', of which each row contains the following information:
- Read name (`read_name`)
- Locus (`locus`)
- Repeat unit (`repeat_unit`)
- STR allele sequence reported by the read (e.g., AAAAAAAA) (`allele`)
- STR allele (=number of repeats) of the reference genome (`reference_STR_allele`)
- Left flanking sequence of STR (`left_flanking_seq`)
- Right flanking sequence of STR (`right_flanking_seq`)
- Flag of the read (`flag`)
- Cell barcode (`CB`)
- Unique molecular identifier (`UMI`)
- Corrected allele (error-corrected allele) (`corrected_allele`)
- Editing distance between the uncorrected allele and the corrected allele (`editing distance`)
- STR allele reported by the read (`read_STR_allele`)

### 3) `scMnT-score` → Scanpy object file
Running `scMnT-score` creates a Scanpy object file with MSI score labeled under `adata.obs['MSI_score']`.

### 4) `scMnT-find` → TSV file
Running `scMnT-find` creates a TSV file with the following columns: 'CT', 'pval', 'delta', 'n_cells'. 
- CT = Cell type
- pval = the p-value of the null hypothesis — that the MSI score distribution of the given cell type originates from the same distribution as that of the normal cell types — as tested by the Kolmogorov–Smirnov test
- delta = Cliff's delta value
- Number of cells used for the Kolmogorov-Smirnov test

## Using scMnT
Generally, the scMnT workflow is as follows:
1. Cell type annotation (e.g., using Scanpy/Seurat)
2. `scMnT` command
3. Downstream analyses (e.g., UMAP inspection of MSI scores and subsequent subclustering).

## FAQ
1. What MSI score threshold should I use to confidently classify MSI cells? 
  - While we understand from the user perspective that a clear-cut threshold value would be nice, we recommend users to first inspect the MSI score distribution among the normal cells (which are surely MSS in the cancer context).
  - Depending on the sequencing configuration and library prep methods, the levels of noise can vary significantly across different data. Thus, we advise users to use the normal cells' MSI scores as reference to decide the threshold. This approach is demonstrated in our tutorial notebook (`~/scMnT/tutorial.ipynb`).

2. TBA

## Tutorial - Identifying MSI status in single-cell resolution using cancer cell data from a study by [Kinker et al.](https://www.nature.com/articles/s41588-020-00726-6)
We provide a tutorial notebook(`~/scMnT/tutorial.ipynb`) that well-represents the full capability of scMnT.
**We highly recommend the users to go through the tutorial (or even use them for their analyses).**

This tutorial uses the cancer cell line dataset from a study by Kinker et al., guiding the users on how to proceed after running `scMnT`, covering how to:
- identify MSI cells
- estimate MSI intensity
- visualize results using UMAP

You are free to reproduce the results using the same data, or you can simply replace the input files with your own to run the tutorial.

## Citation  
If you use scMnT in your research, please cite our [preprint](https://www.biorxiv.org/content/10.1101/2025.05.09.653227v1).
We will update the DOI when we publish our paper in a peer-reviewed journal.

Like NanoMnT, scMnT is under active development, so feel free to make suggestions!
