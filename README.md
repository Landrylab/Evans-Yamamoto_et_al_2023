# Evans-Yamamoto_et_al_2023

This repository contains the computer codes and control files, defining the settings for the codeml program execution, used in the submitted manuscript [**Parallel nonfunctionalization of CK1δ/ε kinase ohnologs following a whole-genome duplication event**](https://doi.org/10.1101/2023.10.02.560513).


## Installation

Please make sure you have appropriate Python, pip, and R before starting.
```sh
Python version >= 3.5
pip    version >= 1.1.0
R      version >= 4.2.2
```

Download scripts by first clone this repository by execuiting the following command in the terminal.

```sh
git clone https://github.com/LandryLab/EVANS-Yamamoto_et_al_2023.git
```



### Dependencies 

- Python <br>
    ```
    numpy  version >=1.19 
    ```

    In the terminal, go to the location of the downloaded folder, and install the dependencies above by executing the following command.<br>
    ```sh
    pip install .
    ```

- R <br>
    ```sh
    XXX version >= 1.0.0
    ```

    To install these packages, execute the following script in the terminal.
    ```sh
    Rscript install_dependencies.r
    ```

### Other programs
- [Jupyterlab](https://jupyter.org/install) <br>
While I provide html codes of jupyter notebooks, installing jupyterlab would benefit to execute the scripts.
Install Jupyterlab by pasting the following in the terminal and press return.
    ```sh
    pip install jupyterlab
    ```

- Commandline BLAST+ <br>
Follow the [instruction manual](https://www.ncbi.nlm.nih.gov/books/NBK569861/) for installation.
Set the PATH of the binary file.

    Execute the following to see if installation is complete.
    ```sh
    blastn -help
    ```
- Commandline MAFFT<br>
- PAML<br>
- 



## Content

This repository contains the following folders.

### 00_Preliminary_analysis
Scripts regarding the preliminary analysis on YGOB data base, idintifying essential genes in _S.cerevisiae_, which are maintained as duplicates in other species.
* YGOB_wgd_essentiality_stats.csv <br>
    Input data, created from Gene order & annotation from the [YGOB database](http://ygob.ucd.ie/Pillars.tab) and gene essentiallity information from the [SGD database](https://www.yeastgenome.org/observable/APO:0000112).
  
- YGOB_stats_R.ipynb/.html <br>
  Script to filter and count species maintaining ohnologs for each gene.
  
- YGOB_ScerEssential_ScerCount1_ZscorePostWGDOver2.csv <br>
  Stats outputed from YGOB_stats_R.



### 01_Ortholog_sequence_retrieval
Scripts regarding **Ortholog sequence retrieval** section in the manuscript.
It contains the following folders and files;

**01_RefSeq_Protein_retrival**<br>
* Saccharomycetaceae_species.csv<br>
    List of _Saccharomycetaceae_ species in the NCBI database.
* download_ncbi_genomes.sh<br>
    Script to download files from NCBI.
* NCBI_download_wrapper.ipynb<br>
    Python notebook to download the genomes and protein files for the species listed in Saccharomycetaceae_species.csv.
* 2023-03-01_NCBI_download_summary.csv<br>
    Intermediate output from NCBI_download_wrapper.ipynb.
* hrr25.faa<br>
    Fasta file containing the _S. cerevisiae_ Hrr25p sequence.
* blast.ipynb  <br>                        
    Python notebook to create BLASTp databases from the downloaded protein files (under ./blastp/db), and perform BLASTp using the _S. cerevisiae_ Hrr25p (output under ./blastp/out).
* blastp<br>
    Folder containing intermediate files for protein blast.
* 2023-03-02_BLASTp_parsed.csv<br>
    Parsed data from the protein blast.
* Saccharomycetaceae_BLASTp_hits.fasta<br>
    Protein fasta file containing the identified orthologs.

**02_Phylogenetic_tree**<br>
* 1672taxa_290genes_bb_1.treefile<br>
    Phylogenetic tree file from [Li et al. (2021) _Current Biology_](https://www.cell.com/current-biology/fulltext/S0960-9822(21)00139-1)
* tree_Li_etal_2021.ipynb<br>
    Script to load and trim the phylogenetic tree, based on a set of species which protein files were downloaded.
* selected_species_tree.txt<br>
    Trimmed tree output from the script.
* sel_species.csv<br>
    Output from the script, with list of spcecies present in the trimmed tree.
* SelectedSpeciesTree_plot.pdf<br>
    Visualized tree output from the script.





### 02_Multiple_Sequence_Alignment_analysis
hoge

### 03_dNdS_analysis
hoge

### 04_Combinatorial_functional_complementation_screening
hoge

### 05_DHFR-PCA_assay
hoge

### 06_GO_enrichment_analysis_of_PPI_partners
hoge

### 07_SH3_domain_motif_analysis
hoge
