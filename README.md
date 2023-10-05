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
    Protein fasta file containing the 206 identified orthologs in the first alignment.
  
* 2023-03-02_NCBI_BLASTp_SGD_hits_parsed.xlsx
    Excel file containing BLASTp results of Saccharomycetaceae_BLASTp_hits.fasta against the [SGD database (_S. cerevisiae_ proteins)](https://www.yeastgenome.org/blast-sgd).
  
* Saccharomycetaceae_Hrr25_summary.csv
    csv file with summary of extracted Hrr25p sequences, removing all false positives. It also includes annotated orthologs in the YGOB database. the 

**02_Phylogenetic_tree**<br>
* 1672taxa_290genes_bb_1.treefile<br>
    Phylogenetic tree file from [Li et al. (2021) _Current Biology_](https://www.cell.com/current-biology/fulltext/S0960-9822(21)00139-1)
  
* tree_Li_etal_2021.ipynb<br>
    R script in jupyter notebook to load and trim the phylogenetic tree, based on a set of species presnet in ../01_RefSeq_Protein_retrival/Saccharomycetaceae_Hrr25_summary.
  
* selected_species_tree.txt<br>
    Trimmed tree output from the script.
  
* sel_species.csv<br>
    Output from the script, with list of spcecies present in the trimmed tree.
  
* SelectedSpeciesTree_plot.pdf<br>
    Visualized tree output from tree_Li_etal_2021.ipynb.

**03_Extended_homolog_search**
* search_homolog.ipynb<br>
    Script to perform BLAST alignments agaisnt all genoe sequences using ./blast/db/HRR25_nuc_nonAligned.fasta as query.
  
* blast<br>
    Folder containing database and outputs from BLAST alignments.
  
* genomes.zip<br>
    Compressed folder with genoome sequences which orthologs are going to be retrieved from. Since this file is too large to upload to github, it is available [here](https://drive.google.com/file/d/1BQuC4V6KXw-DSiV1K4kb3OfkjmJaOb3r/view?usp=drive_link).
  
* blast_hits.csv<br>
    A file containing all BLAST hits, present in ./blast/out.
  
* unique_regions_to_extract.csv<br>
   Unique gene regions parsed from blast_hits.csv
  
* HRR25_homologs_nt_extracted.fna<br>
    Fasta file containing homologs identified from genomic seuquences.
  
* HRR25_merged.fna<br>
    The result from this folder (HRR25_homologs_nt_extracted.fna) was merged with the input for homology search (./blast/db/HRR25_nuc_nonAligned.fasta) to be used for downstream analysis. 

**04_Cleanup_homolog**
* alignment4ORFdetection.ipynb<br>
    Python script to perform MAFFT-linsi and identify ORF regions for HRR25_merged.fna.
  
* HRR25_mafft_linsi.txt<br>
  Output from MAFFT-linsi.
  
* HRR25_homologs_aa_trimmed.fna<br>
    Output from alignment4ORFdetection.ipynb, contiaining protein sequences in fasta format.
  
* HRR25_trimed_aa_info.csv<br>
    Output from alignment4ORFdetection.ipynb, contiaining protein sequences in csv format.
  
* HRR25_homologs_nt_trimmed.fna<br>
    Output from alignment4ORFdetection.ipynb, contiaining nucleotide sequences in fasta format.
  
* HRR25_trimed_nt_info.csv<br>
    Output from alignment4ORFdetection.ipynb, contiaining nucleotide sequences in csv format.

* **TableS1_ListofGenes.xlsx**<br>
The output from 04_Cleanup_homolog was used to create a list of orthologs presented in Supplementary Table 1 (TableS1_ListofGenes.xlsx) of the manuscript.
I assigned each ortholog a unique ID (present in the column GeneID_codeml), since codeml requires identifiers which are short. Using this file, I created which are present in the folder 05_gene_tree_construction.

**05_gene_tree_construction**
* HRR25_geneanalysis_aa.fna and HRR25_geneanalysis_nt.fna<br>
    Fasta files containing the ortholog sequences identified by unique IDs, created from TableS1_ListofGenes.xlsx.
  
* trim_protein.ipynb<br>
    Python notebook to create inputs for TranslatorX, a program to perform alignment based on codons.
  
* HRR25_geneanalysis_aa_trimmed.fna          <br>
    Output from trim_protein.ipynb, where protein sequence is properly annotated (excluding regions after stop codons etc).
  
* HRR25_geneanalysis_nt_translatorXinput.fna <br>
    Output from trim_protein.ipynb, with nucleotide sequences corresponding to HRR25_geneanalysis_aa_trimmed.fna. I use file this for input in TranslatorX.
  
* translatorX_perl<br>
    A folder containing scripts from [TranslatorX](https://translatorx.org/downloads.html)
  
* translatorX_res<br>
    A folder containing results from TranslatorX, using HRR25_geneanalysis_nt_translatorXinput.fna as input.
  
* raxml_res<br>
  A folder containing scripts and results from [raxml-ng](https://github.com/amkozlov/raxml-ng). I created the input file which only contains orthologs from post-WGD species which maintained two orthologs (HRR25_mafft_translatorx.nt_ali_PostWGD_selected.fasta) from the output of TranslatorX (HRR25_mafft_translatorx.nt_ali.fasta). The resulting tree was used manually create the recomciliated tree y replacing the post-WGD species with maintained duplicates with the tree presented in HRR25_mafft_translatorx.nt_ali_PostWGD_selected.fasta.raxml.bestTree. The resulting tree can be found in HRR25_genetree_postWGDGeneTreeIntegrated_ID_M0.txt.

### 02_Multiple_Sequence_Alignment_analysis
Scripts to reproduce Figure 1C of the paper.

* input<br>
    Input for this analysis is the codon based alignment of orthologs, identical to ../05_gene_tree_construction/translatorX_res/HRR25_mafft_translatorx.aa_ali.fasta.
   
* meta_data<br>
    Folder with meta data, includig domain annotations and position information to aid interpretation of the plots.
    
* output<br>
    Folder with outputs, including intermediate files with similarity scores by position.
  
* msa_analysis.ipynb<br>
    A R script in jupyter-notebook, which was used to calculate the similarity score for each residue in orthologs.
  
* plot_similarity.ipynb<br>
    A R script in jupyter-notebook, which was used visualize the data as presented in Figure 1C.

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
