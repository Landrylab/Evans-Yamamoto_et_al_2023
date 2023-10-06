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
    pandas version >=1.3.4
    ```

    In the terminal, go to the location of the downloaded folder, and install the dependencies above by executing the following command.<br>
    ```sh
    pip install .
    ```

- R <br>
    ```sh
    sessioninfo   version >=1.2.2	
    ggplot2       version >=3.4.2
    reshape2      version >=1.4.4
    GGally        version >=2.1.2	
    ggridges      version >=0.5.4
    plyr          version >=1.8.8	
    dplyr         version >=1.1.2
    tidyr         version >=1.3.0
    tidyverse     version >=2.0.0	
    Cairo         version >=1.6.0	
    matrixStats   version >=1.0.0
    forcats       version >=1.0.0
    hardhat       version >=1.3.0	
    gridExtra     version >=2.3	
    ggExtra       version >=0.10.0
    egg           version >=0.4.5
    devtools      version >=2.4.5	
    ggtree        version >=3.6.2	
    castor        version >=1.7.10	
    treeio        version >=1.22.0	
    TreeTools     version >=1.9.2
    stringr       version >=1.5.0
    cowplot       version >=1.1.1	
    ggpubr        version >=0.6.0	
    gggenes       version >=0.5.0	
    ```
    
    To install these packages, execute the following script in the terminal.
    ```sh
    Rscript install_dependencies.r
    ```

### Other programs
- [Jupyterlab](https://jupyter.org/install) <br>
In this repository, most scripts are in jupyter notebook format. Installing jupyterlab would benefit to execute the scripts.
Install Jupyterlab by pasting the following in the terminal and press return.
    ```sh
    pip install jupyterlab
    ```

- Commandline BLAST+ <br>
Follow the [instruction manual](https://www.ncbi.nlm.nih.gov/books/NBK569861/) for installation.
    
- Commandline MAFFT<br>
Visit the [MAFFT website](https://mafft.cbrc.jp/alignment/software/) for installation.

- raxml-ng<br>
Visit the [raxml-ng github page](https://github.com/amkozlov/raxml-ng) for installation.

- PAML<br>
Visit the [PAML github page](https://github.com/abacus-gene/paml) for installation.
I followed the tutorial from [this tutorial paper](https://academic.oup.com/mbe/article/40/4/msad041/7140562?searchresult=1&login=false) and it's [github resource](https://github.com/abacus-gene/paml-tutorial).

- pyphe<br>
Visit the [pyphe page](https://pypi.org/project/pyphe/) for installation. 


## List of content and description

This repository contains the following folders. The folders are numbered in sequencial order for execution.

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
I assigned each ortholog a unique ID (present in the column GeneID_codeml), since codeml requires identifiers which are short. Using this file, I created inputs for downstrream analysis which are present in the folder 05_gene_tree_construction.

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
Scripts to reproduce Figure 1D-G of the paper.

**00_data_preparation**<br>
Data presented in the **raw_file** folder is proccessed using the script **alignment2nogap.ipynb** in order to create fasta files for codeml analysis. Some manual modifications (inserting the header for file format etc) was performed to ensure proper execution of codeml.

**01_codeml**<br>
In this folder, the inputs, control files (*.ctl), log files, and outputs from codeml are shown.

**02_evolution_rate_analysis**<br>
In this folder, intermediate files for generating figures based on codeml output is presented, as well as scripts and visualized output.

* domain_dNdS_heatmap.ipynb<br>
    Script to vizualize the domain based dN/dS values as heatmap.
* evolutionary_rate_analysis_R.ipynb<br>
    Script to analyze branch lengths and assymtry from codeml output (Figure 1E-G).
* Results<br>
    Folder containing all plots


### 04_Combinatorial_functional_complementation_screening
Scripts and output related to combinatorial complementation screening.

* Input<br>
    * Sample information for analysis
    * Image data from S&P imager (Available upon request to the corresponding author)
    * Numeric values extracted from the Image data ([available here](https://drive.google.com/file/d/1RS7EYXSvUvU3izMDShCKwzdJK3wu-c14/view?usp=drive_link))
* Scripts
    * 01_QuantifyAreaFromPlatePicture.ipynb<br>
      Script to extract colony area from each image.
      
    * 02_AUC_computation.ipynb<br>
      Script to compute Area Under the Curve from colony area information.
      
    * 03_parse_auc_data_2_scores_20230828.ipynb<br>
      Script to compute complementation scores, using AUC values in selectio nand non-selection conditions.
      
    * 04_plot_heatmap.ipynb<br>
      Script to plot heatmap from the complementation scores.

* Output<br>
Files generated from the scripts. Plots were used to prepare Figure 2C and Figure 2D of the paper.

### 05_DHFR-PCA_assay
Scripts and output related to the DHFR-PCA screening.

* Input<br>
    * Sample information for analysis
    * Image data from S&P imager (Available upon request to the corresponding author)
    * Numeric values extracted from the Image data (2022-12-09_MTX_Sel2_AUC_data_Cterm.csv)

* Scripts<br>
    * 01_robotpics_analysis.ipynb<br>
      Script to extract colony area from each image.
      
    * 02_AUC_computation.ipynb<br>
      Script to compute PPI scores, using AUC values.
    
    * 03_parse_screening_data.ipynb<br>
      Script to parse screening information and PPI data.
      
    * 04_Analysis.ipynb<br>
      Script to analyze PPI data and output stats.
          
* Output<br>
  Plots and intermediate files generated from the scripts. Plots were used to make Figure 3C and 3D of the paper.

### 06_GO_enrichment_analysis_of_PPI_partners
* Input<br>
    PPI score data (HRR25_orthologs_PPI_screening_parsed_2023-02-17DEY.csv) from 05_DHFR-PCA_assay.

* Scripts<br>
    * 01_data_proccessing.ipynb<br>
        Script to proccess PPI data and meta data for GO enrichment analysis.
    * 02_GO_Analysis.ipynb<br>
        Script to perform GO enrichment analysis on PPI partners.
          
* Output<br>
  Plots and files generated from the scripts. The folder GO_results contians csv files for GO enrichment analysis results for each ortholog's PPI partner, which is combined to one file as seen in GO_aggregated_results.csv. Figures were used to make Figure 3B of the paper.

### 07_SH3_domain_motif_analysis
* Input<br>
    * PPI score data (HRR25_orthologs_PPI_screening_parsed_2023-02-17DEY.csv) from 05_DHFR-PCA_assay.
    * pwm_dir (folder containing SH3 posision weight matrix from [this paper](https://pubmed.ncbi.nlm.nih.gov/26861823/)
    * Protein fasta files of HRR25 orthologs and the yeast proteome for motif search.
    * ID conversion file for SH3 proteins (yeast_sh3_accession_to_GN.txt).

* Scripts<br>
    * 01_motif_search.ipynb<br>
        Script to evaluate SH3 binding motifs in HRR25 orthologs.
    * 02_plot_PPIandSH3Motif.ipynb<br>
        Script to visualize the results.
          
* Output<br>
  Plots and files generated from the scripts. The folder contians a csv file (SH3_PWM_scan_HRR25Orthologs_MSS.csv) with all values from the PWM matches. Plots are as shown in Figure 3D of the paper.

