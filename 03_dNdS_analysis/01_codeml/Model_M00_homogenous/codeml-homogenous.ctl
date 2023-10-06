     seqfile = /home/daney/projects/codeml/hrr25/HRR25_mafft_translatorx_nt_ali_nogaps.phy            * Path to the alignment file
     treefile = /home/daney/projects/codeml/hrr25/HRR25_genetree_postWGDGeneTreeIntegrated_ID_M0.txt * Path to the tree file
      outfile = /home/daney/projects/codeml/hrr25/Model_M00_homogenous/HRR25_codemloutput_M0.txt            * Path to the output file
 
        noisy = 3              * How much rubbish on the screen
      verbose = 1              * More or less detailed report
   seqtype = 1              * Data type
        ndata = 1              * Number of data sets or loci
        icode = 0              * Genetic code 
    cleandata = 0              * Remove sites with ambiguity data?
		
        model = 0          * Models for ω varying across lineages
	NSsites = 0        * Models for ω varying across sites
    CodonFreq = 7          * Codon frequencies
	estFreq = 0        * Use observed freqs or estimate freqs by ML
        clock = 0          * Clock model
        omega = 0.5        * Initial or fixed omega
    fix_omega = 0          * Estimate or fix omega
