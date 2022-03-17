# Sexually divergent development of depression-related brain networks during healthy human adolescence

This repository contains the code for the main analyses of the manuscript "Sexually divergent development of depression-related brain networks during healthy human adolescence" by Lena Dorfschmidt, Richard A.I. Bethlehem, Jakob Seidlitz, František Váša, Simon R. White, Rafael Romero-Garcia, Manfred G. Kitzbichler, Athina Aruldass, Ian M. Goodyer, Peter Fonagy, Peter B. Jones, Raymond J. Dolan, the NSPN Consortium, Petra E. Vértes, and Edward T. Bullmore* (2020).

For details behind these analyses refer to the manuscript: https://doi.org/10.1101/2020.07.06.184473

# Data
The data used here was initially published by [Váša et al. (2020)](https://doi.org/10.6084/m9.figshare.11551602). Please cite them, when using these data. 

All data required to run these analyses can be found at: 

# How to Run
To run this code, download the required data. Most of the scripts read in ouputs from other scripts, so the order in which you run them is essential.

1. Run `maturational_index/maturational_index.R` to generate the result on *sex differences in adolescent brain maturational*. 
2. Run `mdd_decoding/mdd.comparison.R` to run the *co-location analyses with depression*.
4. To run the *enrichment analyses for cell-types, chromosomes and disorders*:

  4.1. Run `PLS1_stats.m` to run the PLS analysis using Allen Human Brain Atlas data.
  
  4.2. Run `PLS2_bootstrap.m` to bootstrap the obtained PLS weights. 
  
  4.3. Now you can run `gene_decoding/gene_enrichment.R` to run the enrichment analyses for cell-types, chromosomes and disorders.

6. Finally, running `neurosynth/wordcloud.R` generates the wordcloud maps in manuscript. 

