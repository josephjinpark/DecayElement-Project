# Survey of the human genome for potential miRNA decay elements (miRNA sponges)

This code from back in the day surveyed the human genome for potential microRNA sponges (decay elements).

It utilizes the miRBase database, https://www.mirbase.org, to identify segments in the genome that can potentially host multiple target sites for human microRNAs.
Pipeline flow:
- Segmentation of the human genome according to size
- Prediction of the secondary structure of the segment using ViennaRNA toolkit, http://rna.tbi.univie.ac.at
- Parsing and annotation of dot-bracket notations for miRNA target sites at open secondary structures (bulges)
- Visualization of the secondary structures and seed sites using svg.

Example:
![miR-1 1](https://github.com/josephjinpark/DecayElement-Project/assets/23091681/26aa2ee1-baa0-4e2d-aeff-3aeb21bd8cdc)

