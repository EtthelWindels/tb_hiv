# Estimating the effects of HIV co-infection on *Mtb* transmission using phylodynamics
This repository contains the code associated with the phylodynamic analyses performed in Windels et al. HIV co-infection is associated with reduced *Mycobacterium tuberculosis* transmissibility in sub-Saharan Africa. *PLOS Pathogens* 20, e1011675 (2024). [https://doi.org/10.1371/journal.ppat.1011675](https://doi.org/10.1371/journal.ppat.1011675)

These analyses aimed at estimating the effect of HIV co-infection on the effective reproductive number for TB and the risk of developing active TB, by fitting a customized birth-death model to genomic sequencing data collected from TB patients in Malawi, South Africa, Tanzania, and Uganda.

The raw sequencing data used in this study are available under the following project accession numbers: 
- Malawi: [PRJEB2358](https://www.ebi.ac.uk/ena/browser/view/PRJEB2358) and [PRJEB2794](https://www.ebi.ac.uk/ena/browser/view/PRJEB2794)
- South Africa: [PRJEB45389](https://www.ebi.ac.uk/ena/browser/view/PRJEB45389) and [PRJNA670836](https://www.ebi.ac.uk/ena/browser/view/PRJNA670836)
- Tanzania: [PRJEB49562](https://www.ebi.ac.uk/ena/browser/view/PRJEB49562)
- Uganda: [PRJEB11460](https://www.ebi.ac.uk/ena/browser/view/PRJEB11460), [PRJNA354716](https://www.ebi.ac.uk/ena/browser/view/PRJNA354716), and [PRJEB64921](https://www.ebi.ac.uk/ena/browser/view/PRJEB64921)

The repository is organized as follows:
- The folder `analyses/` contains the XML files used for the birth-death analyses in BEAST2.
- The folder `data/` contains all metadata.
- The folder `scripts/` contains R scripts used for pre-processing of the alignments and post-processing of the BEAST2 output.