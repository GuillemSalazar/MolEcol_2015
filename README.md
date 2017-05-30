# MolEcol_2015

This GitHub Repository contains code and data included in:

--------
Salazar G., Cornejo-Castillo F.M., Borrull E., Díez-Vives C., Lara E., Vaqué D., Arrieta J.M., Duarte C.M., Gasol J.M. and Acinas S.G. **Particle-association lifestyle is a phylogenetically conserved trait in bathypelagic prokaryotes**. (2015) Molecular Ecology. doi:[10.1111/mec.13419](http://onlinelibrary.wiley.com/doi/10.1111/mec.13419/abstract)

--------

Data is also accessible at [*Dryad*](https://datadryad.org/handle/10255/3/workflow?workflowID=59300&stepID=reviewStep&actionID=reviewAction)

## Files

 - **OTUtable_Salazar_etal_2015_Molecol.txt**: Abundance table for bathypelagic prokaryotes from Malaspina 2010 Expedition.  
 Abundance table containing the number of reads for the Operational Taxonomic Units (OTUs) of particle-attached (PA) and free-living (FL) prokaryotes detected in 30 globally distributed stations from the Malaspina 2010 Expedition. Corresponds to 16SrDNA amplicons Illumina-based sequencing. Table is rarefied to 10,617 reads/sample. Columns are the samples; Rows are the OTUs. Taxonomy is attached after the last sample.  

 - **OTUtable_Salazar_etal_2015_Molecol_norarefac.txt**: Abundance table for bathypelagic prokaryotes from Malaspina 2010 Expedition (without rarefaction).  
 Abundance table containing the raw number of reads, i.e. before the rarefaction. 
  
 - **Metadata_Salazar_etal_2015_Molecol.txt**: Auxiliary data for the samples.
 Auxiliary data containing information on the station, size-fraction, date of sampling, ocean, depth and coordinates.
 
 - **FinalTree_Salazar_etal_2015_Molecol.nwk**: Phylogeny of the OTUs representative sequences
 - **Alignment_Salazar_etal_2015_Molecol.phylip**: Alignment used for the phylogeny contruction (containing the OTU representative sequences and the closest SILVA v.115 sequences)
 - **RScript_Salazar_etal_2015_MolEcol.R**: R script reproducing analyses in the publication.
 - **OTUs.fasta**: Fasta file of OTUs' representative sequences (matching the ones in *OTUtable_Salazar_etal_2015_Molecol.txt*)

## Authors

Guillem Salazar

## Citation

[![DOI](https://zenodo.org/badge/18788/GuillemSalazar/MolEcol_2015.svg)](https://zenodo.org/badge/latestdoi/18788/GuillemSalazar/MolEcol_2015)

## Software Versions

*R*: version 3.2.1  
*RStudio*: Version 0.99.467

## Contact Info

Guillem Salazar (salazar@icm.csic.es)
