---
title: 'epimutacions: a package to detect epimutations in the context of rare diseases'
tags:
  - biohackeu20
  - DNA methylation
  - Rare diseases
authors:
  - name: Carlos Ruiz-Arenas
    orcid: 0000-0002-6014-3498
    affiliation: 1
  - name: Carles Hernandez-Ferrer
    orcid: 0000-0002-8029-7160
    affiliation: 2
  - name: Alejandro Caceres
    affiliation: 3
  - name: James Baye
    orcid: 0000-0002-0078-3688
    affiliation: 4
  - name: Leire Abarrategui
    orcid: 0000-0002-1175-038X
    affiliation: 5
  - name: Lordstrong Akano
    orcid: 0000-0002-1404-0295
    affiliation: 6
  - name: Pavlo Hrab
    affiliation: 7
    orcid: 0000-0002-0742-8478
  - name: Raquel Manzano
    affiliation: 8
    orcid: 0000-0002-5124-8992
  - name: Margherita Mutarelli
    orcid: 0000-0002-2168-5059
    affiliation: 9
affiliations:
  - name: Centro de Investigación Biomédica en Red de Enfermedades Raras (CIBERER), Barcelona, Spain
    index: 1
  - name: Centro Nacional de Análisis Genómico (CNAG-CRG), Center for Genomic, Regulation
    index: 2
  - name: ISGlobal, Institute for Global Health, Dr Aiguader 88, 08003 Barcelona, Spain
    index: 3
  - name: Wellcome/MRC Cambridge Stem Cell Institute, University of Cambridge, Cambridge CB2 0AW, UK
    index: 4
  - name: Faculty of Medical Sciences, Newcastle University, Newcastle-Upon-Tyne, UK
    index: 5
  - name: College of Medicine, University of Ibadan
    index: 6
  - name: Department of Genetics and Biotechnology, Biology faculty, Ivan Franko National University of Lviv
    index: 7
  - name: Cancer Research UK Cambridge Institute
    index: 8
  - name: Institute of Applied Sciences and Intelligent Systems (ISASI-CNR)
    index: 9
date: 13 November 2020
bibliography: report.bib
event: biohackeu20
---

# Background

Rare diseases are pathologies with low prevalence (< 1 per 2,000 people). Most of these pathologies have an onset during childhood and a strong genetic etiology. Consequently, rare disease diagnosis has relied on identifying genetic and genomic mutations that can cause the disease. Although these variants have provided a diagnosis for many patients and families, around 60% of the cases still remained undiagnosed [@Lionel2018]. Aberrant methylation can be an underlying cause of undiagnosed patients, either as a primary event (a.k.a. epimutation) or as a functional consequence of chromatin dysregulation by genetic or environmental agents (a.k.a. episignature). Epimutations are the cause of some rare
diseases, such as Prader-Willi, Angelman or Beckwith-Wiedemann syndromes [@Aref-Eshghi2019; @Inoue2017] and some human malformations [@Serra-Juhe2015]. Syndrome-specific episignatures are increasingly defined as biomarkers for a growing number of disorders [@Aref-Eshghi2020]. Therefore, tools to detect epimutations and episignatures should be made available to the rare disease community and included in standardized analysis workflows.

During the European Biohackathon 2020, we have been working on implementing an R package called `epimutacions` to enable analyzing epimutations to the rare disease community. We have been working on three main areas: (1) implementation of algorithms to detect epimutations; (2) implementation of reporting functions; (3) testing and validation of the package.

# Hackathon Results

## Implementation of algorithms to detect epimutations

We have implemented five different algorithms to detect epimutations. Two of these algorithms were previously described in the literature [@Barbosa2018; @Aref-Eshghi2019], while the other three are adaptations from the algorithm of Aref-Eshghi et al. These algorithms have as input a proband whose DNA methylation is compared with a reference group.

### Barbosa 

We called Barbosa to the method described in [@Barbosa2018]. Briefly, the algorithm checks for each CpG, if the proband’s measurement is an outlier. Then, it calls an epimutation to those regions where 3 contiguous CpGs are outliers and they are separated by less than 500 base pairs.

### MANOVA

We called MANOVA to the implementation of [@Aref-Eshghi2019]. This algorithm starts by running bumphunter, a method to detect DMRs [@Peters2015], comparing the proband versus the reference group. This first step results in a list of DMRs, regions where the proband has different methylation than the reference group. Then, each DMR is tested for outlier using a MANOVA (Multivariate ANalysis Of VAriance). DMRs are then subsetted based on the F-statistic magnitude and mean difference between proband and reference group DNA methylation.

### Additional implementations

From [@Aref-Eshghi2019] implementation, we proposed three different approaches. These three approaches run the bumphunter step but test DMR for outliers using different approaches: (1) mlm, that uses a multivariate linear model; (2) iso.forest, that uses isolation forest; and (3) Mahdist.MCD, which uses robust Mahalanobis distance.

### Main function implementation

Our package `epimutacions` has a main function (`epimutations`) that is able to run the desired algorithm to detect epimutations to our dataset. `epimutations` has as input a `GenomicRatioSet` and returns a data.frame with the detected epimutations and statistics returned by the methods. `epimutations` is composed of modular functions that encapsulate each step of the method. This design is very flexible and facilitates adding additional algorithms in the future.

## Reporting functions

We have worked on reporting user-friendly outputs. Our work was focused on adding biological information that might improve the interpretation of the results. Thus, we implemented a function that maps the epimutations to genes and adds any overlapping ENSEMBL regulatory region. This function also checks if the mapped genes are in OMIM. We also have started to work on a plotting function to represent the epimutation. 

## Testing and package validation

We have developed some toy objects, which can be used to test the package and exemplify how the functions worked. A subgroup of the team applied the package to new datasets in order to test it in a real environment and propose improvements to the implementation subgroup. 

# Conclusion

We have a functional version of an R package to detect epimutations from microarray DNA methylation data. This package will help researchers from the rare diseases’ field to incorporate epimutations discovery in their workflows.

# Future work

## Implementation

We are planning to submit the package to Bioconductor. To this end, the code should be modified to comply with Bioconductor requirements. In addition, tests and examples will be added as well as a fully functional vignette. 

## Reporting
Additional work is required to improve plotting functions. The annotation function can be improved to accept additional annotation packages and add more information about the genes. 

## Package evaluation

During the biohackathon, we focused our efforts on testing whether the algorithms returned their expected value. However, a comparison and evaluation of the methods’ performance should be done. In this context, different approximations to preprocess the data should be considered.


# Acknowledgements

We acknowledge the organizers of Europe Biohackathon 2020 for their support.

# References

