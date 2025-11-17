---
editor_options: 
  markdown: 
    wrap: 72
bibliography: references.bib
---

Please send issues and questions to:
[enrico.e.armelloni\@gmail.com](mailto:enrico.e.armelloni@gmail.com){.email}

# Description

The repository contains the data and code supporting the manuscript
Sguotti C., Armelloni, E. N., Lizzi S., Scarcella, G.. From
overexploited to overshoot: how spatial protection, reduced fishing
pressure, and environmental conditions triggered a positive tipping
point in the Adriatic Common sole. Submitted for consideration to ICES
JMS. The code and the data here reported serves to replicate the
analysis.

## Analysis summary

| Research Question | Dataset | Script | Analysis |
|----|----|----|----|
| Q1; Q2/M1; Q2/M2 | model_2d_input.csv | analysis_2d | Tipping point detection protocol: change point analysis; state driver plots, bimodality, linear model; threshold gam. |
| Q2/M3 |  model_3d_input.csv | analysis_3d_model_pa; analysis_3d_model_pos | Spatial analysis: gamm, hotspot persistence |

Description of data, scripts and methodologies used to answers the
research questions from Sguotti et al. (submitted). For details upon the
research questions refer to the manuscript text.

## R scripts folder

Contains the code used to perform the analysis, as described in Table 1
of the present repo, and the images included in the manuscript.

The scripts to replicate results are named with logic order (step1,
step2 etc.). Each file contains a short note describing the purpose of
the code. *supporting_functions.R* file contains custom functions that
are used in the other scripts. *Appendix_S1.R* contains the code used
for the simulation study described in Appendix_A

-   analysis_2d: code to replicate analysis based on 2d dataset and all
    the figures showing results;

-   analysis_3d_model_pa; analysis_3d_model_pos: code to replicate
    spatial analysis

-   analysis_3d_plots: code to generate all the images showing results
    of spatial analysis

-   figure1: generates the subpanels of figure 1

-   HighstatLibV11: code from [@zuur2017] containing functions to
    perform data inspection.
