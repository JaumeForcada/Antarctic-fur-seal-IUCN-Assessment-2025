---
output:
  html_document: default
  pdf_document: default
---
# Antarctic-fur-seal-IUCN-Assessment-2025 - CODE

Functions used to obtain fur seal abundant estimates, trends and population
reduction for the Antarctic fur seal 2025 IUCN Red List assessment

## Description
Analysis functions in R and BUGS code:

1. R script **SSBFemale_IPM_2001-2025.R** fits an IPM (Integrated Population Model) to FEMALE count, demographic, and mark-recapture data collected at the Special Syudy Beach (SSB) of Bird Island from 2000/2001
to 2024/2025. Generates data outputs: **AFSFemaleSimulations200125.csv**, **AFSFemaleSimulations200125.Rds**, and **SSBFemaleIPM200125.RData**

2. R script **SSBMale_IPM_1995-2025.R** fits an IPM (Integrated Population Model) to MALE count, demographic, and mark-recapture data collected at the Special Syudy Beach (SSB) of Bird Island from 1994/1995
to 2024/2025. Generates data outputs: **SSBMaleIPM19942025.RData**, **AFSMaleSimulations19952025.csv**, and **AFSMaleSimulations19952025.Rds**.

3. R script **SSBMatureFemale_IPM_2001-2025.R** fits an IPM (Integrated Population Model) to MATURE FEMALE counts, demographic, and mark-recapture data collected at the Special Syudy Beach (SSB) of Bird Island from 1978/1979
to 2024/2025. Generates following data outputs: **AFSFemaleSimulations19842025.Rds**, **AFSFemaleSimulations19842025.csv**, and **SSBFemaleIPM19842025.RData**

4. R script **AFSFemale_vital_rates.R** to derive estimates of vital rates (survival, recruitment, fecundity, immigration), life tables, and generation time (length) for Antarctic fur seal females

5. R script **AFSMale_vital_rates.R** to derive estimates of vital rates (survival, recruitment (becoming territorial), immigration), life tables, and generation time (length) for Antarctic fur seal males

6. R script **SSBFemalePMCountsForTrendAnalysis.R** fits Bayesian GAMMs to SSB female daily counts for prediction of SSB count value on the day a survey count is obtained. 

7. R script  **SSBMaleAMCountsForTrendAnalysis.R** does the same for AM territorial male counts at SSB.

8. R script **BIMatureFemaleAbundance** estimates abundance of mature females (breeding and non-breeding) at Bird Island from complete island counts and correction factors from SSB derived from the female IPM.

9. R script **BIMatureMaleAbundance.R** estimates abundance of mature males (territorial and non-territorial) at Bird Island from complete island counts and correction factors from SSB derived from the male IPM.

10. R script $\texttt{PupProductionMultipliers.R}$ estimates multipliers of total pup production to obtain estimates of mature population size and total population size. Outputs include $\texttt{PtoM.Rds}$ and $\texttt{PtoT.Rds}$ which are vectors with 15,000 simulations.

5. Support function **glbpf.R** is sourced to fit a low pass Gaussian filter to summarise trends.

6. Support function **addTrans.R** is sourced to add colour transparency in summary plots.

7. Suport function **estBetaParams** to estimate shape parameters for Beta distribution priors.

8. Function **ffs.obs.R** fits non-linear least-squares models with logistic structure to estimate median peak pupping, variation and total pup production to a time series of cumulative daily counts of new pups born.
$\texttt{rjags}$

## Getting Started

### Dependencies

* Requires [R](https://cran.r-project.org/)

* Requires [JAGS software](https://sourceforge.net/projects/mcmc-jags/) installed to fit MCMC models. IPMs are coded in BUGS language, which is called from the R scripts using appropriate R packages, as described in the analysis scripts.

### Reproducing the analysis

* How to run the program
* Step-by-step bullets
```
code blocks for commands
```

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors

Contributors names and contact info

ex. Dominique Pizzie  
ex. [@DomPizzie](https://twitter.com/dompizzie)

## Version History

* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)