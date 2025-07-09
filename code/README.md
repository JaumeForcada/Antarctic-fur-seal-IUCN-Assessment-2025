---
output:
  html_document: default
  pdf_document: default
---
# Antarctic-fur-seal-IUCN-Assessment-2025 - CODE

R Scripts used to obtain Antarctic fur seal abundant estimates, trends, and population change and/or reduction for the Antarctic fur seal 2025 as quantitative support of the IUCN Red List Assessment of the species for 2025.

## Description
Analysis functions in R and BUGS code:

1. R script **SSBFemale_IPM_2001-2025.R** fits an IPM (Integrated Population Model) to FEMALE count, demographic, and mark-recapture data collected at the Special Syudy Beach (SSB) of Bird Island from 2000/2001
to 2024/2025. Generates data outputs: **AFSFemaleSimulations200125.csv**, **AFSFemaleSimulations200125.Rds**, and **SSBFemaleIPM200125.RData**

2. R script $\texttt{SSBMale_IPM_1995-2025.R}$ fits an IPM (Integrated Population Model) to MALE count, demographic, and mark-recapture data collected at the Special Syudy Beach (SSB) of Bird Island from 1994/1995
to 2024/2025. Generates data outputs: **SSBMaleIPM19942025.RData**, **AFSMaleSimulations19952025.csv**, and **AFSMaleSimulations19952025.Rds**.

3. R script $\texttt{SSBMatureFemale_IPM_1984-2025.R}$ fits an IPM (Integrated Population Model) to MATURE FEMALE counts, demographic, and mark-recapture data collected at the Special Syudy Beach (SSB) of Bird Island from 1978/1979
to 2024/2025. Generates following data outputs: **AFSFemaleSimulations19842025.Rds**, **AFSFemaleSimulations19842025.csv**, and **SSBFemaleIPM19842025.RData**

4. R script $\texttt{AFSFemale_vital_rates.R}$ to derive estimates of vital rates (survival, recruitment, fecundity, immigration), life tables, and generation time (length) for Antarctic fur seal females

5. R script $\texttt{AFSMale_vital_rates.R}$ to derive estimates of vital rates (survival, recruitment (becoming territorial), immigration), life tables, and generation time (length) for Antarctic fur seal males

6. R script $\texttt{SSBFemalePMCountsForTrendAnalysis.R}$ fits Bayesian GAMMs to SSB female daily counts for prediction of SSB count value on the day a survey count is obtained. 

7. R script  $\texttt{SSBMaleAMCountsForTrendAnalysis.R}$ does the same for AM territorial male counts at SSB.

8. R script $\texttt{BIMatureFemaleAbundance}$ estimates abundance of mature females (breeding and non-breeding) at Bird Island from complete island counts and correction factors from SSB derived from the female IPM.

9. R script $\texttt{BIMatureMaleAbundance.R}$ estimates abundance of mature males (territorial and non-territorial) at Bird Island from complete island counts and correction factors from SSB derived from the male IPM.

10. R script $\texttt{PupProductionMultipliers.R}$ estimates multipliers of total pup production to obtain estimates of mature population size and total population size. Outputs include $\texttt{PtoM.Rds}$ and $\texttt{PtoT.Rds}$ which are vectors with 15,000 simulations.

11. R script $\texttt{SouthGeorgiaPopulationAssessment.R}$ projects the abundance estimates for South Georgia obtained in 2007/9 to total estimates in 2022, based on survey data obtained on that year.

12. R script $\texttt{CriterionA2Assessment.R}$ sources abundance simulations for South Georgia in 2007/9 AND 2022, generated with script 
$\texttt{SouthGeorgiaPopulationAssessment.R}$; obtains abundance simulations for other subpopulations from pup production data; estimates population change and reduction under criteria A1/A2.

13. R script $\texttt{AFSDynamicsSigny.R}$ reproduces the population dynamics analysis of Signy Island population index, including synchrony analysis.

Support R functions called by the some of the scripts listed above:

1. Support function $\texttt{glbpf.R}$ is sourced to fit a low pass Gaussian filter to summarise trends.

2. Support function $\texttt{addTrans.R}$ is sourced to add colour transparency in summary plots.

3. Suport function $\texttt{estBetaParams}$ to estimate shape parameters for Beta distribution priors.

4. Function $\texttt{ffs.obs.R}$ fits non-linear least-squares models with logistic structure to estimate median peak pupping, variation and total pup production to a time series of cumulative daily counts of new pups born.
$\texttt{rjags}$

5. Function $\texttt{getFemales.R}$ estimates mature female abundance from a single or vector of female count(s) during the helicopter survey of South Georgia from 2009.
$\texttt{rjags}$

6. Function $\texttt{getFemales07.R}$ estimates mature female abundance from a single or vector of female count(s) during the helicopter survey of South Georgia from 2007.
$\texttt{rjags}$

7. Function $\texttt{getMales.R}$ estimates mature male abundance from a single or vector of male count(s) during the helicopter survey of South Georgia from 2009.
$\texttt{rjags}$

8. Function $\texttt{getMales07.R}$ estimates mature male abundance from a single or vector of male count(s) during the helicopter survey of South Georgia from 2007.
$\texttt{rjags}$

## Getting Started

### Dependencies

* Requires [R](https://cran.r-project.org/)

* Requires [JAGS software](https://sourceforge.net/projects/mcmc-jags/) installed to fit MCMC models. IPMs are coded in BUGS language, which is called from the R scripts using appropriate R packages, as described in the analysis scripts.

### Reproducing the analysis

* Run sequentially the R scripts 1 to 12 above, and optionally 13.
Some of the scripts which estimate abundance at SSB-BI including $\texttt{SSBFemale_IPM_2001-2025.R}$, $\texttt{SSBMale_IPM_1995-2025.R}$, and $\texttt{SSBMatureFemale_IPM_1984-2025.R}$ take days to complete. It is recommended to skip those unless required, as the data outputs from the scripts are available here.

Some of the data files are not available at the github repository due to large size. Those are included in the linked Zenodo repository, and include:

1. **AFSDynamicsSigny.RData**

2. **AFSFemaleEarlySimulations19842025.Rds**

3. **AFSFemaleSimulations200122Best.Rds**

4. **AFSFemaleSimulations200125.csv**

5. **AFSFemaleSimulations200125.Rds**

6. **AFSFemaleSimulations19842025.csv**

7. **AFSFemaleSimulations19842025.Rds**

8. **AFSMaleSimulations19952007Best.Rds**

9. **AFSMaleSimulations19952025.csv**

10. **AFSMaleSimulations19952025.Rds**


## Help

Any advise for common problems or issues:
jfor@bas.ac.uk

## Author

Dr Jaume Forcada, 
(https://www.bas.ac.uk/profile/jfor/)

## Version History
* 0.1
    * Initial Release

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments
Contributions:
  <span style="text-decoration:underline">Data support and project management</span> -- All the Antarctic fur seal data from Bird Island, South Georgia, and Signy Island, South Orkney Islands used in this analysis are part of the British Antarctic Survey (BAS) core science programme; it was generated by zoological field assistants and scientists working for BAS since 1978; and is curated by the BAS Polar Data Centre.
  
  <span style="text-decoration:underline">Genetic recapture data</span> -- Bird Island samples were analysed and contributed by the J. Hoffman lab at University of Bielefeld, Germany, especially by J. Hoffman and A. Paijmans in collaboration with BAS.
  
  <span style="text-decoration:underline">Additional data and project management</span> --  South Georgia mainland, South Sandwich Islands, and Signy Island: Bucktrout, P., Collins, M., Coleman, J., Dunn, M., Dickens, J., Fenney, N., Fox, A., Hollyman, P., and Ratcliffe, N., Wood, A. G.
  
  <span style="text-decoration:underline">Data support from other locations</span> -- de Bruyn, N., Delord, K., Guinet, C., Goldsworthy, S., Jordaan, R., Krause, D. J., Lea, M.-A., and Lowther, A. Below, section 'Pup count data for other subpopulations' lists additional references when data were obtained from the literature.
