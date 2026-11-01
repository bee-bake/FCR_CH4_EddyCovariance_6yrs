# FCR_CH4_EddyCovariance_6yrs
Scripts for the data processing and analysis of six years of eddy covariance CH₄ flux measurements at Falling Creek Reservoir at Vinton, Virginia, USA.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23864704.svg)](https://doi.org/10.5281/zenodo.23864704)

## Overview

This repository contains data processing scripts, statistical analyses, and visualization workflows for six years of methane (CH₄) flux measurements at Falling Creek Reservoir (FCR), Virginia, USA.

The project uses six years of near-continuous non-gap-filled eddy covariance measurements to investigate temporal variability in reservoir methane emissions and the environmental factors controlling CH₄ fluxes, with particular emphasis on the interactions between water temperature and dissolved oxygen availability.

## Research Objectives

The primary objectives are to:

- Characterize temporal variability in reservoir methane fluxes using eddy covariance measurements.
- Investigate the relationships between methane emissions, water temperature, dissolved oxygen (DO), and other environmental variables.
- Evaluate whether dissolved oxygen conditions modify the relationship between water temperature and methane fluxes.
- Identify important environmental predictors of methane flux variability using machine-learning approaches.

## Repository Contents

The repository is organized into the following components:

- **R/** – R scripts for data processing, statistical analyses, and figure generation.
- **Figure/** – Publication-quality figures generated from the analyses.
- **README.md** – Project description and documentation.
- **LICENSE** – MIT License for software reuse.
- **CITATION.cff** – Citation metadata for the repository.

*Note: The directory structure may be updated as the project develops.*

## Data Sources

The analyses use eddy covariance methane flux measurements and associated environmental observations from Falling Creek Reservoir.

Environmental variables include:

- Water temperature at multiple depths.
- Dissolved oxygen concentration and percent saturation.
- Wind speed and atmospheric conditions.
- Chlorophyll-a and other water-quality measurements.

Raw and supporting datasets may be maintained separately from this repository because of file size, data-management, or distribution restrictions. Refer to the relevant data repositories and publications for data access.

## Methods and Analytical Workflow

The analytical workflow includes:

1. Processing and quality control of half-hourly eddy covariance methane flux measurements.
2. Aggregation of methane flux observations into daily estimates.
3. Integration of methane fluxes with environmental and water-quality measurements.
4. Statistical modeling to evaluate temperature–oxygen interactions.
5. Machine-learning analysis using XGBoost to identify important environmental predictors.
6. SHAP analysis to interpret model predictions.
7. Permutation tests and sensitivity analyses to evaluate the temperature–oxygen interaction.
8. Generation of figures for scientific interpretation and publication.

## Software Requirements

The analyses are primarily conducted in R. See the individual scripts for the information about required packages for machine learning, SHAP analysis, statistical modeling, package dependencies and software requirements.

## Reproducibility
To reproduce the analyses:

1. Clone or download this repository.
2. Install the required R packages.
3. Obtain the required raw and supporting datasets from their respective sources.
4. Update file paths in the scripts as needed.
5. Run the scripts in the order indicated by their filenames or documentation.

## Citation
If you use the code or analytical workflows in this repository, please cite the archived version through Zenodo.

**Zenodo DOI:** To be added following archival.

## License
This repository's software is distributed under the MIT License. See the `LICENSE` file for details.

## Acknowledgments
This project began in February 2020 and was shaped by creative, interdisciplinary collaboration among researchers across institutions and international borders. We thank present and past reservoir group members for their assistance with field collection and laboratory analysis. Additionally, we thank the Western Virginia Water Authority (WVWA), especially Jamie Morris, for long-term access to field sites and logistical support. We gratefully acknowledge funding from U.S. National Science Foundation Grants DEB-1753639, DEB-2330211, RISE- 2438447, DEB-2213550, EF- 2318861, and DEB-2327030. We also acknowledge Discovery Grant support to Johnson provided by the Natural Sciences and Engineering Research Council of Canada (NSERC), RGPIN-2020-06252. 
