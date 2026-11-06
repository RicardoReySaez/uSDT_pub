# A Multilevel Signal Detection Theory Approach to Unconscious Processing Research

This repository contains the data, R scripts, results, and Quarto
sources needed to reproduce the manuscript and its Supplemental Material:

> Rey-Sáez, R., Franco-Martínez, A., Vadillo, M. A., & Garre-Frutos, F.
> (2026). *A multilevel signal detection theory approach to unconscious
> processing research*. Manuscript in preparation.

The large files (fitted Bayesian models and raw simulation outputs) are
archived on Zenodo:
[doi.org/10.5281/zenodo.24019945](https://doi.org/10.5281/zenodo.24019945).

## The `{uSDT}` R package

All models in this work are fitted with `{uSDT}`, an R package for
hierarchical signal detection theory models of unconscious processing.
If you want to apply the method to your own data, start there:

- Website and documentation: <https://ricardoreysaez.github.io/uSDT/>
- Source code: <https://github.com/RicardoReySaez/uSDT>
- CRAN: <https://CRAN.R-project.org/package=uSDT>

```r
install.packages("uSDT")                 # released version
# pak::pak("RicardoReySaez/uSDT@dev")    # development version
```

Some scripts in this repository were last run with the development
version of `{uSDT}` (0.1.0.9000) and include the line to install it.

## Repository structure

```
uSDT_pub/
├── uSDT_pub.Rproj        RStudio project; all scripts run from this folder
├── Data/                 Trial-level data of the three empirical studies
│   ├── Greenwald 2017/   Greenwald & De Houwer (2017), from osf.io/ydeza
│   ├── Vadillo 2024/     Vadillo, Malejka, & Shanks (2025), from osf.io/jp3gx
│   └── Wentura 2025/     Wentura, Rohr, & Kiefer (2025), from osf.io/9t4pk
├── R scripts/
│   ├── Empirical illustrations/   Reanalyses of the three empirical studies
│   ├── Extra analysis/            Reanalysis of all 18 data sets of Greenwald & De Houwer (2017)
│   ├── Manuscript figures/        Figures 1-4 (illustrative figures)
│   ├── Power simulation/          Figure 5 (power by number of trials)
│   ├── Simulation study/          Simulation Studies 1 (validation) and 2 (robustness)
│   └── R functions/               Helper functions sourced by the scripts above
├── Results/
│   ├── Empirical illustrations/   Summaries of the empirical models (.rds)
│   └── Simulation studies/
│       ├── Processed/             Tables and long data read by the manuscript
│       └── Raw/                   Raw SimDesign outputs (on Zenodo, not on GitHub)
├── Fits/                 Fitted Bayesian models (on Zenodo, not on GitHub)
├── Figures/
│   ├── Manuscript/       Figure_1.png ... Figure_11.png
│   └── Supplemental/     Figure_S1.png ... Figure_S11.png
└── Manuscript/
    ├── Manuscript.qmd     Main text (Quarto + apaquarto)
    ├── supplemental.qmd   Supplemental Material
    ├── uSDT_bib.bib       References
    └── _extensions/       apaquarto extension used to render both documents
```

### Scripts and the outputs they produce

| Script | Produces |
|---|---|
| `Empirical illustrations/1. greenwald2017_prepare_data.R` | Downloads the OSF archive and builds `Data/Greenwald 2017/greenwald2017_long.csv` |
| `Empirical illustrations/2. greenwald2017_exp1c_analysis.R` | Figures 8 and 9; results and fit of the masked conditioning study |
| `Empirical illustrations/3. vadillo2025_exp2_analysis.R` | Figure 10; results and fit of the contextual cuing study |
| `Empirical illustrations/4. wentura2025_exp2a_analysis.R` | Figure 11; results and fit of the affective priming study |
| `Extra analysis/greenwald2017_series1_analysis.R` | Tables of the three hypotheses in the 18 data sets of Series 1 (footnote in the masked conditioning example); needs script 1 first |
| `Manuscript figures/manuscript_figures.R` | Figures 1-4 |
| `Power simulation/power_trials_direct.R` | Figure 5 |
| `Simulation study/uSDT_validation_sim.R` | Raw outputs of Simulation Study 1 |
| `Simulation study/uSDT_validation_analysis.R` | Figures 6, 7, and S1-S8; processed results of Simulation Study 1 |
| `Simulation study/uSDT_robustness_sim.R` | Raw outputs of Simulation Study 2 |
| `Simulation study/uSDT_robustness_analysis.R` | Figures S9-S11; processed results of Simulation Study 2 |

## How to reproduce the results

1. Open `uSDT_pub.Rproj` in RStudio (or set the working directory to
   this folder). Every script uses paths relative to the project root.
2. **To render the manuscript only**, nothing else is needed: the
   figures and the processed results are already in the repository.
   Render `Manuscript/Manuscript.qmd` and `Manuscript/supplemental.qmd`
   with [Quarto](https://quarto.org/); the apaquarto extension is
   included in `Manuscript/_extensions/`.
3. **To rerun the analyses**, run the scripts in the table above. The
   empirical scripts must run in their numbered order. The simulation
   scripts (`*_sim.R`) are computationally expensive; to skip them,
   download the raw outputs from Zenodo into
   `Results/Simulation studies/Raw/` and run only the `*_analysis.R`
   scripts.
4. To inspect the fitted Bayesian models without refitting them,
   download them from Zenodo into `Fits/`.

The results were obtained with R 4.5.0 and Quarto 1.10.

## Contact

Ricardo Rey-Sáez (<ricardo.rey@uam.es>) and Francisco Garre-Frutos
(<francisco.garre@uam.es>), Autonomous University of Madrid.
