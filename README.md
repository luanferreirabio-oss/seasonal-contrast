# From River to Soil — data and analysis code

Data and analysis scripts supporting the manuscript:

> **From River to Soil: Seasonal Contrast and Cation Fractionation in a Floodplain Cacao System, Lower Tocantins, Pará**
> Ferreira LDS, Cordeiro LC, Berrêdo da Silva JFB, Vilhena MPSP, Cárdenas Alegría OV, Ramos RTJ, Sousa MPA
> Submitted to *Journal of Soils and Sediments*.

The study integrates three compartments of the lower Tocantins River (Pará, Brazil): river water sampled hourly through tidal cycles at two stations in two seasons, soil chemistry from six cacao-cultivating floodplain islands, and the genetic potential for cation transport obtained by shotgun metagenomics of those islands.

---

## Contents

### `data/`

| File | Rows | What it holds |
|---|---|---|
| `01_water_51_samples.csv` | 51 | Hourly surface-water samples, two stations, two campaigns, 15 retained variables. Values below the limit of quantification are recorded as zero, with no imputation. Two outlying anion records from campaign 1 were excluded, as described in the Methods. |
| `02_soil_chemistry_6islands_2depths.csv` | 12 | Soil chemical attributes, 6 islands × 2 depths (0–10 and 10–20 cm). Exchangeable aluminium is **not** included, for the reasons given in the Methods and in the Discussion's limitations. |
| `03_spearman_correlations_105_pairs.csv` | 105 | Spearman correlations between the 15 water constituents, computed on the pooled set (n = 51) and separately within each campaign, with p values and Benjamini-Hochberg q values. |
| `04_guild_composition_curated_COG2024.csv` | 66 | Composition of the four cation transport guilds: COG identifier, functional category, gene, curated COG-2024 description, and whether the group was detected. |
| `05_COG_abundance_matrix_by_island.csv` | 2746 | Full abundance matrix, 2746 orthologous groups × 6 islands, in copies per million. |

**Island identifiers.** The data files carry the field codes `P1`–`P6`, which map to the islands as follows:

| Field code | Island |
|---|---|
| P1 | Angapijó |
| P2 | Conceição |
| P3 | Santana |
| P4 | Santaninha |
| P5 | São Joaquim |
| P6 | Tauaré |

Always join datasets by island **name**, not by code or position.

### `scripts/`

| Script | What it does |
|---|---|
| `10_dados.py` | Builds and internally checks the water database |
| `11_agua_final.py` | Detection percentages, intratidal CV, station and campaign contrasts, FDR correction |
| `12_pca_corr.py` | Principal component analysis and pooled vs within-campaign correlations |
| `20_solo.py` | Soil attributes, molar ratios, Gapon coefficient, apparent selectivities |
| `30_guildas.py` | Guild assembly, abundances, exact p values by complete enumeration |


---

## Raw sequencing reads

Raw reads are deposited in the Sequence Read Archive (SRA) of the National Center for Biotechnology Information (NCBI) under BioProject **PRJNA1224407**. They are not duplicated here.

## Genetic heritage

Access to genetic heritage was registered in the Brazilian National System for the Management of Genetic Heritage and Associated Traditional Knowledge (SisGen) under number A2AD73B, in accordance with Law 13.123/2015.

## Reproducing the analyses

```bash
python -m pip install pandas numpy scipy scikit-learn matplotlib
python scripts/10_dados.py
python scripts/11_agua_final.py
python scripts/12_pca_corr.py
python scripts/20_solo.py
python scripts/30_guildas.py
```

The scripts read from `data/` and were written for Python 3.11. Paths at the top of each script may need adjusting to your local layout.

## A note on the soil data

The soil chemistry comes from a single consolidated analytical record, without the bench readings behind it. Two consequences are stated in the manuscript and repeated here so that anyone reusing these data is aware of them:

1. Exchangeable aluminium was excluded because, across all units, the recorded value is a fixed multiple of the calcium content, with a residual smaller than the decimal place to which results were reported. The column is derived, not independently measured.
2. The between-island dispersion of calcium and magnesium in this record is very low and should be treated as indicative, not as an estimate of real variability between islands. The conclusions of the paper rest on order-of-magnitude contrasts (3.9–19.9× in Ca/K between water and soil), which are unaffected.


