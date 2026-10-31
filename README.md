# Crime Near Transit: A Spatial Analysis of Metrorail Stations in Washington, DC

A project investigating the statistical relationship between public transit infrastructure and crime patterns in Washington, DC.

**Authors:** Isaac Damon, David McLone, Sky Reznik, and Ian Sargent
Department of Mathematics and Statistics, University of Vermont

## Overview

This project looks at how Washington, DC's Metrorail infrastructure relates to crime at the neighborhood level, using census tract panel data from 2020 to 2024. It has two parts:

1. **Statistical modeling.** Two fixed-effects negative binomial regressions, one for **violent** crime counts and one for **property** crime counts, across 12,024 tract-months. Each model relates crime rates to distance from the nearest Metro station, station characteristics (above-ground or underground, transfer, line corridor), and tract-level socioeconomic controls.
2. **Crime-weighted shortest path.** A walking-route algorithm built on the OpenStreetMap pedestrian network. It lets a user accept a longer walk in exchange for passing fewer reported crimes.

### Key findings

- Tracts **closer to Metro stations** have higher expected violent and property crime rates.
- Tracts whose nearest station is **above ground** have about **52% higher** expected crime rates than tracts whose nearest station is below ground. This holds for both crime types.
- Tracts nearest a **transfer station** have **64.8% lower** expected violent crime rates.
- Tracts on the **Green/Yellow** and **Blue/Orange/Silver** corridors have significantly higher violent crime rates than tracts on the Red Line (+72.3% and +46.8%).
- Unemployment, vacancy, and racial composition are also meaningful predictors.
- Across the tested routes, the crime-weighted router cut crime exposure by **23% to 100%**. The cost was between **5.73 m and 524.86 m** of extra walking.

All findings are associations, not causal effects.

## Hypotheses

The hypotheses draw on routine activity theory (Cohen & Felson, 1979) and earlier work linking transit to criminal opportunity (Phillips & Sandler, 2015; Groff, 2014):

- **H1:** Tracts closer to Metro stations have higher crime counts. *(Supported.)*
- **H2:** Proximity to Metro is more strongly associated with property crime than with violent crime. *(Mixed results.)*
- **H3:** Transfer stations and underground stations have stronger associations with crime than non-transfer and above-ground stations. *(Not supported. Above-ground stations had higher crime rates, and transfer stations had lower violent crime rates.)*

## Data

All data are publicly available and restricted to the District of Columbia boundary.

| Source | Description | Location in repo |
|---|---|---|
| [DC Open Data: Crime Incidents](https://opendata.dc.gov/datasets/DCGIS::crime-incidents-in-2020/about) | Geocoded incidents with offense, report date, method, and tract. Years 2017–2025 are included; 2020–2024 are used for modeling. | `Data/Crime/` |
| [DC Open Data: Metro Stations in DC](https://opendata.dc.gov/datasets/DCGIS::metro-stations-in-dc/about) | The 40 stations inside DC: name, line, address, and coordinates. Above-ground status, transfer status, and line corridor were labeled by hand. | `Data/Metro/Metros_Ground_Transfer.geoJson` |
| U.S. Census ACS 5-Year Estimates (2020–2024) | Tract-level income, poverty, unemployment, race/ethnicity, tenure, vacancy, and education. Pulled with `tidycensus`. | `Data/Census/` |

Model-ready panels are in `Data/Combined/`. `nb_data.csv` is the tract × year × month dataset, N = 12,024.

### Variable construction

- **Crime classification** follows FBI UCR categories:
  - *Violent:* homicide, assault with a dangerous weapon, sexual abuse, and robbery.
  - *Property:* burglary, theft (general), theft from auto, motor vehicle theft, and arson.
- **Aggregation:** incidents are spatially joined to 2020 census tracts with `sf`, then counted by tract × year × month. Tract-months with no incidents are set to 0. The totals are 20,655 violent crimes and 125,637 property crimes.
- **Metro proximity** is measured from each tract's centroid:
  - the Euclidean distance in meters to the nearest DC station;
  - the number of stations within 400 m;
  - the number of stations within 800 m.
- **Station characteristics** describe the nearest station: whether it is above ground, whether it is a transfer station, and its line corridor. The corridors are Green–Yellow, Blue–Orange–Silver, and Red (the reference category).

## Statistical Methodology

### Negative binomial regression

The crime counts are overdispersed (dispersion ≈ 2.83), so negative binomial regression is used instead of Poisson. Separate models are fit for violent and property crime with `fixest::fenegbin`. The exploratory models used `MASS::glm.nb`.

```math
\begin{aligned}
\ln\left(\frac{\mu_{itm}}{\text{Pop}_{itm}}\right) = {} & \beta_0 + \beta_1 \text{Metro}_{itm} + \beta_2 \text{SES}_{itm} + \beta_3 \text{Race}_{itm} + \beta_4 \text{Housing}_{itm} \\
& + \beta_5 (\text{Buffer}_{400} \times \text{Buffer}_{800})_{itm} + \beta_6 (\text{Dist} \times \text{PctBlack})_{itm} \\
& + \alpha_i + \gamma_t + \delta_m + \varepsilon_{itm}
\end{aligned}
```

| Term | What it contains |
|---|---|
| **Metro** | distance to the nearest station, above-ground or underground, transfer status, line corridor |
| **SES** | median household income, poverty rate, unemployment rate |
| **Race** | percent White, Black non-Hispanic, and Hispanic/Latino |
| **Housing** | percent renters, vacancy rate |
| **αᵢ, γₜ, δₘ** | tract, year, and month fixed effects, which absorb seasonality and COVID-19 disruptions to ridership and policing |
| **ln(Pop)** | population offset, so the coefficients describe crime *rates* |

The 400 m × 800 m buffer interaction is in the violent crime model only. The distance × % Black interaction is in the property crime model only. Standard errors are clustered by census tract.

### Feature selection

- **LASSO** (`mpath`, 10-fold CV), with the NB dispersion θ held fixed from an unpenalized fit:
  - violent model: kept 27 of 28 predictors (λ = 0.0003, θ = 3.02);
  - property model: kept 26 of 28 predictors (λ = 0.0008, θ = 2.74).
- **Random forest** (`randomForest`) variable importance, measured by %IncMSE. Distance to the nearest Metro ranked **#1** for both crime types.

### Diagnostics

Fit was checked with hanging rootograms (Kleiber & Zeileis, 2016):

- **Violent model:** fits well overall, with slight under-prediction above 10 crimes per tract-month.
- **Property model:** over-predicts counts of 4–12 and under-predicts counts above 18.

### Results: incidence rate ratios

| Variable | Violent IRR (SE) | Property IRR (SE) |
|---|---|---|
| Nearest Metro distance (m) | 0.9997\* (0.0001) | 0.9995\*\*\* (0.0001) |
| Stations within 400 m | 0.623\* (0.210) | 1.176 (0.143) |
| Stations within 800 m | 0.925 (0.069) | 1.160\* (0.069) |
| Above ground | 1.521\* (0.107) | 1.527\* (0.108) |
| Transfer station | 0.352\* (0.212) | 0.771 (0.221) |
| Green/Yellow line | 1.723\* (0.120) | 1.025 (0.112) |
| Blue/Orange/Silver line | 1.468\* (0.098) | 1.155 (0.128) |
| Median household income | 1.000\* (0.000) | 1.000\*\* (0.000) |
| Poverty rate | 1.002 (0.004) | 1.003 (0.004) |
| Unemployment rate | 1.013 (0.009) | 1.013 (1.010) |
| % Black non-Hispanic | 1.016\*\*\* (0.002) | 0.997 (0.003) |
| Percent renters | 1.008 (0.003) | 1.002 (0.003) |
| Vacancy rate | 1.009 (0.007) | 1.017\* (0.008) |
| 400 m × 800 m buffers | 1.510\*\*\* (0.007) | n/a |
| Metro distance × % Black | n/a | 1.001\*\*\* (0.000) |
| McFadden's pseudo-R² | .1054 | .0608 |
| χ² | 24.86\*\*\* | 20.46\*\*\* |

*Notes:* \*p < .05, \*\*p < .01, \*\*\*p < .001. The reference corridor is the Red Line. Both models include tract, year, and month fixed effects and a log-population offset, with SEs clustered by tract. "n/a" means the term is not in that model.

## Crime-Weighted Shortest Path

The algorithm is in [`Scripts/Analysis/Flow Net/safe_path copy.ipynb`](Scripts/Analysis/Flow%20Net/safe_path%20copy.ipynb). It is written in Python with OSMnx, NetworkX, GeoPandas, and Folium.

1. **Build the network.** Download the DC walking network with OSMnx, centered on the user's origin and destination.
2. **Filter by time shift.** The user gives a departure time, and only crimes from the same shift are used. The shifts are:

   | Shift | Hours |
   |---|---|
   | Graveyard | 00:00–03:59 |
   | Morning | 04:00–11:59 |
   | Afternoon | 12:00–15:59 |
   | Evening | 16:00–20:59 |
   | Night | 21:00–23:59 |

3. **Snap crimes to edges.** Each crime within `buffer_dist = 41.5` m of an edge is assigned to that edge. One crime can count toward several edges, for example both sides of a street.
4. **Weight the edges.** Each edge gets a crime-penalized weight:

   $$\text{CrimeWeight}_e = \ell_e \cdot \left(1 + \alpha \cdot \tilde{d}_e\right)$$

   Here $\ell_e$ is the edge length and $\tilde{d}_e$ is its normalized crime density. The parameter α ≥ 0 sets how strongly the router avoids crime; α = 0 gives the ordinary shortest path.
5. **Route.** Run Dijkstra's algorithm on `crime_weight`, and separately on `length` as a baseline.
6. **Score.** Crime exposure is the number of crimes along the route multiplied by (travel time ÷ shift length), at a walking speed of 1.31 m/s. The result is reported as the percent change between the shortest route and the crime-weighted route.

### Route results

**Route 1** runs from Lincoln Memorial to Farragut station, departing at 17:45, with α = 100. The crime-weighted route is 30.97 m longer than the shortest route and reduces crime exposure by **81.57%**.

**Routes 2–8** run from Foggy Bottom–GWU station to Georgetown University, departing at 01:30. The shortest route is 2,885.56 m.

| α | Crime-weighted route (m) | Distance increase | Crime exposure reduction |
|---|---|---|---|
| 254 | 3,410.42 | 18.19% | 100.00% |
| 200 | 3,165.58 | 9.70% | 97.53% |
| 100 | 3,115.93 | 7.98% | 96.36% |
| 60 | 3,052.74 | 5.79% | 94.06% |
| 20 | 3,007.61 | 4.23% | 90.63% |
| 10 | 2,963.73 | 2.71% | 82.69% |
| 2 | 2,891.29 | 0.20% | 23.44% |

Crime reduction rises quickly at small α and levels off at large α. The recommended range is α ≈ 20–200, with 100 as a good default for this route.

## Repository Structure

```
Data/
  Census/       ACS 5-year tract estimates (raw and cleaned)
  Combined/     model-ready panels (master.*, nb_data.csv, nb_data_by_shift.csv)
  Crime/        DC crime incidents 2017–2025 (GeoJSON) and yearly aggregates
  Metro/        DC stations with ground/transfer labels, regional lines, daily ridership
Scripts/
  Fetching/     Census_ACS_Fetch.Rmd: ACS pull via tidycensus
                Metro_Transfers_Ground.Rmd: station ground/transfer/corridor labels
  Cleaning/     Updated_DC_Cleaning.Rmd: builds nb_data.csv (tract × year × month panel)
                New master dataset.Rmd: incident-level master dataset
  Analysis/     Negative_Binomial_Modeling.Rmd: NB models, LASSO, random forest, rootograms
                Violent_vs_Property_Crime_EDA.Rmd: exploratory analysis and early glm.nb models
                Shift.Rmd: crime frequency by time-of-day shift
                time_crime.ipynb, plots.ipynb: temporal and spatial EDA (Python)
                Flow Net/: crime-weighted shortest path notebooks
```

## Reproducing the Analysis

**R** (≥ 4.x). Packages are loaded with `pacman`:

```r
install.packages(c("pacman", "tidyverse", "sf", "tidycensus", "fixest", "MASS",
                   "mpath", "randomForest", "janitor", "lubridate"))
```

The ACS fetch needs a free [Census API key](https://api.census.gov/data/key_signup.html), set with `tidycensus::census_api_key()`.

**Python** (3.10+):

```bash
pip install geopandas pandas numpy networkx osmnx folium matplotlib seaborn scikit-learn
```

Run the scripts in this order:

1. `Fetching/`
2. `Cleaning/Updated_DC_Cleaning.Rmd`
3. `Analysis/Negative_Binomial_Modeling.Rmd`

The modeling scripts read `nb_data.csv` directly from this repository's GitHub URL, so step 3 also works on its own. The routing notebook prompts for an origin, a destination, and a departure time.

## Limitations

- **Associational only.** The results do not show that Metro proximity causes crime, and they do not show that a suggested route is safe.
- **DC-specific.** The findings may not generalize to other transit systems.
- **Time window.** Tract boundaries changed in 2020, so the study covers only 2020–2024, which largely overlaps the COVID-19 pandemic.
- **Omitted variables.** Police activity, land use, surveillance, and ridership are not modeled.
- **Independence assumption.** Spatial and temporal clustering of crime may violate it.
- **Straight-line distance.** Proximity is Euclidean, so it ignores walking paths and barriers such as rivers and highways.
- **Unexplained variation.** McFadden's R² is modest, so much of the variation in crime is not explained by the models.

## Contributors

All four authors helped conceptualize the project and contributed to every part of the work.

| Contributor | Analysis | Writing |
|---|---|---|
| **Isaac Damon** | Data cleaning and feature selection | Figures and tables |
| **Ian Sargent** | Data cleaning and feature selection | Introduction and methods |
| **David McLone** | Negative binomial models | Results |
| **Sky Reznik** | Crime-weighted path algorithm | Crime-weighted path sections |

## Acknowledgments

Thanks to Professor Buzas for feedback and guidance throughout this capstone project. Generative AI (Anthropic's Claude) was used for data cleaning and LaTeX manuscript preparation. All interpretations are the authors' own.

## License

[MIT](LICENSE) © 2026 Ian Sargent
