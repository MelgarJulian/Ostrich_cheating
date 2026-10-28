# Frequency-dependent cheating in a cooperative breeding bird

Data, analysis code and manuscript sources for a study of frequency-dependent
cheating in cooperatively breeding ostriches, using the Klein Karoo experimental
population (2012–2018).

Everything in the manuscript and supplementary information is generated from
this repository: the tables and figures are built by knitting the two `.Rmd`
files, and the numbers quoted in the text are inline R expressions rather than
typed values.

## Running order

The scripts form a chain, and each one depends on the output of the one before.

| # | Script                                        | Reads                                                                             | Writes                                                                                                                                                                                            |
| - | --------------------------------------------- | --------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1 | `Scripts/Cheating - Data preparation.R`       | `R_Input_files/*.csv`                                                             | `R_Output_files/Cheating_data.RData`                                                                                                                                                              |
| 2 | `Scripts/Cheating - Model data.R`             | `R_output_files/Cheating_data.RData`                                              | `R_Output_files/Cheating_Model_Data.RData`                                                                                                                                                        |
| 3 | `Scripts/Cheating_MCMCglmmModels.R`           | both of the above                                                                 | `R_Output_files/CheatingModels.RData`                                                                                                                                                             |
| 4 | `Scripts/Frequency_dependent_cheating_SI.Rmd` | `R_output_files/CheatingModels.RData`, `R_output_files/Cheating_Model_Data.RData` | `Scripts/Manuscripts/Frequency_dependent_cheating_SI_[date].docx` , `Scripts/Manuscripts/Frequency_dependent_cheating_SupplementaryTables.xlsx`, and `R_Output_files/Cheating_SInformation.RData` |
| 5 | `Scripts/Frequency_dependent_cheating_MS.Rmd` | `R_Output_files/Cheating_SInformation.RData`                                      | `Scripts/Manuscripts/Frequency_dependent_cheating_MS_[date].docx`                                                                                                                                 |

**The Extended Data (Frequency_dependent_cheating_SI.Rmd) must be knitted before the manuscript.** 
Step 4 writes the `.RData` that step 5 loads, so knitting the manuscript against a stale one silently
reports the previous run's numbers.

Step 3 takes about 23 hours, most of it in `m_laying_window`, which has 16258
rows and accounts for roughly two thirds of the total. The fitted models are
committed, so steps 4 and 5 can be run without repeating it.

These five scripts are the whole pipeline. `Scripts/` holds nothing else that
runs.

## Naming

Models and the data frames they are fitted to are named for what they estimate,
never for their position in the manuscript, so that the paper can be
restructured without renaming code:

- models are `m_<topic>_<structure>_<component>`, for example `m_eggs_quad_c2`
- data frames are `df_<topic>`, so `m_eggs` is fitted to `df_eggs`

The full list, and the mapping from each model to its manuscript number and
supplementary table, is at the top of
`Scripts/Cheating_MCMCglmmModels.R`. The corresponding list of data
frames is at the top of `Scripts/Cheating - Model data.R`.
Manuscript and table numbers are recorded only in those two guides, because they
change whenever the paper is reorganised.

## Reproducing the results

Every step is reproducible. Running steps 1 and 2 rebuilds each data frame
identically to the version committed here, and every model is seeded, so
refitting reproduces the published estimates exactly.

All models use one set of MCMC settings, defined once near the top of
`Cheating_MCMCglmmModels.R` and reported in the manuscript methods:

```r
nitts <- 2000000
burn  <- 200000
thin  <- 500
```

This retains 3600 samples per chain and 10800 across the three chains. The
values can be read back off any fitted model with `attr(model$VCV, "mcpar")`.

Each of the 100 fits is preceded by its own `set.seed()`, so a single model can
be re-run on its own and reproduce exactly. The seed encodes the model as
`topic * 100 + structure * 10 + chain`, with topics numbered in supplementary
table order, so `m_eggs_quad_c2` uses 532. The value sits at the call site, so
it travels with the model if the script is reordered.

## Environment

`session-info.md` records the version of R and of all 55 packages the pipeline
loads, from the machine that produced the committed results. There is no `renv`
lockfile, so that file is the only record of the environment. There is no
`.Rprofile` either, so R starts cleanly and the scripts can be run directly.

## Input data

`R_Input_files/` holds the raw records the pipeline reads: parentage
assignments, egg records, incubation bout durations, body mass, the
female–female and male–female social interaction datasets, copulations, and the
experimental treatment assignments. These are the only non-derived files in the
repository; everything in `R_Output_files/` is produced by the scripts above.
