# BayesCountFDR: Bayesian local false discovery rate for sparse count data

Code for

> Gauran, I. I. M., Park, J., Rattsev, I., Peterson, T. A., Kann, M. G., and
> Park, D. (2022). Bayesian local false discovery rate for sparse count data
> with application to the discovery of hotspots in protein domains. *Annals of
> Applied Statistics* 16(3):1459-1475. doi:10.1214/21-AOAS1551.

The method finds the unusually large counts in a sparse vector of counts, for
example the positions of a protein domain that carry more mutations than the
rest (hotspots). The counts follow a two-group mixture: a zero-inflated
generalized Poisson (ZIGP) null and an alternative that starts above an
unknown cutoff. `bayes_local_fdr()` in `R/bayes_fdr.R` samples the posterior
by Gibbs sampling and returns the local false discovery rate of every count;
the alternative can be nonparametric (NP), semiparametric (SP), Poisson (P)
or generalized Poisson (GP). `R/eb_fdr.R` holds the empirical Bayes methods
the paper compares it with. Each table and figure of the paper that comes
from code has one script in `scripts/`, and `nix flake check` builds all of
them and compares them with the paper. Figure 2 of the paper is a picture of
a protein structure and has no code.

The code exists twice, in R (`R/`, `scripts/`) and in Python (`python/`).
The Python version is a second implementation that gives the same numbers as
R: the same random draws, the same Gibbs chains, every csv cell equal and
every figure equal mark by mark (`checks/PYTHON.md`). R is the reference.
A third version, `rcpp/`, is the R code with its sampler written in C++ (see
Rcpp version); it gives the same numbers as R too (`checks/RCPP.md`).

```r
source("R/setup.R")
x <- read_domain("cd00083", PATHS$domains)$counts   # mutation counts of 93 positions
set.seed(1234)
fit <- bayes_local_fdr(x, model = "NP", gamma = 1)
data.frame(count = fit$values, positions = fit$frequency, fdr = round(fit$fdr, 3))
eb_rejections(x, level = 0.05)                      # the empirical Bayes count
```

## Layout

| Path | Content |
|---|---|
| `R/setup.R` | Loader sourced by every script: settings of the paper, paths, argument parser |
| `R/distributions.R` | ZIGP density and sampler, shifted Poisson and generalized Poisson densities, Dirichlet and multinomial samplers |
| `R/bayes_fdr.R` | `bayes_local_fdr()`, the Gibbs sampler with its four alternatives, and `SAMPLER_VERSIONS`, the versions of the sampler used by the runs of the paper |
| `R/eb_fdr.R` | Empirical Bayes: EM fit of the ZIGP null, choice of the cutoff, the two-stage and one-stage procedures and Storey's procedure |
| `R/simulate.R` | Simulated data (`rmixture()`), the settings of the simulation runs, and the simulation study run in pieces and summarized |
| `R/real_data.R` | Reader of the protein domain files and the analysis of a domain (Table 3) |
| `R/plots.R` | Bar plots, line plots and the shaded table |
| `scripts/` | One script per table or figure |
| `python/localfdr/` | The Python version of `R/`, one module per R file, with `rng.py` (R's random number generators), `rmath.py` and `rbase.py` (R's arithmetic, and BLAS and LAPACK through R's own library) |
| `python/scripts/` | The Python version of each script, same options and output files |
| `rcpp/src/` | The C++ version of `R/distributions.R`, of the Gibbs sampler of `R/bayes_fdr.R` and of `rmixture()`, one file per R file, with `distributions_core.h` (the code the files share) and `rbase.h` (R's arithmetic, and BLAS and LAPACK called as R calls them) |
| `rcpp/R/load.R` | Loader of the Rcpp version: compiles `rcpp/src/` and puts each C++ function in place of the R function of the same name |
| `rcpp/scripts/` | The Rcpp version of each script: the R script run on the C++ functions, same options and output files |
| `tests/`, `python/tests/`, `rcpp/tests/` | Unit tests of every R function (testthat), every Python function (pytest) and every C++ function (testthat, against its R function) |
| `archive/` | Summaries of the original simulation runs that the scripts read (`simulation_summaries.csv`), the tables behind the published Figures S1 and S3 (`figS1_plot_input.csv`, `figS3_full_format.csv`), the published figure files (`figures/`, named as the figures the scripts draw), and `reduce_archive.R`, which made all of these from the run folders |
| `data/domains/` | Mutation counts of the six protein domains |
| `checks/` | The numbers printed in the paper (`paper_values.csv`), the comparison scripts and their result (`RESULTS.md`), `check_equivalence.R`, which compares the functions of `R/` with the original code, and `versions.sh`, which prints the software versions; `check_python.sh`, which runs every R and Python script and compares their outputs (`PYTHON.md`); `check_rcpp.sh`, which does the same for the R and Rcpp scripts (`RCPP.md`) |
| `.github/workflows/ci.yml` | Runs `nix flake check` on every push and pull request |
| `flake.nix`, `flake.lock` | Pinned environment (R and Python with their packages, Rcpp and a C++ compiler), one command per output in each version, the unit tests and the checks |

## Running

With [Nix](https://nixos.org) (flakes enabled), from the repository root:

```sh
nix run .#table3        # one output: table1 to table3, fig1, figS1 to figS3
nix run .#all           # every output
nix flake check         # build everything, compare with the paper and R with Python and Rcpp, run the tests
nix run .#py-table3     # the same output from the Python code (py-all for every one)
nix run .#rcpp-table3   # the same output from the Rcpp code (rcpp-all for every one)
nix run .#test          # unit tests of the R, the Rcpp and the Python code
nix run .#versions      # print the software versions
nix run .#help          # every command (also: nix run . -- --help)
```

Tables are written to `results/` (csv) and figures to `figures/` (pdf), each
under the name of its table or figure in the paper (see Outputs).
`nix flake check` compares every cell of Tables 1 to 3 with the paper and
every figure with the published one; `checks/RESULTS.md` holds the comparison.

Table 3 and Figures 1 and S2 are always computed, from the data in
`data/domains/` or from simulated data. Tables 1 and 2 and Figures S1 and S3
are simulation studies of 1000 replications per setting; by default their
scripts read the summaries of the original runs in `archive/`. To run a
study here instead, add `--rerun`:

```sh
nix run .#table1 -- --rerun --runs=1:100 --cores=8
```

`nix develop` opens a shell with the same R, where the scripts can be run
directly: `Rscript scripts/table3_protein_domains.R [--cores=4]`. Each script
lists its options in its header. The settings are those of the paper: local
fdr level 0.05, 1000 replications, 5000 stored iterations after a burn-in of
1000 in the simulations, 10000 after 1000 for the protein domains.

## Outputs

Each table and figure of the paper, the script that builds it, what it reads,
the files it writes (figures under `figures/`) and how `nix flake check`
checks them.

| Paper | Script | Input | Output | Check |
|---|---|---|---|---|
| Table 1 | `scripts/table1_simulation.R` | `archive/simulation_summaries.csv` | `results/Table1_simulation.csv` | every cell against the paper (`checks/paper_values.csv`) |
| Table 2 | `scripts/table2_prior_sensitivity.R` | `archive/simulation_summaries.csv` | `results/Table2_prior-sensitivity.csv` | every cell against the paper (`checks/paper_values.csv`) |
| Table 3 | `scripts/table3_protein_domains.R` | `data/domains/<domain>.txt` | `results/Table3_protein-domains.csv` | every cell against the paper (`checks/paper_values.csv`) |
| Figure 1 | `scripts/fig1_domain_barplots.R` | `data/domains/<domain>.txt` | `Fig1A_mutation-counts_cd00054.pdf`<br>`Fig1B_mutation-counts_cd00083.pdf`<br>`Fig1C_mutation-counts_cd00590.pdf`<br>`Fig1D_mutation-counts_smart00179.pdf`<br>`Fig1E_mutation-counts_smart00209.pdf`<br>`Fig1F_mutation-counts_cd00063.pdf`<br>`results/Fig1_frequencies.csv` | every bar of each PDF against the published file of the same name in `archive/figures/` |
| Figure S1 | `scripts/figS1_fdr_tpr_well_separated.R` | `archive/simulation_summaries.csv` | `FigS1A_FDR_well-separated_setting16.pdf`<br>`FigS1B_TPR_well-separated_setting16.pdf`<br>`FigS1C_FDR_well-separated_setting18.pdf`<br>`FigS1D_TPR_well-separated_setting18.pdf`<br>`results/FigS1_values.csv` | the six lines of each PDF against the published file of the same name; the plotted values against `archive/figS1_plot_input.csv` |
| Figure S2 | `scripts/figS2_histograms.R` | none (simulated data, seed 1439) | `FigS2A_simulated-counts_N50_lambda2.pdf` to `FigS2I_simulated-counts_N200_lambda4.pdf` (N = 50, 100, 200 for lambda = 2, 3, 4)<br>`results/FigS2_frequencies.csv` | every bar of each PDF against the published file of the same name |
| Figure S3 | `scripts/figS3_fdr_table.R` | `archive/simulation_summaries.csv` | `results/FigS3_fdr_table.csv`<br>`FigS3_fdr_table.pdf` | every cell against `archive/figS3_full_format.csv`, the sheet the published table shows |

With `--rerun`, Tables 1 and 2 and Figures S1 and S3 simulate their data and
store one file per setting, method and replication under
`results/replications/`, so a study can be run in pieces and continued:
`--runs=1:100` runs replications 1 to 100, and the table is written once
replications 1 to 1000 are all there. Replications shared by two outputs
(Table 1 and Table 2 with the uniform prior, Table 1 and Figure S3) are run
once. The comparisons with the paper are in `checks/RESULTS.md`.

Commands and run times:

| Paper | Command | Default | With `--rerun` |
|---|---|---|---|
| Table 1 | `nix run .#table1` | 1 s | about 3.5 h, 8 cores (estimate) |
| Table 2 | `nix run .#table2` | 1 s | about 5 h, 8 cores (estimate) |
| Table 3 | `nix run .#table3` | 1 min, 2 cores (always computed from the data) | same |
| Figure 1 | `nix run .#fig1` | 1 s (always computed from the data) | same |
| Figure S1 | `nix run .#figS1` | 1 s | about 9 h, 8 cores (estimate) |
| Figure S2 | `nix run .#figS2` | 1 s (always computed) | same |
| Figure S3 | `nix run .#figS3` | 1 s | about 23 h, 8 cores (estimate) |

Times are from an 8 core Apple M2 laptop; estimates are the processor time of
one replication of every setting of the output, times 1000, divided by 8.

Tables 1 and 2 and Figures S1 and S3 are built from the summaries of the
original simulation runs of February 2021 (`archive/`); `--rerun` recomputes
them with the seeds, settings and sampler
version of each run (`SIMULATION_RUNS` in `R/simulate.R`). The protein domain
analysis of Table 3 uses the sampler version of 2018 (`SAMPLER_VERSIONS` in
`R/bayes_fdr.R`). `checks/check_equivalence.R` runs the original functions
and the functions of `R/` on the same cases and compares their results. On a
different machine or R build, recomputed results can differ.

## Rcpp version

`rcpp/` is the R code with the parts that take the time written in C++
(Rcpp). `rcpp/src/` holds the functions of `R/distributions.R` (the ZIGP
density and draws, the shifted densities, the Dirichlet and multinomial
samplers), the whole Gibbs sampler of `R/bayes_fdr.R` (`bayes_local_fdr()`
and every function it calls, for the four alternatives and the three versions
of the sampler) and `rmixture()` of `R/simulate.R`. Each C++ function has the
name and the arguments of its R function. `rcpp/R/load.R` compiles them and
puts them in place of the R functions, so everything else runs as before and
calls C++: the scripts in `rcpp/scripts/` run the scripts of `scripts/` this
way, with the same options and output files.

```r
source("R/setup.R")
source("rcpp/R/load.R")     # compiles rcpp/src/ on the first call
x <- read_domain("cd00083", PATHS$domains)$counts
set.seed(1234)
fit <- bayes_local_fdr(x, model = "NP", gamma = 1)   # the chain of the R code, run in C++
```

What stays in R: the empirical Bayes methods (`R/eb_fdr.R`), the settings
(`sampler_settings()`, `SAMPLER_VERSIONS`), the reader of the domain files,
the driver of the simulation studies, the tables and the plots. Inside a
chain two things are still asked of R: the covariance matrix of the chain
behind the adapted proposal (`var()`, once every 100 iterations) and the two
small functions stored in a sampler version.

The outputs equal those of R in every value (`checks/RCPP.md`): the C++ code
draws from R's random number generator in the order of the R code, calls R's
own mathematical functions and the BLAS and LAPACK library R is linked to,
and adds its sums as R does. `rcpp/tests/` holds the unit tests: every C++
function against its R function under fixed seeds, with identical values and
the same state of the random number generator afterwards, up to whole seeded
chains.

Run times, 2 cores: Table 3 takes 54 s in R and 5 s with Rcpp; two
replications of every setting take 100 s and 10 s for Table 1 and 119 s and
7 s for Table 2. Compiling the C++ files takes about 8 s, once; the
compiled files are kept in the temporary folder of the machine (or in the
folder named by `LOCALFDR_RCPP_CACHE`).

## Tests

Every function in `R/` has a unit test (`tests/testthat/test-<file>.R`, 68
functions) and every function of `python/localfdr/` has one
(`python/tests/test_<module>.py`), on small inputs with a known answer: a
density that sums to one, a sampler's mean over many draws, a short chain that
gives the same draws twice under one seed. The Python tests marked `against_r`
also compare random streams, densities and whole Gibbs chains with R, draw by
draw. The tests of the Rcpp version (`rcpp/tests/testthat/test-rcpp.R`)
compare each of the 22 C++ functions with its R function. `nix run .#test`
runs the three suites; `nix flake check` runs them with the output checks,
and GitHub runs `nix flake check` on every push.

Each simulated data set and each chain is seeded where it is drawn (seed plus
replication number in the simulation studies, 1234 for each domain and model
in Table 3), so a script gives the same output run alone or after another.

## Data

- **Protein domain mutation counts** (`data/domains/<domain>.txt`): six
  domains, cd00054 (90 positions), cd00083 (93), cd00590 (175), smart00179
  (182), smart00209 (188) and cd00063 (193), from the Domain Mapping of
  Disease Mutations database (DMDM). The first line of a file holds the
  number of positions N; every further line is `position:count:mutations`
  for positions 0 to N, where position 0 holds no mutation. `read_domain()`
  returns the counts of positions 1 to N. Used by Table 3 and Figure 1.
- **Archived runs** (`archive/simulation_summaries.csv`): for each setting of
  the original simulation runs and each method, the mean and standard
  deviation over 1000 replications of the number of rejections R, the FDR
  and the TPR, as the runs wrote them. `archive/reduce_archive.R` made this
  file, `figS1_plot_input.csv` and `figS3_full_format.csv` from the output
  folders of the runs.
- **Published figures** (`archive/figures/`): the figure files of the paper
  that come from code, under the names the scripts give their figures; the
  published Figure S3 is a picture of a table (`FigS3_fdr_table.png`).

## Software versions

`flake.lock` pins the nixpkgs revision, and with it the version of R, Python
and every package (there is no `renv.lock`). `nix run .#versions` (or
`bash checks/versions.sh` inside `nix develop`) prints them; `nix flake check`
compares this section with its output.

```text
nixpkgs  NixOS/nixpkgs f3aa7d3d72c2bc66e0f30f4262d7ec2e33c9eab2

Python   3.14.7
  numpy        2.5.2
  pandas       3.0.4
  matplotlib   3.11.1
  numba        0.67.0
  pytest       9.1.1

R        4.6.1
  boot         1.3-32
  brio         1.1.5
  callr        3.8.0
  class        7.3-23
  cli          3.6.6
  cluster      2.1.8.2
  codetools    0.2-20
  cpp11        0.5.5
  crayon       1.5.3
  desc         1.4.3
  diffobj      0.3.6
  evaluate     1.0.5
  farver       2.1.2
  foreign      0.8-91
  fs           2.1.0
  ggplot2      4.0.3
  glue         1.8.1
  gtable       0.3.6
  isoband      0.3.0
  jsonlite     2.0.0
  KernSmooth   2.23-26
  labeling     0.4.3
  lattice      0.22-9
  lifecycle    1.0.5
  magrittr     2.0.5
  MASS         7.3-65
  Matrix       1.7-5
  mgcv         1.9-4
  nlme         3.1-169
  nnet         7.3-20
  openxlsx     4.2.8.1
  otel         0.2.0
  pkgbuild     1.4.8
  pkgload      1.5.3
  praise       1.0.0
  processx     3.9.0
  ps           1.9.3
  R6           2.6.1
  RColorBrewer 1.1-3
  Rcpp         1.1.2
  rlang        1.3.0
  rpart        4.1.27
  rprojroot    2.1.1
  S7           0.2.2
  scales       1.4.0
  spatial      7.3-18
  stringi      1.8.7
  survival     3.8-6
  testthat     3.3.2
  vctrs        0.7.3
  viridisLite  0.4.3
  waldo        0.6.2
  withr        3.0.3
  zip          3.0.0

OpenBLAS 0.3.33
```
