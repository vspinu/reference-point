

### Code for the paper Baillon et. al. _Searching for the Reference Point_ Management Science 2018

Data is availabe in [`./data`](./data/README.md) sub-folder.

In order to replicate the results do:

- Clone the repository with

```sh
git clone --recursive https://github.com/vspinu/reference-point.git
```

- Install needed packages:

```R
install.packages(c("Rcpp", "devtools", "coda", "abind", 
                   "knitr", "ggplot2", "reshape2", "plyr", "knitr", "data.table", "xtable",
                   "methods", "xtable", "ggplot2", "purrr", "magrittr", "data.table",
                   "dplyr", "ggplot2", "lattice", "Hmisc"))
```

- Run the simulations and generate the pdf report:

```sh
./build.R 
```

Pass additional parameters to customize for alternative models

```sh
./build.R SEED I_nonLSKR I_LSKR W_nonLSKR W_LSKR util_type weight_type
```
For example

```sh
# sec 2 model 1               : ./build.R 17 50000 T T T T pow prelec1 (same as ./build.R)
# sec 3 model 1 with exp util : ./build.R 17 50000 T T T T exp_log prelec1
# sec 4 model 1 with prelec2  : ./build.R 17 50000 T T T T pow prelec2
# sec 5 model 1 with ibeta    : ./build.R 17 50000 T T T T pow ibeta
# sec 6 model 2               : ./build.R 17 50000 F T T T pow prelec1
# sec 7 model 3               : ./build.R 17 50000 T T T F pow prelec1
# sec 8 model 4               : ./build.R 17 50000 F T T F pow prelec1
```

!!!NOTE:!!!

  The simulation is using the PBM package which was written prior to the advent
  of STAN and is not under the active development. It's highly advisable to use
  STAN for similar task. Decision functions written in C++ (available in
  `decision-function.cpp`) might be useful.


