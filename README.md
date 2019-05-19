

Code for the paper Baillon et. al. _Searching for the Reference Point_ Management Science 2018

Data is availabe in [`./data`](./data/) sub-folder.

### Replicate the results

- Clone the repository with

```sh
git clone --recursive https://github.com/vspinu/reference-point.git
```

- Install needed packages:

```R
install.packages(c("Rcpp", "devtools", "coda", "abind", 
                   "knitr", "ggplot2", "plyr", "dplyr", "data.table", "xtable",
                   "purrr", "magrittr", "data.table", "Hmisc"))
```

- Run the simulations and generate the pdf report:

```sh
./build.R 

## Pass additional parameters for alternative models
#
#    ./build.R SEED I_nonLSKR I_LSKR W_nonLSKR W_LSKR util_type weight_type
#
# where
#
#             SEED: random seed
# I_nonLSKR/I_LSKR: TRUE of FALSE for I term (consumption utility) for non-LS&KR
#                   or LS&KR models respectively
# W_nonLSKR/W_LSKR: TRUE of FALSE for weighting term for non-LS&KR or LS&KR
#                   models respectively
#        util_type: utility type - "pow" or "exp_log"
#      weight_type: 1 parameter weighting - "prelec1", "pow1", "linlog1", "TK"
#                   2 parameter weighting - "prelec2", "pow2", "linlog2", "ibeta"

## For example

# sec 2 model 1               : ./build.R 17 50000 T T T T pow prelec1 (same as ./build.R)
# sec 3 model 1 with exp util : ./build.R 17 50000 T T T T exp_log prelec1
# sec 4 model 1 with prelec2  : ./build.R 17 50000 T T T T pow prelec2
# sec 5 model 1 with ibeta    : ./build.R 17 50000 T T T T pow ibeta
# sec 6 model 2               : ./build.R 17 50000 F T T T pow prelec1
# sec 7 model 3               : ./build.R 17 50000 T T T F pow prelec1
# sec 8 model 4               : ./build.R 17 50000 F T T F pow prelec1
```

You must have `pdflatex` on your machine to be able to produce the pdf report.

### Use on your own data

If you can put your data in the required format set the environment variable
`RP-DATA-CSV` to your data file and run `./build.R` as above. 

```sh
export RP-DATA-CSV="/path/to/your-data.csv"
./build.R
```

Data requirements are:

- Must contain columns IdSubject, IdQuestion and Preference. IdSubject and
  IdQuestion can be of arbitrary type but they will be internally converted to
  integer indexes ranging from 1 to max number of subjects or
  questions. `Preference` must take integer values 1 (A prospect) and 2 (B prospect).
- Outcomes and probabilities for prospect A must be named x1a, x2a ... and p1a,
  p2a,... . Similarly, for prospect B they must be x1b, x2b ... and p1b,
  p2b,... .
- Prospects of various length are allowed. Structural missing values must be
  blanks (or NA), not 0!
- Outcomes of each prospects must be in increasing order. That is, x1a < x2a <
  ..., x1b < x2b < ....

Besides the pdf report, an RDS archive containing MCMC simulations will be
stored in `mcmc.rds` file. You can use it for your own inferences.

```
> mc <- readRDS("mcmc.rds")
> str(mc)
List of 5
 $ mean_pop: 'mcmc' num [1:6500, 1:4] 0.473 0.496 0.518 0.465 0.494 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:4] "mu[alpha]" "mu[gamma]" "mu[lambda]" "mu[xi]"
  ..- attr(*, "mcpar")= num [1:3] 6501 13000 1
 $ sd_pop  : 'mcmc' num [1:6500, 1:4] 0.208 0.232 0.281 0.219 0.234 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:4] "sigma[alpha]" "sigma[gamma]" "sigma[lambda]" "sigma[xi]"
  ..- attr(*, "mcpar")= num [1:3] 6501 13000 1
 $ ind     : 'mcmc' num [1:6500, 1:556] 0.228 0.228 0.228 0.228 0.228 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:556] "l1" "l2" "l3" "l4" ...
  ..- attr(*, "mcpar")= num [1:3] 6501 13000 1
 $ RP_pop  : 'mcmc' num [1:6500, 1:6] 0.306 0.275 0.327 0.232 0.297 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:6] "Status Quo" "MaxMin" "MinMax" "X at Max P" ...
  ..- attr(*, "mcpar")= num [1:3] 6501 13000 1
 $ RP_ind  : 'mcmc' num [1:6500, 1:139] 6 6 6 6 6 6 6 6 6 6 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:139] "RP1" "RP2" "RP3" "RP4" ...
  ..- attr(*, "mcpar")= num [1:3] 6501 13000 1

```


!!!NOTE!!!

    The simulation uses the PBM package which was written prior to the advent of
    STAN and is no longer developed. It's highly advisable to use STAN for
    similar tasks in the future. Decision functions written in C++ (available in
    `decision-function.cpp`) might be useful in that regard.


