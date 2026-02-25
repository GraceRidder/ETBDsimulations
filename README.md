
## ETBD Simulaions

## Installation

You can install the development version of **ETBDsimulations** from GitHub using the `devtools` package.

### 1. Install `devtools` (if not already installed)

```r
install.packages("devtools")
```

### 2. Install ETBDsimulations from GitHub

```r
devtools::install_github("GraceRidder/ETBDsimulations")
```

### 3. Load the package

```r
library(ETBDsim)
```

## Example Usage

Below is an example of how to run the `simETBD()` function in R and inspect the outputs.

```r

# Run ETBD simulation
res1 <- simETBD(
  t = 200,                        # Number of time steps
  DIST = "NORM",                  # Species abundance distribution
  JmaxV = scaled_samples[e, 1],   # Maximum number of individuals
  NegExpEx = TRUE,                # Population size dependent extinction
  exparm0 = scaled_samples[e, 3], # Extinction parameter 0
  exparm1 = scaled_samples[e, 4], # Extinction parameter 1
  consp = scaled_samples[e, 5],   # Constant probability of speciation (if enabled)
  ExpSp = TRUE,                   # Size dependent speciation (at least one speciation mode must be TRUE)
  spparm1 = scaled_samples[e, 6], # Speciation parameter 1
  spparm0 = scaled_samples[e, 7], # Speciation parameter 0
  conex = scaled_samples[e, 5],   # Constant probability of extinction (if enabled)
  splitparm = scaled_samples[e, 2] # Heritability parameter
)

## Results

res1$tree          # Final Newick tree
res1$trees         # All Newick trees across time steps
res1$matrix_list   # Final population sizes
res1$matrix_lists  # Population sizes across all time steps
```
