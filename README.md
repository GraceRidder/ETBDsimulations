
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

# Run ETBD simulation with both population size dependent speciation and extinction 
res1 <- simETBD(
  t = 200,  # Number of time steps
  DIST = "NORM",  # Species abundance distribution
  JmaxV = 2000,   # Maximum number of individuals
  NegExpEx = TRUE,  # Population size dependent extinction
  exparm0 = -0.76, # Extinction parameter 0
  exparm1 = -0.56, # Extinction parameter 1
  consp = s.21,   # Constant probability of speciation (if size dependent speciation is disabled)
  ExpSp = TRUE,   # Size dependent speciation (at least one speciation mode must be TRUE)
  spparm1 =  -0.05, # Speciation parameter 1
  spparm0 = 0.95, # Speciation parameter 0
  conex = .21,   # Constant probability of extinction (if size dependent extinctioin is disabled)
  splitparm = scaled_samples[e, 2] # Heritability parameter
)

## Results

res1$tree          # Final Newick tree
res1$trees         # All Newick trees across time steps
res1$matrix_list   # Final population sizes
res1$matrix_lists  # Population sizes across all time steps
```
