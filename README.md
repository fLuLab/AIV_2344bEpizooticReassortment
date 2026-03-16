# The Winners Take It All? Evolutionary Success of H5Nx Reassortants in the 2020–2024 Panzootic

## **Table of Contents**
-   [Overview](#Overview)
-   [Downloading this repository](#Downloading-this-Repository-into-R)
-   [Dependencies](#Dependencies)
-   [Data Curation](#Data-Curation)
-   [Statistical Models](#Statistical-Models)
-   [Licence](#Licence)
<br /> 

## **Overview**

This repository contains code used for the data handling and statistical analysis of Clade 2.3.4.4b High Pathogenicity Avian Influenza Virus (HPAIV) circulating during the 2020-2024 Panzootic. Specifically, we analysed the emergence, persistence and drivers of unique reassortants worldwide, up to May 2024. Methods are described in detail at <https://www.biorxiv.org/content/10.1101/2025.07.19.665680v1>.
<br /> 
<br /> 

## **Downloading this Repository into R**

There are a few ways to download this repository and work with it in **R**.

### Option 1: Using `usethis::create_from_github()`

If you use the **usethis** package, you can pull the repository directly into RStudio:

``` r
# Install usethis if you don’t have it yet
install.packages("usethis")

# Replace 'username/repo' with the GitHub repo name
usethis::create_from_github("J-Baxter/AIV_2344bEpizooticReassortment")
```

### Option 2: Using `git clone` + RStudio

1.  Open a terminal (or Git Bash).\
2.  Run:

``` bash
git clone https://github.com/J-Baxter/AIV_2344bEpizooticReassortment.git
```

### Option 3: Download as ZIP

If you don’t want to use Git:

1.  On the GitHub repo page, click the green **Code** button → *Download ZIP*.

2.  Extract the folder to your computer.

3.  Open RStudio → *File* → *Open Project* → choose the `.Rproj` file inside the folder.
<br /> 
<br /> 

## **Dependencies** 

These statistical analyses were fitted in R version 4.5.1. The 'number of reassortants model' was fitted using Stan v2.36 via cmdstanr v0.9.0, and the remainder were fitted using BRMS v2.22.0. We summarised model outputs and caculated average marginal effects using tidybayes v3.0.7 and marginaleffects v0.25.1. All models have been tested on:

1.  Apple M4 Max 16-core CPU with 48 Gb RAM

2.  AMD Ryzen 9 7950X 16 Core CPU with 96 Gb RAM.

In each case, the runtime for any model was less than 10 minutes.
<br /> 
<br /> 

## **Data Curation**

Scripts used to assist in the curation of sequence and location data are included within the [data_curation](scripts/data_curation/) sub directory. Required helper functions are sourced from [funcs](scripts/funcs/).
<br /> 
<br /> 

## **Statistical Models** 

We fitted three statistical models to quantify patterns of reassortant emergence across continents and to understand the drivers of reassortant spatial diffusion. Each model has three associated scripts:
one fitting the model, one to run model evaluations and one for interpretation. All models are contained within the [statistical_models](scripts/statistical_models/) sub directory. Evaluation and
interpretation plots are produced in \*\_model_evaluation and \*model_interpretation scripts, however these may differ slightly from the final published plots (located in [scripts/figure_scripts](scripts/figure_scripts)).

The 'number of reassortants model' has additional scripts to describe the model in Stan (located in [scripts/statistical_models/stan_models](scripts/statistical_models/stan_models)) and for pre-processing.
<br /> 
<br /> 

### 1. Number of Reassortants Model

A mixture model comprised of three components, inspired by previously developed ecological models. For each year-month observation, $`i\in\{1,2,...,I\}`$, taken in continent, $`j\in\{\text{africa}, \text{asia}, \text{americas}, \text{europe}\}`$, let $`y_{ij}\in\mathbb{Z}_{\geq0}`$ be the observed number of reassortants. We assume $y_{ij}$ can be modelled as a mixture of three components: a detection model, an abundance model, and a zero-inflation model.\

First, we consider that only a proportion, $`p_{ij}\in(0,1)`$, of true (latent) reassortants, $`N_{ij}\in\mathbb{Z}_{\geq y_{ji}}`$, are ultimately observed:

```math
y_{ij}|N_{ij} \sim \mathrm{Binomial}(N_{ij},p_{ij})
```

Second, we model the true number of reassortants, $`N_{ij}\in\mathbb{Z}_{\geq y_{ji}}`$, as a discrete latent variable that follows a Poisson distribution:

```math
N_{ij} \sim \mathrm{Poisson}(\lambda_{ij})
```

where $`\lambda_{ij}`$ is the expected number of reassortants per observation on the log scale.

Third, we consider that ecological or epidemiological conditions may not always be conducive for reassortment/reassortant emergence. We assume this process is fundamentally distinct from a structural absence of reassortment (i.e situations where reassortment/reassortant emergence is feasible but does not occur). We model a Bernoulli zero-inflation component, $`z_{ij}\in\{0,1\}`$, parametrised by a continent-specific probability that a conditions are not permissive for reassortment/reassortant emergence, $`\theta_{i}`$:\

```math
z_{ij} \sim \mathrm{Bernoulli}(\theta_{\text{continent}[j]})
```
<br /> 
<br /> 

### 2. Reassortant Class Model

We estimated the probability that a novel reassortant, is assigned to one of the following classes: minor, moderate, major. We assumed that the probability a reassortant is assigned a given class follows a cumulative distribution, with classes increasing from minor to moderate to major. We modelled each class as the discretisation of a latent (unobserved) continuous variable, via threshold parameters which partition the distribution.
```math
 y_{i}=\begin{cases}
\text{minor}\qquad\qquad\text{if}\;\tau_{\text{minor}-1}<\tilde y_i\leq\tau_{\text{minor}},\\
\text{moderate}\qquad\text{if}\;\tau_{\text{minor}}<\tilde y_i\leq\tau_{\text{moderate}},\\
\text{major}\qquad\qquad\text{if}\;\tau_{\text{moderate}}<\tilde y_i<\tau_{\text{moderate}+1},
\end{cases}

```
<br /> 
<br /> 

### 3. Diffusion Model

We fitted a mixed model to predict the weighted diffusion coefficients calculated from our phylogeographic analysis for each novel reassortant. We restricted our analysis to reassortants with a clade size greater than 1, since we cannot confidently distinguish between reassortants that truly exist at a single locus and reassortants with limited (but non-zero) circulation and incomplete sampling. For all reassortants with non-zero, we assumed a gamma distribution parametrised such that,
```math
y_{i} \sim \mathrm{Gamma}(\kappa_{i},\theta_{i})
```
<br /> 
<br /> 

## **Licence** 
This code is shared under the **GPL-3.0 licence**.
