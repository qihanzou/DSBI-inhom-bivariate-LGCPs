# dsbi-inhom-bivariate-lgcp

Code and simulations for deep simulation-based inference in inhomogeneous bivariate log-Gaussian Cox process models, with an application to the gorilla dataset.

## Model

We consider a bivariate log-Gaussian Cox process with point types (p \in {1,2}). The observed multitype point pattern is denoted by $X = {X_1, X_2},$ where $X_p$ contains the observed locations of points of type $p$.

For each type, the spatially varying intensity is defined by

$$\log \{\Lambda_p(\mathbf{s})\} = \mathbf{Z}(\mathbf{s})^\top \boldsymbol{\beta}_p + Y(\mathbf{s}) + U_p(\mathbf{s}), \quad p \in {1,2},$$

where:
* $\mathbf{s}$ denotes a spatial location;
* $\mathbf{Z}(\mathbf{s})$ is a vector of spatial covariates;
* $\boldsymbol{\beta}_p$ is the type-specific vector of first-order parameters;
* $Y(\mathbf{s})$ is a shared latent spatial field that induces dependence between the two point processes; and
* $U_p(\mathbf{s})$ is a type-specific latent spatial field that captures spatial variation unique to process $p$.

This decomposition allows the model to represent both spatial dependence shared across the two point types and spatial variation specific to each type.
