# Deep Simulation-Based Inference for Inhomogeneous Bivariate Log-Gaussian Cox Processes


## Main Model

We consider a bivariate log-Gaussian Cox process with point types $p \in \lbrace 1,2 \rbrace$. The observed multitype point pattern is denoted by $X = \lbrace X_1, X_2 \rbrace,$ where $X_p$ contains the observed locations of points of type $p$.

For each type, the spatially varying intensity is defined by

$$\log \lbrace \Lambda_p(\mathbf{s}) \rbrace = \mathbf{Z}(\mathbf{s})^\top \boldsymbol{\beta}_p + Y(\mathbf{s}) + U_p(\mathbf{s}), \quad p \in \lbrace 1,2 \rbrace,$$

where:
* $\mathbf{s}$ denotes a spatial location;
* $\mathbf{Z}(\mathbf{s})$ is a vector of spatial covariates;
* $\boldsymbol{\beta}_p$ is the type-specific vector of first-order parameters;
* $Y(\mathbf{s})$ is a shared latent spatial field that induces dependence between the two point processes; and
* $U_p(\mathbf{s})$ is a type-specific latent spatial field that captures spatial variation unique to process $p$.
