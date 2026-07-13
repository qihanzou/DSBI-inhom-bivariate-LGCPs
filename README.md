# dsbi-inhom-bivariate-lgcp
Code and simulations for deep simulation-based inference in inhomogeneous bivariate log-Gaussian Cox process models.

We consider a bivariate log-Gaussian Cox process for types $p = 1,2$, we have two types of point patterns and $X=\{X_1, X_2\}$ and the following intensity structure, 
$$\log\{\Lambda_p(\bm s)\} = \bm Z(\bm s)\bm \beta_p + Y(\bm s) + U_p(\bm s),$$
where $Y(\bm s)$ is the shared field that affects all types of processes and $U_p(\bm s)$ individual field, which only affects individually the type $p$ process.
