## Hurdle Truncated Log-Normal Distribution for Zero-Inflated Heavy-Tailed Data

**htln** provides tools for working with a **hurdle truncated log-normal distribution**,  
a simple yet effective model for **zero-inflated and heavy-tailed data**, particularly
suited for next–generation sequencing (NGS) counts after applying the transformation
`log(x + 1)`.

The package implements:

- **`mle_htln()`** – maximum likelihood parameter estimation  
- **`dhtln()`** – density  
- **`qhtln()`** – quantile function  
- **`rhtln()`** – random generation  

All functions operate on the log-transformed scale  
\[
Z = \\log(Y + 1),
\]
where the continuous component is modeled as a **truncated Normal distribution**.

---

## 📦 Installation

```r
# Install the development version from GitHub
# install.packages("remotes")
remotes::install_github("Fuschi/htln")
```

---

## 🔍 Motivation

NGS count data and biological measurements often exhibit:

- many zeros 
- heavy-tailed positive values 
- strong skewness 

A **hurdle truncated log-normal model** tries to capture these properties by combining:

1. **A hurdle component**

$$
P(Y = 0) = \phi
$$

2. **A continuous component** on the log-scale

$$
Z = \log(Y + 1) \sim \mathrm{TruncNorm}(\mu, \sigma^2; [0, b])
$$

where the upper bound $b$ is automatically chosen to stabilise numerical optimisation.

