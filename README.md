# Elastic Distance-based Clustering

A general-purpose elastic-distance k-means clustering algorithm for
(multi-axis) functional / time-series data, implemented in R.

The algorithm assigns each curve to the cluster whose center curve has the
smallest elastic distance (phase + amplitude variation, via
[`fdasrvf::elastic.distance`](https://cran.r-project.org/package=fdasrvf)),
then iteratively updates cluster centers until convergence — similar in
spirit to k-means, but using elastic distance instead of Euclidean distance.

## Files

- `elastic_kmeans.R` — the core algorithm (`elastic_kmeans()` and its helper
  functions)
- `example.R` — runnable example on synthetic single-axis and multi-axis data

## Input format

`data` is a **list of matrices**, one matrix per axis/dimension. Each matrix
has shape `n_time_points x n_samples` (rows = time points of one curve,
columns = samples). For single-axis data, pass a list of length 1.

```r
# single-axis: 30 time points, 20 samples
data <- list(matrix_30x20)

# multi-axis (e.g. x/y/z accelerometer axes): same shape per axis
data <- list(axis_x, axis_y, axis_z)
```

## Usage

```r
source("elastic_kmeans.R")

result <- elastic_kmeans(
  data,
  K = 2,                # number of clusters
  weight = NULL,        # NULL = equal weight across all Dx/Dy components
  normalize = TRUE,     # global min-max normalize distances before weighting
  parallel = FALSE,     # set TRUE to use parallel::mclapply
  n_cores = 2,          # only used when parallel = TRUE
  criteria = 0.01,      # convergence threshold on mean center change
  max_iter = 500,
  seed = 1,
  verbose = TRUE
)

result$cluster_label   # integer vector, length = n_samples
result$centers         # list of K x n_time_points matrices, one per axis
result$iterations      # number of iterations run
result$converged       # TRUE/FALSE
```

## Arguments

| Argument | Description |
|---|---|
| `data` | List of `n_time_points x n_samples` matrices, one per axis |
| `K` | Number of clusters |
| `weight` | Numeric vector of length `2 * length(data)` to combine each axis's Dx/Dy elastic distances. Defaults to equal weights |
| `normalize` | Whether to globally min-max normalize each distance component across clusters before applying `weight` |
| `parallel` | Whether to parallelize distance computation with `parallel::mclapply` |
| `n_cores` | Number of cores to use when `parallel = TRUE` |
| `criteria` | Convergence threshold on the mean absolute change in cluster centers |
| `max_iter` | Maximum number of iterations |
| `seed` | Optional random seed for initial center selection |
| `verbose` | Whether to print progress per iteration |

## Output

A list with:
- `cluster_label`: cluster assignment for each sample
- `centers`: list (one per axis) of `K x n_time_points` matrices — the final
  cluster center curves
- `iterations`: number of iterations actually run
- `converged`: whether the algorithm converged before `max_iter`

## Dependencies

- [`fdasrvf`](https://cran.r-project.org/package=fdasrvf) — elastic distance computation
- `parallel` (base R) — optional parallel processing

## Example

See `example.R` for a runnable example with synthetic single-axis and
multi-axis data, demonstrating default/custom weights, with/without
normalization, and serial/parallel execution.

```bash
Rscript example.R
```

## Reference
Liang, YT., Wang, C. Motif clustering and digital biomarker extraction for free-living physical activity analysis. BioData Mining 18, 8 (2025). https://doi.org/10.1186/s13040-025-00424-1
