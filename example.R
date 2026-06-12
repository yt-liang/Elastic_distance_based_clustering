# Example usage of elastic_kmeans() on synthetic data:
# single-axis and multi-axis clustering with default/custom weights,
# with and without normalization, serial and parallel.

source("elastic_kmeans.R")

set.seed(123)

n_time <- 30
n_per_group <- 10
t <- seq_len(n_time)

make_curve <- function(shape, noise_sd) {
  base <- if (shape == "sine") sin(2 * pi * t / n_time) else cos(2 * pi * t / n_time)
  base + rnorm(n_time, sd = noise_sd)
}

# group A: sine-like curves, group B: cosine-like curves
group_a <- replicate(n_per_group, make_curve("sine", 0.1))
group_b <- replicate(n_per_group, make_curve("cosine", 0.1))
true_label <- rep(c(1, 2), each = n_per_group)

axis1 <- cbind(group_a, group_b)               # 30 x 20
axis2 <- cbind(2 * group_a, 0.5 * group_b)      # second axis, different scale/shape

cat("\n==== Single-axis ====\n")

cat("\n-- default weights, normalize = TRUE, parallel = FALSE --\n")
res1 <- elastic_kmeans(list(axis1), K = 2, seed = 1, verbose = FALSE)
print(table(res1$cluster_label, true_label))

cat("\n-- default weights, normalize = TRUE, parallel = TRUE --\n")
res2 <- elastic_kmeans(list(axis1), K = 2, seed = 1, parallel = TRUE, n_cores = 2, verbose = FALSE)
print(table(res2$cluster_label, true_label))

cat("\n==== Multi-axis (2 axes) ====\n")

cat("\n-- default (equal) weights, normalize = TRUE, parallel = FALSE --\n")
res3 <- elastic_kmeans(list(axis1, axis2), K = 2, seed = 1, verbose = FALSE)
print(table(res3$cluster_label, true_label))

cat("\n-- custom weight (emphasize axis 2), normalize = FALSE --\n")
res4 <- elastic_kmeans(list(axis1, axis2), K = 2, weight = c(0.1, 0.1, 0.4, 0.4),
                        normalize = FALSE, seed = 1, verbose = FALSE)
print(table(res4$cluster_label, true_label))

cat("\n-- default weights, normalize = TRUE, parallel = TRUE --\n")
res5 <- elastic_kmeans(list(axis1, axis2), K = 2, seed = 1, parallel = TRUE, n_cores = 2, verbose = FALSE)
print(table(res5$cluster_label, true_label))
