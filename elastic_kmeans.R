# Elastic-distance k-means for (multi-axis) functional data
#
# `data` is a list of matrices, one per axis/dimension, each of shape
# (n_time_points x n_samples). A single axis is just a list of length 1.

elastic_dist_pair <- function(x, ref, t) {
  d <- fdasrvf::elastic.distance(as.numeric(x), as.numeric(ref), t)
  c(Dx = d$Dx, Dy = d$Dy)
}

compute_distances <- function(data, ref_mat_list, t, parallel, n_cores) {
  axis_num <- length(data)
  K <- nrow(ref_mat_list[[1]])
  n_samples <- ncol(data[[1]])

  apply_fun <- if (parallel) {
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = n_cores)
  } else {
    lapply
  }

  lapply(seq_len(K), function(k) {
    per_axis <- lapply(seq_len(axis_num), function(a) {
      ref <- ref_mat_list[[a]][k, ]
      res <- apply_fun(seq_len(n_samples), function(i) {
        elastic_dist_pair(data[[a]][, i], ref, t)
      })
      do.call(rbind, res)
    })
    do.call(cbind, per_axis)
  })
}

normalize_distances <- function(dist_list) {
  all_dist <- do.call(rbind, dist_list)
  col_min <- apply(all_dist, 2, min)
  col_max <- apply(all_dist, 2, max)
  col_range <- col_max - col_min
  col_range[col_range == 0] <- 1
  lapply(dist_list, function(m) {
    sweep(sweep(m, 2, col_min, "-"), 2, col_range, "/")
  })
}

assign_clusters <- function(data, ref_mat_list, weight, t, normalize, parallel, n_cores) {
  dist_list <- compute_distances(data, ref_mat_list, t, parallel, n_cores)
  if (normalize) dist_list <- normalize_distances(dist_list)

  # dist_store: K x n_samples matrix of weighted combined distances
  dist_store <- do.call(rbind, lapply(dist_list, function(m) as.numeric(m %*% weight)))
  cluster_label <- apply(dist_store, 2, which.min)

  list(cluster_label = cluster_label, dist_store = dist_store)
}

#' Relocate points into any empty clusters (sklearn-style)
#'
#' For each empty cluster, steal the sample that currently fits its own
#' assigned cluster worst (i.e. has the largest distance to its own cluster
#' center), and reassign it to the empty cluster. If there are multiple empty
#' clusters, the top-N worst-fit points (across all samples) are chosen at
#' once so no two empty clusters compete for the same point.
#'
#' @return The (possibly modified) `cluster_label` vector.
relocate_empty_clusters <- function(cluster_label, K, dist_store) {
  n_samples <- length(cluster_label)
  counts <- tabulate(cluster_label, nbins = K)
  empty_clusters <- which(counts == 0)

  if (length(empty_clusters) == 0) return(cluster_label)

  # distance of each sample to the center of its own assigned cluster
  own_dist <- dist_store[cbind(cluster_label, seq_len(n_samples))]

  ord <- order(own_dist, decreasing = TRUE)
  n_needed <- length(empty_clusters)

  if (n_needed > n_samples) {
    stop("More empty clusters than samples available to relocate; ",
         "reduce K relative to the number of samples.")
  }

  chosen <- ord[seq_len(n_needed)]
  cluster_label[chosen] <- empty_clusters
  cluster_label
}

update_means <- function(data, cluster_label, K, dist_store) {
  axis_num <- length(data)
  n_time <- nrow(data[[1]])

  cluster_label <- relocate_empty_clusters(cluster_label, K, dist_store)

  centers <- lapply(seq_len(axis_num), function(a) {
    out <- matrix(0, nrow = K, ncol = n_time)
    for (k in seq_len(K)) {
      idx <- which(cluster_label == k)
      if (length(idx) == 0) {
        stop(sprintf(
          "Cluster %d is still empty after relocation; K may be too large for this data.",
          k))
      }
      out[k, ] <- if (length(idx) > 1) {
        rowMeans(data[[a]][, idx, drop = FALSE])
      } else {
        data[[a]][, idx]
      }
    }
    out
  })

  list(centers = centers, cluster_label = cluster_label)
}

#' Elastic-distance k-means clustering for multi-axis functional data
#'
#' @param data List of (n_time_points x n_samples) matrices, one per axis.
#' @param K Number of clusters.
#' @param weight Numeric vector of length 2 * length(data) used to combine the
#'   per-axis Dx/Dy elastic distances. Defaults to equal weights.
#' @param normalize Whether to globally min-max normalize each distance
#'   component across clusters before applying `weight`.
#' @param parallel Whether to parallelize distance computations with
#'   `parallel::mclapply`.
#' @param n_cores Number of cores to use when `parallel = TRUE`.
#' @param criteria Convergence threshold on the mean absolute change in
#'   cluster centers.
#' @param max_iter Maximum number of iterations.
#' @param seed Optional random seed for initial center selection.
#' @param verbose Whether to print progress per iteration.
#'
#' @return A list with `cluster_label`, `centers` (list of K x n_time_points
#'   matrices, one per axis), `iterations`, and `converged`.
elastic_kmeans <- function(data, K, weight = NULL, normalize = TRUE,
                            parallel = FALSE,
                            n_cores = max(1, parallel::detectCores() - 1),
                            criteria = 0.01, max_iter = 500,
                            seed = NULL, verbose = TRUE) {
  if (!is.list(data)) stop("`data` must be a list of matrices (one per axis).")

  axis_num <- length(data)
  n_time <- nrow(data[[1]])
  n_samples <- ncol(data[[1]])
  for (a in seq_along(data)) {
    if (nrow(data[[a]]) != n_time || ncol(data[[a]]) != n_samples) {
      stop("All axes in `data` must have the same dimensions.")
    }
  }

  if (is.null(weight)) {
    weight <- rep(1 / (axis_num * 2), axis_num * 2)
  } else if (length(weight) != axis_num * 2) {
    stop("`weight` must have length 2 * length(data).")
  }

  if (!is.null(seed)) set.seed(seed)
  t_index <- seq_len(n_time)

  sample_id <- sample(n_samples, K)
  centers <- lapply(data, function(m) t(m[, sample_id, drop = FALSE]))

  cluster_label <- rep(NA_integer_, n_samples)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    assigned <- assign_clusters(data, centers, weight, t_index, normalize, parallel, n_cores)
    updated <- update_means(data, assigned$cluster_label, K, assigned$dist_store)

    mean_diff <- mean(abs(unlist(updated$centers) - unlist(centers)))
    label_changed <- !identical(updated$cluster_label, cluster_label)

    if (verbose) {
      message(sprintf("iter %d: mean center diff = %.5f, labels changed = %s",
                       iter, mean_diff, label_changed))
    }

    cluster_label <- updated$cluster_label
    centers <- updated$centers

    if (mean_diff <= criteria && !label_changed) {
      converged <- TRUE
      break
    }
  }

  list(cluster_label = cluster_label, centers = centers,
       iterations = iter, converged = converged)
}
