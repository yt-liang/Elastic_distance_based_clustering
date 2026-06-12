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

update_means <- function(data, cluster_label, K, dist_store) {
  axis_num <- length(data)
  n_time <- nrow(data[[1]])

  lapply(seq_len(axis_num), function(a) {
    out <- matrix(0, nrow = K, ncol = n_time)
    for (k in seq_len(K)) {
      idx <- which(cluster_label == k)
      if (length(idx) > 1) {
        out[k, ] <- rowMeans(data[[a]][, idx, drop = FALSE])
      } else if (length(idx) == 1) {
        out[k, ] <- data[[a]][, idx]
      } else {
        # empty cluster: re-seed with the point farthest from the other centers
        dist_other <- dist_store[-k, , drop = FALSE]
        far_id <- which(dist_other == max(dist_other), arr.ind = TRUE)[1, 2]
        out[k, ] <- data[[a]][, far_id]
      }
    }
    out
  })
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
    new_centers <- update_means(data, assigned$cluster_label, K, assigned$dist_store)

    mean_diff <- mean(abs(unlist(new_centers) - unlist(centers)))
    label_changed <- !identical(assigned$cluster_label, cluster_label)

    if (verbose) {
      message(sprintf("iter %d: mean center diff = %.5f, labels changed = %s",
                       iter, mean_diff, label_changed))
    }

    cluster_label <- assigned$cluster_label
    centers <- new_centers

    if (mean_diff <= criteria && !label_changed) {
      converged <- TRUE
      break
    }
  }

  list(cluster_label = cluster_label, centers = centers,
       iterations = iter, converged = converged)
}
