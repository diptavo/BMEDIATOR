jg02_default_variances <- function(role) {
    data.frame(
        v_x = unname(c(A = 0.0400, B = 0.0004, C = 0.0400)[role]),
        v_m = unname(c(A = 0.0025, B = 0.0400, C = 0.0400)[role]),
        v_y = 0.0025
    )
}

jg02_block_ld <- function(block, rho) {
    n <- length(block)
    result <- diag(n)
    for (current in unique(block)) {
        index <- which(block == current)
        distance <- abs(outer(seq_along(index), seq_along(index), "-"))
        result[index, index] <- rho^distance
    }
    result
}

jg02_simulate <- function(seed, blocks = 12L, variants_per_block = 3L,
                          ld_rho = 0.40, a = 0, b = 0, c_path = 0.20,
                          lambda = 0, q = 0, eta = 0,
                          sampling_rho = c(xm = 0, xy = 0, my = 0),
                          reported_sampling_rho = sampling_rho,
                          se = 0.02, true_variance_scale = 1,
                          reported_variance_scale = true_variance_scale,
                          orientation_accuracy = 1) {
    set.seed(seed)
    n <- blocks * variants_per_block
    variant <- sprintf("v%04d", seq_len(n))
    ld_block <- sprintf("block%03d", rep(seq_len(blocks), each = variants_per_block))
    role <- rep(c("A", "B", "C"), length.out = n)
    variance <- jg02_default_variances(role)
    variance[c("v_x", "v_m", "v_y")] <-
        variance[c("v_x", "v_m", "v_y")] * true_variance_scale
    true_orientation <- sample(c(-1, 1), n, replace = TRUE)
    orientation_correct <- rbinom(n, 1, orientation_accuracy) == 1
    reported_orientation <- ifelse(orientation_correct, true_orientation,
                                   -true_orientation)
    ld <- jg02_block_ld(ld_block, ld_rho)

    g <- rnorm(n, sd = sqrt(variance$v_x))
    d <- rnorm(n, sd = sqrt(variance$v_m))
    e <- rnorm(n, sd = sqrt(variance$v_y))
    h_by_block <- if (q > 0) rbinom(blocks, 1, q) else rep(0, blocks)
    h <- rep(h_by_block, each = variants_per_block)
    true_m <- a * g + d
    true_y <- c_path * g + b * true_m + h * lambda * d +
        eta * true_orientation + e

    trait_correlation <- matrix(c(
        1, sampling_rho[["xm"]], sampling_rho[["xy"]],
        sampling_rho[["xm"]], 1, sampling_rho[["my"]],
        sampling_rho[["xy"]], sampling_rho[["my"]], 1
    ), 3, 3, byrow = TRUE)
    if (min(eigen(trait_correlation, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
        stop("sampling correlation matrix must be positive definite")
    }
    error <- matrix(0, n, 3)
    trait_factor <- t(chol(trait_correlation))
    for (current in unique(ld_block)) {
        index <- which(ld_block == current)
        ld_factor <- t(chol(ld[index, index, drop = FALSE]))
        z <- matrix(rnorm(length(index) * 3), nrow = length(index), ncol = 3)
        error[index, ] <- ld_factor %*% z %*% t(trait_factor) * se
    }

    data <- data.frame(
        variant = variant,
        ld_block = ld_block,
        role = role,
        beta_x = as.vector(ld %*% g) + error[, 1],
        se_x = se,
        beta_m = as.vector(ld %*% true_m) + error[, 2],
        se_m = se,
        beta_y = as.vector(ld %*% true_y) + error[, 3],
        se_y = se,
        v_x = variance$v_x * reported_variance_scale / true_variance_scale,
        v_m = variance$v_m * reported_variance_scale / true_variance_scale,
        v_y = variance$v_y * reported_variance_scale / true_variance_scale,
        orientation = reported_orientation,
        orientation_probability = orientation_accuracy,
        rho_xm = reported_sampling_rho[["xm"]],
        rho_xy = reported_sampling_rho[["xy"]],
        rho_my = reported_sampling_rho[["my"]],
        stringsAsFactors = FALSE
    )
    list(data = data, ld = ld)
}

jg02_write_fixture <- function(fixture, input, ld_path) {
    write.table(fixture$data, input, sep = "\t", quote = FALSE,
                row.names = FALSE)
    ld_table <- data.frame(variant = fixture$data$variant, fixture$ld,
                           check.names = FALSE)
    names(ld_table)[-1] <- fixture$data$variant
    write.table(ld_table, ld_path, sep = "\t", quote = FALSE,
                row.names = FALSE)
}

jg02_log_mvn <- function(observed, mean, covariance) {
    factor <- chol(covariance)
    residual <- observed - mean
    solved <- forwardsolve(t(factor), residual)
    -0.5 * (length(observed) * log(2 * pi) +
            2 * sum(log(diag(factor))) + sum(solved^2))
}

jg02_loglik <- function(data, ld, a, b, c_path, lambda, q, eta) {
    log_sum_exp <- function(value) {
        maximum <- max(value)
        maximum + log(sum(exp(value - maximum)))
    }
    values <- numeric()
    for (current in unique(data$ld_block)) {
        index <- which(data$ld_block == current)
        r <- ld[index, index, drop = FALSE]
        n <- length(index)
        kg <- r %*% diag(data$v_x[index], n) %*% r
        kd <- r %*% diag(data$v_m[index], n) %*% r
        ky <- r %*% diag(data$v_y[index], n) %*% r
        observed <- c(data$beta_x[index], data$beta_m[index], data$beta_y[index])
        sampling <- matrix(0, 3 * n, 3 * n)
        for (i in seq_len(n)) {
            for (j in seq_len(n)) {
                sampling[i, j] <- r[i, j] * data$se_x[index[i]] * data$se_x[index[j]]
                sampling[n + i, n + j] <-
                    r[i, j] * data$se_m[index[i]] * data$se_m[index[j]]
                sampling[2 * n + i, 2 * n + j] <-
                    r[i, j] * data$se_y[index[i]] * data$se_y[index[j]]
                sampling[i, n + j] <- data$rho_xm[index[i]] * r[i, j] *
                    data$se_x[index[i]] * data$se_m[index[j]]
                sampling[n + j, i] <- sampling[i, n + j]
                sampling[i, 2 * n + j] <- data$rho_xy[index[i]] * r[i, j] *
                    data$se_x[index[i]] * data$se_y[index[j]]
                sampling[2 * n + j, i] <- sampling[i, 2 * n + j]
                sampling[n + i, 2 * n + j] <- data$rho_my[index[i]] * r[i, j] *
                    data$se_m[index[i]] * data$se_y[index[j]]
                sampling[2 * n + j, n + i] <- sampling[n + i, 2 * n + j]
            }
        }
        component <- function(h, actual_orientation) {
            bg <- c_path + a * b
            bd <- b + h * lambda
            covariance <- sampling
            covariance[seq_len(n), seq_len(n)] <-
                covariance[seq_len(n), seq_len(n)] + kg
            covariance[seq_len(n), n + seq_len(n)] <-
                covariance[seq_len(n), n + seq_len(n)] + a * kg
            covariance[n + seq_len(n), seq_len(n)] <-
                t(covariance[seq_len(n), n + seq_len(n)])
            covariance[seq_len(n), 2 * n + seq_len(n)] <-
                covariance[seq_len(n), 2 * n + seq_len(n)] + bg * kg
            covariance[2 * n + seq_len(n), seq_len(n)] <-
                t(covariance[seq_len(n), 2 * n + seq_len(n)])
            covariance[n + seq_len(n), n + seq_len(n)] <-
                covariance[n + seq_len(n), n + seq_len(n)] + a^2 * kg + kd
            covariance[n + seq_len(n), 2 * n + seq_len(n)] <-
                covariance[n + seq_len(n), 2 * n + seq_len(n)] +
                a * bg * kg + bd * kd
            covariance[2 * n + seq_len(n), n + seq_len(n)] <-
                t(covariance[n + seq_len(n), 2 * n + seq_len(n)])
            covariance[2 * n + seq_len(n), 2 * n + seq_len(n)] <-
                covariance[2 * n + seq_len(n), 2 * n + seq_len(n)] +
                bg^2 * kg + bd^2 * kd + ky
            directional_mean <- as.vector(r %*% actual_orientation)
            mean <- c(rep(0, 2 * n), eta * directional_mean)
            jg02_log_mvn(observed, mean, covariance)
        }
        configurations <- expand.grid(rep(list(c(FALSE, TRUE)), n))
        orientation_terms <- numeric(nrow(configurations) * 2)
        position <- 0L
        for (configuration in seq_len(nrow(configurations))) {
            correct <- as.logical(configurations[configuration, ])
            probability <- data$orientation_probability[index]
            log_weight <- sum(ifelse(correct, log(probability), log1p(-probability)))
            actual_orientation <- data$orientation[index] * ifelse(correct, 1, -1)
            position <- position + 1L
            orientation_terms[position] <- log1p(-q) + log_weight +
                component(0, actual_orientation)
            position <- position + 1L
            orientation_terms[position] <- log(q) + log_weight +
                component(1, actual_orientation)
        }
        values <- c(values, log_sum_exp(orientation_terms))
    }
    sum(values)
}
