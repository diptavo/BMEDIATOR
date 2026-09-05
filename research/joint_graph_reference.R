#!/usr/bin/env Rscript

# Transparent reference implementation of docs/JOINT_GRAPH_MODEL_V0_1.md.
# Base R only. The finite quadrature is exact for the frozen discrete reference
# approximation and is intentionally kept separate from production code.

jg_logsumexp <- function(x) {
    maximum <- max(x)
    if (!is.finite(maximum)) return(maximum)
    maximum + log(sum(exp(x - maximum)))
}

jg_logaddexp <- function(x, y) {
    maximum <- max(x, y)
    if (!is.finite(maximum)) return(maximum)
    maximum + log(exp(x - maximum) + exp(y - maximum))
}

jg_gh7 <- function(sd) {
    nodes <- c(
        -2.6519613568352335,
        -1.6735516287674714,
        -0.8162878828589647,
         0.0,
         0.8162878828589647,
         1.6735516287674714,
         2.6519613568352335
    )
    weights <- c(
        0.0009717812450995192,
        0.05451558281912703,
        0.4256072526101278,
        0.8102646175568073,
        0.4256072526101278,
        0.05451558281912703,
        0.0009717812450995192
    )
    list(
        value = sqrt(2) * sd * nodes,
        log_weight = log(weights / sqrt(pi))
    )
}

jg_default_config <- function() {
    list(
        pi_xm = 0.25,
        pi_my = 0.10,
        pi_pleio = 0.10,
        prior_sd_a = 0.70,
        prior_sd_b = 0.70,
        prior_sd_c = 0.175,
        prior_sd_lambda = 0.70,
        q_values = c(0.15, 0.35, 0.60),
        q_log_weights = rep(-log(3), 3),
        vx = c(A = 0.0400, B = 0.0004, C = 0.0400),
        vm = c(A = 0.0025, B = 0.0400, C = 0.0400),
        vy = 0.0025
    )
}

jg_validate_data <- function(data) {
    required <- c("variant", "role", "beta_x", "se_x", "beta_m", "se_m",
                  "beta_y", "se_y")
    missing <- setdiff(required, names(data))
    if (length(missing)) {
        stop("missing joint-graph columns: ", paste(missing, collapse = ", "))
    }
    if (!nrow(data)) stop("joint-graph input contains no variants")
    if (any(!data$role %in% c("A", "B", "C"))) {
        stop("joint-graph roles must be A, B, or C")
    }
    numeric_columns <- required[required != "variant" & required != "role"]
    for (column in numeric_columns) {
        if (any(!is.finite(data[[column]]))) {
            stop("non-finite value in column ", column)
        }
    }
    if (any(data$se_x <= 0 | data$se_m <= 0 | data$se_y <= 0)) {
        stop("joint-graph standard errors must be positive")
    }
    invisible(TRUE)
}

jg_log_mvn3_zero <- function(x1, x2, x3, s11, s12, s13, s22, s23, s33) {
    determinant <- s11 * (s22 * s33 - s23 * s23) -
        s12 * (s12 * s33 - s13 * s23) +
        s13 * (s12 * s23 - s13 * s22)
    if (!is.finite(determinant) || determinant <= 0) return(-Inf)

    i11 <- (s22 * s33 - s23 * s23) / determinant
    i12 <- (s13 * s23 - s12 * s33) / determinant
    i13 <- (s12 * s23 - s13 * s22) / determinant
    i22 <- (s11 * s33 - s13 * s13) / determinant
    i23 <- (s12 * s13 - s11 * s23) / determinant
    i33 <- (s11 * s22 - s12 * s12) / determinant

    quadratic <- i11 * x1 * x1 + i22 * x2 * x2 + i33 * x3 * x3 +
        2 * i12 * x1 * x2 + 2 * i13 * x1 * x3 + 2 * i23 * x2 * x3
    -0.5 * (3 * log(2 * pi) + log(determinant) + quadratic)
}

jg_component_loglik <- function(data, a, b, c_path, lambda, h, config) {
    vx <- unname(config$vx[data$role])
    vm <- unname(config$vm[data$role])
    bg <- c_path + a * b
    bd <- b + h * lambda

    s11 <- vx + data$se_x^2
    s12 <- a * vx
    s13 <- bg * vx
    s22 <- a^2 * vx + vm + data$se_m^2
    s23 <- a * bg * vx + bd * vm
    s33 <- bg^2 * vx + bd^2 * vm + config$vy + data$se_y^2

    vapply(seq_len(nrow(data)), function(index) {
        jg_log_mvn3_zero(
            data$beta_x[index], data$beta_m[index], data$beta_y[index],
            s11[index], s12[index], s13[index], s22[index], s23[index],
            s33[index]
        )
    }, numeric(1))
}

jg_loglik <- function(data, a, b, c_path, lambda, q, config) {
    clean <- jg_component_loglik(data, a, b, c_path, lambda, 0, config)
    if (q <= 0 || lambda == 0) return(sum(clean))
    shared <- jg_component_loglik(data, a, b, c_path, lambda, 1, config)
    sum(mapply(function(first, second) {
        jg_logaddexp(log1p(-q) + first, log(q) + second)
    }, clean, shared))
}

jg_states <- function() {
    data.frame(
        state = c("S000", "S100", "S010", "S110",
                  "S001", "S101", "S011", "S111"),
        z_xm = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L),
        z_my = c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L),
        z_pleio = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
        stringsAsFactors = FALSE
    )
}

jg_state_log_prior <- function(z_xm, z_my, z_pleio, config) {
    bernoulli <- function(z, probability) {
        if (z == 1L) log(probability) else log1p(-probability)
    }
    bernoulli(z_xm, config$pi_xm) +
        bernoulli(z_my, config$pi_my) +
        bernoulli(z_pleio, config$pi_pleio)
}

jg_parameter_grid <- function(active, quadrature) {
    if (active) quadrature else list(value = 0, log_weight = 0)
}

jg_integrate_state <- function(data, state, config, quadrature) {
    a_grid <- jg_parameter_grid(state$z_xm == 1L, quadrature$a)
    b_grid <- jg_parameter_grid(state$z_my == 1L, quadrature$b)
    c_grid <- quadrature$c_path
    lambda_grid <- jg_parameter_grid(state$z_pleio == 1L, quadrature$lambda)
    if (state$z_pleio == 1L) {
        q_values <- config$q_values
        q_log_weights <- config$q_log_weights
    } else {
        q_values <- 0
        q_log_weights <- 0
    }

    count <- length(a_grid$value) * length(b_grid$value) *
        length(c_grid$value) * length(lambda_grid$value) * length(q_values)
    log_terms <- numeric(count)
    a_terms <- numeric(count)
    b_terms <- numeric(count)
    c_terms <- numeric(count)
    lambda_terms <- numeric(count)
    q_terms <- numeric(count)
    position <- 0L

    for (ia in seq_along(a_grid$value)) {
        for (ib in seq_along(b_grid$value)) {
            for (ic in seq_along(c_grid$value)) {
                for (il in seq_along(lambda_grid$value)) {
                    for (iq in seq_along(q_values)) {
                        position <- position + 1L
                        a <- a_grid$value[ia]
                        b <- b_grid$value[ib]
                        c_path <- c_grid$value[ic]
                        lambda <- lambda_grid$value[il]
                        q <- q_values[iq]
                        log_terms[position] <-
                            a_grid$log_weight[ia] + b_grid$log_weight[ib] +
                            c_grid$log_weight[ic] + lambda_grid$log_weight[il] +
                            q_log_weights[iq] +
                            jg_loglik(data, a, b, c_path, lambda, q, config)
                        a_terms[position] <- a
                        b_terms[position] <- b
                        c_terms[position] <- c_path
                        lambda_terms[position] <- lambda
                        q_terms[position] <- q
                    }
                }
            }
        }
    }

    log_evidence <- jg_logsumexp(log_terms)
    weights <- exp(log_terms - log_evidence)
    list(
        log_evidence = log_evidence,
        mean_a = sum(weights * a_terms),
        mean_b = sum(weights * b_terms),
        mean_c = sum(weights * c_terms),
        mean_lambda = sum(weights * lambda_terms),
        mean_q = sum(weights * q_terms),
        mean_indirect = sum(weights * a_terms * b_terms)
    )
}

jg_fit <- function(data, config = jg_default_config()) {
    jg_validate_data(data)
    quadrature <- list(
        a = jg_gh7(config$prior_sd_a),
        b = jg_gh7(config$prior_sd_b),
        c_path = jg_gh7(config$prior_sd_c),
        lambda = jg_gh7(config$prior_sd_lambda)
    )
    states <- jg_states()
    fits <- vector("list", nrow(states))
    log_joint <- numeric(nrow(states))

    for (index in seq_len(nrow(states))) {
        fits[[index]] <- jg_integrate_state(data, states[index, ], config, quadrature)
        log_joint[index] <- fits[[index]]$log_evidence + jg_state_log_prior(
            states$z_xm[index], states$z_my[index], states$z_pleio[index], config
        )
    }

    log_normalizer <- jg_logsumexp(log_joint)
    posterior <- exp(log_joint - log_normalizer)
    state_pp <- setNames(posterior, states$state)
    weighted <- function(field, subset = rep(TRUE, nrow(states))) {
        denominator <- sum(posterior[subset])
        if (!(denominator > 0)) return(NaN)
        sum(posterior[subset] * vapply(fits[subset], `[[`, numeric(1), field)) /
            denominator
    }

    xm <- states$z_xm == 1L
    my <- states$z_my == 1L
    pleio <- states$z_pleio == 1L
    two_path <- xm & my
    result <- c(
        setNames(state_pp, paste0("PP_", names(state_pp))),
        PP_XM = sum(posterior[xm]),
        PP_global_MY = sum(posterior[my]),
        PP_nonaligned_P = sum(posterior[pleio]),
        PP_two_path = sum(posterior[two_path]),
        PP_two_path_plus_P = state_pp[["S111"]],
        mean_a_given_XM = weighted("mean_a", xm),
        mean_b_given_MY = weighted("mean_b", my),
        mean_c = weighted("mean_c"),
        mean_lambda_given_P = weighted("mean_lambda", pleio),
        mean_q_given_P = weighted("mean_q", pleio),
        mean_indirect_given_two_path = weighted("mean_indirect", two_path),
        log_evidence = log_normalizer
    )
    attr(result, "state_log_evidence") <- setNames(
        vapply(fits, `[[`, numeric(1), "log_evidence"), states$state
    )
    result
}

jg_simulate <- function(seed, n_per_role, a, b, c_path, lambda = 0, q = 0,
                        se = 0.02, config = jg_default_config()) {
    set.seed(seed)
    roles <- rep(c("A", "B", "C"), each = n_per_role)
    variants <- sprintf("v%04d", seq_along(roles))
    vx <- unname(config$vx[roles])
    vm <- unname(config$vm[roles])
    g <- rnorm(length(roles), sd = sqrt(vx))
    d <- rnorm(length(roles), sd = sqrt(vm))
    direct_y <- rnorm(length(roles), sd = sqrt(config$vy))
    h <- if (q > 0) rbinom(length(roles), 1, q) else rep(0, length(roles))
    true_m <- a * g + d
    true_y <- c_path * g + b * true_m + h * lambda * d + direct_y
    data.frame(
        variant = variants,
        role = roles,
        beta_x = g + rnorm(length(roles), sd = se),
        se_x = se,
        beta_m = true_m + rnorm(length(roles), sd = se),
        se_m = se,
        beta_y = true_y + rnorm(length(roles), sd = se),
        se_y = se,
        stringsAsFactors = FALSE
    )
}

jg_write_result <- function(result, output) {
    row <- as.data.frame(as.list(result), check.names = FALSE)
    row <- data.frame(
        model_version = "JG-0.1",
        identification_scope = "CONDITIONAL_ON_NO_EXACT_ALIGNED_PLEIOTROPY",
        row,
        check.names = FALSE
    )
    write.table(row, output, sep = "\t", quote = FALSE, row.names = FALSE,
                na = "nan")
}

jg_cli <- function(arguments) {
    if (length(arguments) != 2L) {
        stop("usage: Rscript research/joint_graph_reference.R INPUT.tsv OUTPUT.tsv")
    }
    data <- read.delim(arguments[[1]], check.names = FALSE,
                       stringsAsFactors = FALSE)
    jg_write_result(jg_fit(data), arguments[[2]])
}

if (sys.nframe() == 0L) {
    jg_cli(commandArgs(trailingOnly = TRUE))
}
