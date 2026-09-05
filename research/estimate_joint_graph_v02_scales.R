#!/usr/bin/env Rscript

jg02_scale_columns <- c(
    "role", "beta_x", "se_x", "beta_m", "se_m", "beta_y", "se_y",
    "rho_xm", "rho_xy", "rho_my"
)

jg02_validate_scale_data <- function(data) {
    missing <- setdiff(jg02_scale_columns, names(data))
    if (length(missing)) stop("missing columns: ", paste(missing, collapse = ", "))
    if (any(!data$role %in% c("A", "B", "C"))) stop("role must be A, B, or C")
    numeric_columns <- setdiff(jg02_scale_columns, "role")
    if (any(!vapply(data[numeric_columns], function(x) all(is.finite(x)), logical(1)))) {
        stop("scale input contains non-finite numeric values")
    }
    if (any(data$se_x <= 0 | data$se_m <= 0 | data$se_y <= 0)) {
        stop("standard errors must be positive")
    }
    counts <- table(factor(data$role, levels = c("A", "B", "C")))
    if (any(counts < 30L)) {
        stop("external scale estimation requires at least 30 variants per role")
    }
    invisible(TRUE)
}

jg02_mvn2_nll <- function(parameters, data) {
    a <- parameters[[1]]
    vx <- exp(parameters[2:4])
    vm <- exp(parameters[5:7])
    role_index <- match(data$role, c("A", "B", "C"))
    s11 <- vx[role_index] + data$se_x^2
    s12 <- a * vx[role_index] + data$rho_xm * data$se_x * data$se_m
    s22 <- a^2 * vx[role_index] + vm[role_index] + data$se_m^2
    determinant <- s11 * s22 - s12^2
    if (any(!is.finite(determinant)) || any(determinant <= 0)) return(1e100)
    quadratic <- (s22 * data$beta_x^2 - 2 * s12 * data$beta_x * data$beta_m +
                  s11 * data$beta_m^2) / determinant
    sum(log(2 * pi) + 0.5 * log(determinant) + 0.5 * quadratic)
}

jg02_mvn3_zero_nll <- function(observed, covariance) {
    determinant <- det(covariance)
    if (!is.finite(determinant) || determinant <= 0) return(1e100)
    0.5 * (3 * log(2 * pi) + log(determinant) +
           drop(crossprod(observed, solve(covariance, observed))))
}

jg02_outcome_nll <- function(parameters, data, a, vx, vm) {
    c_path <- parameters[[1]]
    b <- parameters[[2]]
    vy <- exp(parameters[[3]])
    role_index <- match(data$role, c("A", "B", "C"))
    total <- 0
    for (index in seq_len(nrow(data))) {
        current_vx <- vx[role_index[index]]
        current_vm <- vm[role_index[index]]
        bg <- c_path + a * b
        covariance <- matrix(c(
            current_vx + data$se_x[index]^2,
            a * current_vx + data$rho_xm[index] * data$se_x[index] * data$se_m[index],
            bg * current_vx + data$rho_xy[index] * data$se_x[index] * data$se_y[index],
            a * current_vx + data$rho_xm[index] * data$se_x[index] * data$se_m[index],
            a^2 * current_vx + current_vm + data$se_m[index]^2,
            a * bg * current_vx + b * current_vm +
                data$rho_my[index] * data$se_m[index] * data$se_y[index],
            bg * current_vx + data$rho_xy[index] * data$se_x[index] * data$se_y[index],
            a * bg * current_vx + b * current_vm +
                data$rho_my[index] * data$se_m[index] * data$se_y[index],
            bg^2 * current_vx + b^2 * current_vm + vy + data$se_y[index]^2
        ), 3, 3, byrow = TRUE)
        total <- total + jg02_mvn3_zero_nll(
            c(data$beta_x[index], data$beta_m[index], data$beta_y[index]),
            covariance
        )
    }
    total
}

jg02_estimate_scales <- function(data) {
    jg02_validate_scale_data(data)
    roles <- c("A", "B", "C")
    corrected_x <- pmax(data$beta_x^2 - data$se_x^2, 1e-8)
    initial_vx <- vapply(roles, function(role) {
        max(mean(corrected_x[data$role == role]), 1e-6)
    }, numeric(1))
    corrected_cross <- data$beta_x * data$beta_m -
        data$rho_xm * data$se_x * data$se_m
    initial_a <- sum(corrected_cross) / max(sum(corrected_x), 1e-8)
    residual_m <- data$beta_m - initial_a * data$beta_x
    residual_m_error <- data$se_m^2 + initial_a^2 * data$se_x^2 -
        2 * initial_a * data$rho_xm * data$se_x * data$se_m
    initial_vm <- vapply(roles, function(role) {
        index <- data$role == role
        max(mean(residual_m[index]^2 - residual_m_error[index]), 1e-6)
    }, numeric(1))
    first <- optim(
        c(initial_a, log(initial_vx), log(initial_vm)),
        jg02_mvn2_nll,
        data = data,
        method = "BFGS",
        control = list(maxit = 1000, reltol = 1e-10)
    )
    if (first$convergence != 0L || !is.finite(first$value)) {
        stop("external X/M scale optimization failed")
    }
    a <- first$par[[1]]
    vx <- exp(first$par[2:4])
    vm <- exp(first$par[5:7])

    design <- cbind(data$beta_x, data$beta_m)
    initial_slopes <- tryCatch(
        solve(crossprod(design), crossprod(design, data$beta_y)),
        error = function(error) c(0, 0)
    )
    residual_y <- data$beta_y - drop(design %*% initial_slopes)
    initial_vy <- max(mean(residual_y^2 - data$se_y^2), 1e-6)
    second <- optim(
        c(initial_slopes[[1]], initial_slopes[[2]], log(initial_vy)),
        jg02_outcome_nll,
        data = data,
        a = a,
        vx = vx,
        vm = vm,
        method = "BFGS",
        control = list(maxit = 1000, reltol = 1e-10)
    )
    if (second$convergence != 0L || !is.finite(second$value)) {
        stop("external outcome scale optimization failed")
    }
    vy <- exp(second$par[[3]])
    data.frame(
        role = roles,
        v_x = vx,
        v_m = vm,
        v_y = vy,
        estimated_a = a,
        estimated_b = second$par[[2]],
        estimated_c = second$par[[1]],
        n_role = as.integer(table(factor(data$role, levels = roles))),
        xm_convergence = first$convergence,
        y_convergence = second$convergence,
        stringsAsFactors = FALSE
    )
}

jg02_scale_cli <- function(arguments) {
    if (length(arguments) != 2L) {
        stop("usage: Rscript research/estimate_joint_graph_v02_scales.R INPUT.tsv OUTPUT.tsv")
    }
    data <- read.delim(arguments[[1]], stringsAsFactors = FALSE,
                       check.names = FALSE)
    result <- jg02_estimate_scales(data)
    write.table(result, arguments[[2]], sep = "\t", quote = FALSE,
                row.names = FALSE)
}

if (sys.nframe() == 0L) {
    jg02_scale_cli(commandArgs(trailingOnly = TRUE))
}
