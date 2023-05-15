#' @title set_likelihood
#' @description
#' @param bbama_model
#' @param likelihood
#' @param fixed
#' @param y1
#' @param y0
#' @param pt1
#' @param pt0
#' @param est
#' @param se
#' @param ll
#' @param point
#' @param value
#' @param n_vals
#' @param min_val
#' @param max_val
#' @param extreme_min_val
#' @param extreme_max_val
#' @param extrapolate_formula
#' @param sep
#' @param ...
#' @return
#' @examples
#' @export

set_likelihood <- function(bbama_model, likelihood = c("poisson", "normal", "grid"),
                           fixed = NULL,
                           y1 = NULL, y0 = NULL, pt1 = NULL, pt0 = NULL,
                           est = NULL, se = NULL,
                           ll = NULL, point = NULL, value = NULL,
                           n_vals = NULL, min_val = NULL, max_val = NULL,
                           extreme_min_val = log(1/100), extreme_max_val = log(100),
                           extrapolate_formula = value ~ poly(point, 4), sep = ";", ...) {
    if (class(bbama_model) != "bbama_model") stop("Must use the function bbama_data() first.")

    bbama_model$stan_data$dist_main <- switch(likelihood, poisson = 1, normal = 2, grid = 3)
    need_args <- switch(likelihood, poisson = c("y1", "y0", "pt1", "pt0"),
                        normal = c("est", "se"), grid = c("ll", "point", "value"))
    got_args <- as.list(match.call())[need_args]

    if (likelihood == "grid") {
        if (all(sapply(got_args, is.null))) stop("Missing required argument for this likelihood")
        if ("ll" %in% names(got_args)) {
            ll_dfs <- eval(substitute(ll), bbama_model$data)
            if (all(c("point", "value") %in% names(got_args))) {
                lapply(ll_dfs, function(.x) data.frame(point = eval(substitute(point), ll),
                                                       value = eval(substitute(value), ll)))
            } else if (all(c("point", "value") %in% names(ll_dfs[[1]]))) {
                warning("Provide names for variables in `ll` dataframes")
            }
        } else if ("value" %in% names(got_args)) {
            ll_values <- lapply(strsplit(eval(substitute(value), bbama_model$data), sep), as.numeric)
            if ("point" %in% names(got_args)) {
                ll_points <- lapply(strsplit(eval(substitute(point), bbama_model$data), sep), as.numeric)
            } else { # using min, max val
                n_vals_each <- sapply(ll_values, length)
                if (anyNA(as.logical(sapply(n_vals_each, function(.x) all.equal(.x, mean(n_vals_each)))))) {
                    stop("The same number of values aren't provided; must provide points")
                }
                ll_points <- seq(min_val, max_val, length.out = n_vals_each[1])
                ll_points <- rep(list(ll_points), length(ll_values))
            }
            ll_dfs <- Map(data.frame, point = ll_points, value = ll_values)
            if (bbama_model$stan_data$M == 1) {
                ll_dfs <- list(do.call(rbind, ll_dfs))
            }
        }
        # see if these hold for all
        # TODO: this is just going to use the grid of the first one
        # instead use the densest grid...
        vals_evaled <- ll_dfs[[1]]$point
        n_vals <- length(vals_evaled)
        min_val <- min(vals_evaled)
        max_val <- max(vals_evaled)

        # if likelihood is monotone for any one need to extrapolate
        # or else there will likely be sampling errors
        to_extrapolate <- sapply(ll_dfs, function(.x) which.max(.x$value) %in% c(1, nrow(.x)))
        # if the grid is not regular, interpolate first
        to_interpolate1 <- sapply(ll_dfs, function(.x) nrow(.x) != n_vals)
        to_interpolate2 <- as.logical(sapply(ll_dfs, function(.x) all.equal(min(.x$point), min_val)))
        to_interpolate3 <- as.logical(sapply(ll_dfs, function(.x) all.equal(max(.x$point), max_val)))
        # may be na
        to_interpolate <- any(c(to_interpolate1, to_interpolate2, to_interpolate3))

        if (any(to_extrapolate)) {
            warning("Likelihood monotone; appoximating likelihood function between ",
                    round(extreme_min_val, 2), " and ", round(extreme_max_val, 2), ". Consider changing arguments ",
                    "`extreme_min_val` and `extreme_vax_val` if that is not sufficient.")
            by <- vals_evaled[2] - vals_evaled[1]
            min_val <- extreme_min_val
            max_val <- extreme_max_val
            vals_evaled <- seq(min_val, max_val, by = by)
            n_vals <- length(vals_evaled)

            l_mods <- lapply(ll_dfs, function(.x) lm(extrapolate_formula, data = .x))
            pred_vals <- lapply(l_mods, predict, newdata = data.frame(point = vals_evaled))
            new_ll <- lapply(pred_vals, function(.x) data.frame(point = vals_evaled, value = .x))
        } else if (!is.na(to_interpolate) & to_interpolate) {
            new_ll <- lapply(ll_dfs, function(.x) setNames(data.frame(approx(.x$point, .x$value,
                                                                             xout = vals_evaled,
                                                                             rule = 2)), c("point", "value")))
        } else new_ll <- ll_dfs

        need_data <- list(ll = sapply(new_ll, function(.x) .x$value),
                          vals_evaled = vals_evaled,
                          n_vals = n_vals,
                          min_val = min_val,
                          max_val = max_val)
    } else {
        if (any(sapply(got_args, is.null))) stop("Missing required argument for this likelihood")

        need_data <- lapply(substitute(got_args), eval, bbama_model$data)

        if (likelihood == "poisson") {
            need_data$y_star <- need_data$pt1*(need_data$y0/need_data$pt0)
            need_data$y_int <- need_data$y1
            need_data <- need_data[c("y_int", "y_star")]
        } else {
            names(need_data) <- c("y", "s_0")
        }
    }
    bbama_model$stan_data <- modifyList(bbama_model$stan_data, need_data)
    bbama_model$likelihood_call <- match.call()
    bbama_model$likelihood = likelihood
    bbama_model$has_likelihood <- 1

    if (!is.null(fixed)) {
        bbama_model$stan_data$fixed_main <- 1
        bbama_model$stan_data$p <- as.matrix(model.matrix(fixed, bbama_model$data)[,-1])
        bbama_model$stan_data$K <- ncol(bbama_model$stan_data$p)
    }

    ca <- match.call()
    ca$bbama_model <- NULL
    bbama_model$call <- c(bbama_model$call, deparse(ca, width.cutoff = 500))

    return(bbama_model)
}
