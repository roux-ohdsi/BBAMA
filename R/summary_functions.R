#' @title get_draws
#' @description
#' @param res
#' @return
#' @importFrom tidybayes gather_draws
#' @examples
#' @export

get_draws <- function(res) {
    if (!inherits(res, "bbama_fit")) stop("Not the right kind of object")
    if (res$stanfit@par_dims$lambda < 1) {
        suppressWarnings(draws <- tidybayes::gather_draws(res$stanfit, mu, tau, theta_0[.id]))
    } else {
        suppressWarnings(draws <- tidybayes::gather_draws(res$stanfit, mu, tau, lambda[j], eta[j], theta_0[.id]))
        draws$j <- NULL
    }
    mod <- attr(res, "bbama_model")
    draws$.id <- factor(draws$.id, levels = mod$data$.site, labels = mod$data$.id)
    class(draws) <- c("bbama_draws", class(draws))
    draws
}

#' @title get_summary
#' @description
#' @param res
#' @param measures
#' @return
#' @importFrom posterior summarise_draws quantile2
#' @export
#' @examples

get_summary <- function(res, ...) {
    if (inherits(res, "bbama_fit")) res <- get_draws(res)
    a <- list(...)
    if (length(a) == 0) {
    df <- posterior::summarise_draws(res, mean, sd, median,
                 ~posterior::quantile2(.x, probs = c(0.025, 0.975)))
    } else df <- posterior::summarise_draws(res, ...)
    df$variable <- NULL
    df
}

#' @title get_mcse
#' @description
#' @param res
#' @return
#' @importFrom posterior default_mcse_measures
#' @examples
#' @export

get_mcse <- function(res) {
    get_summary(res, posterior::default_mcse_measures())
}

#' @title get_convergence
#' @description
#' @param res
#' @return
#' @importFrom posterior default_convergence_measures
#' @examples
#' @export

get_convergence <- function(res) {
    get_summary(res, posterior::default_convergence_measures())
}

#' @method print bbama_model
#' @importFrom pillar style_subtle
#' @export
print.bbama_model <- function(x, ...) {
    cat(pillar::style_subtle("# A BBAMA model called with:\n"))
    cat(x$call, sep = " |>\n  ")
    return(invisible(x))
}

#' @method print bbama_fit
#' @importFrom pillar style_subtle
#' @importFrom utils getFromNamespace
#' @export
print.bbama_fit <- function(x, ...) {
    print_stan <- utils::getFromNamespace("print.stanfit", "rstan")
    to_print <- capture.output(print_stan(x$stanfit, ...))
    inds <- which(to_print == "")
    to_print[1:(inds[1] - 1)] <- pillar::style_subtle(paste("#", to_print[1:(inds[1] - 1)]))
    to_print[(inds[2] + 1):length(to_print)] <- pillar::style_subtle(paste("#", to_print[(inds[2] + 1):length(to_print)]))
    cat(to_print, sep = "\n")
    return(invisible(x))
}

#' @method print bbama_draws
#' @importFrom pillar style_subtle
#' @export
print.bbama_draws <- function(x, ...) {

    niterations <- length(unique(x$.iteration))
    nchains <- length(unique(x$.chain))
    v <- unique(x$.variable)
    ids <- levels(as.factor(x$.id))
    header <- paste0("# A bbama_draws object with ", nchains, " chains of ",
                     niterations, " iterations each\n", "# with values for variables ",
                     paste0(v, collapse = ", "), "\n# for sources ", paste0(ids, collapse = ", "), "\n")
    cat(pillar::style_subtle(header)) # make it match the tibble printing
    class(x) <- class(x)[class(x) != "bbama_draws"]
    print(x)
    return(invisible(x))
}

#' @method summary bbama_model
#' @export
summary.bbama_model <- function(object, ...) {
    p <- pillar::style_subtle
    object$stan_data$dist_main <- switch(object$stan_data$dist_main, "1" = "poisson", "2" = "normal", "3" = "grid-approximated")
    object$stan_data$dist_bias <- switch(object$stan_data$dist_bias, "1" = "poisson", "2" = "normal", "3" = "grid-approximated")
    has_fixed_main <- object$stan_data$K > 0
    has_fixed_bias <- object$stan_data$L > 0
    to_print <- paste0(p("# A BBAMA model:\n"),
                       ifelse(object$has_likelihood, paste0(p("#  - combines estimates from "), object$stan_data$M, p(" sources assuming a "),
                                                            object$stan_data$dist_main, p(" likelihood\n")), ""),
                       ifelse(object$has_negative_control_likelihood, paste0(p("#  - adjusts for bias using "), object$stan_data$N,
                                                                             p(" negative controls assuming a "), object$stan_data$dist_bias, p(" likelihood\n")), ""),
                       ifelse(has_fixed_main, paste0(p("#  - includes "), object$stan_data$K, p(" fixed effects for the main effect\n")), ""),
                       ifelse(has_fixed_bias, paste0(p("#  - includes "), object$stan_data$L, p(" fixed effects for the bias effect\n")), ""),
                       p("#  - uses as priors for the overall effect: mu ~ normal("), object$stan_data$mu_mean, p(", "), object$stan_data$mu_sd,
                       p("^2) and tau ~ normal("), object$stan_data$tau_mean, p(", "), object$stan_data$tau_sd, p("^2)\n"),
                       ifelse(object$has_negative_control_likelihood, paste0(p("#  - uses as priors for the overall bias: lambda ~ normal("),
                                                                             object$stan_data$lambda_mean, p(", "), object$stan_data$lambda_sd,
                                                                             p("^2) and eta ~ normal("), object$stan_data$eta_mean, p(", "), object$stan_data$eta_sd,
                                                                             p("^2)\n#    and for the data source-specific bias: gamma ~ normal("), object$stan_data$gamma_mean, p(", "), object$stan_data$gamma_sd, p("^2)")), ""))
    cat(to_print)
    # maybe rename these at some point so they're readable
    return(invisible(list(model = stanmodels$bbama, data = object$stan_data)))
}

#' @method summary bbama_fit
#' @export
summary.bbama_fit <- function(object, ...) {
    get_summary(object, "mean", "mcse_mean", "rhat", "ess_bulk", "ess_tail")
}

#' @method summary bbama_draws
#' @export
summary.bbama_draws <- function(object, ...) {
    get_summary(object)
}




