#' @title bbama_sample
#' @description
#' @param bbama_model
#' @param pars
#' @param ...
#' @return
#' @examples
#' @export
#' @importFrom rstan sampling

# check for missing data here?
bbama_sample <- function(bbama_model, pars = c("mu", "tau", "lambda", "eta", "theta_0"), ...) {
    # if no changes have been made to the parameter choice
    # automatically add in if there are fixed effects
    if (identical(pars, c("mu", "tau", "lambda", "eta", "theta_0")) &
        bbama_model$stan_data$fixed_main == 1) {
        pars <- c(pars, "alpha")
        pars[pars == "theta_0"] <- "theta_transformed"
    }

    bbama_data <- bbama_model$stan_data
    # TODO: automatically add alpha etc par if fixed
    # nothing as of right now allows for different gamma params; all will be length 1
    if (length(bbama_data$gamma_mean) == 1) bbama_data$gamma_mean <- rep(bbama_data$gamma_mean, bbama_data$M)
    if (length(bbama_data$gamma_sd) == 1) bbama_data$gamma_sd <- rep(bbama_data$gamma_sd, bbama_data$M)
    to_array <- names(bbama_data) %in% c("y_int", "x_int", "y_star", "x_star",
                                         "y", "x", "s_0", "s_j", "ll", "vals_evaled",
                                         "p", "q", "gamma_mean", "gamma_sd")
    array_elems <- bbama_data[to_array & sapply(bbama_data, function(.x) length(.x) > 0)]
    bbama_data[to_array & sapply(bbama_data, function(.x) length(.x) > 0)] <- lapply(array_elems, as.array)
    samples <- rstan::sampling(stanmodels$bbama, data = bbama_data, pars = pars, ...)
    structure(list(stanfit = samples),
              class = "bbama_fit",
              bbama_model = bbama_model)
}
