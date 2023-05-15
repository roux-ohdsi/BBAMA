#' @title set_priors
#' @description
#' @param bbama_model
#' @param ...
#' @return
#' @examples
#' @export

set_priors <- function(bbama_model, ...) {
    priors <- list(...)
    priors_ch <- lapply(priors, function(x) strsplit(deparse1(x) , "~")[[1]])
    new_params <- lapply(priors_ch, function(x) prior_to_data(x[1], x[2]))

    bbama_model$stan_data <- modifyList(bbama_model$stan_data, as.list(unlist(new_params)))

    ca <- match.call()
    ca$bbama_model <- NULL
    bbama_model$call <- c(bbama_model$call, deparse(ca, width.cutoff = 500))

    return(bbama_model)
}

#' @title prior_to_data
#' @description
#' @param parameter
#' @param prior
#' @return
#' @examples

prior_to_data <- function(parameter, prior) {
    parameter <- match.arg(trimws(parameter), c("mu", "tau", "lambda", "eta", "gamma"))
    is_normal <- grepl("normal\\s*\\(\\s*\\d+\\s*\\,\\s*\\d+\\s*\\)", prior, ignore.case = TRUE)
    if (!is_normal) stop("Non-normal priors are not yet supported. Please use 'parameter ~ normal(mean, sd)'.")
    digits <- regmatches(prior, gregexpr("(\\d+)", prior))[[1]]
    new_params <- list(as.numeric(digits[1]), as.numeric(digits[2]))
    names(new_params) <- c(paste0(parameter, "_mean"), paste0(parameter, "_sd"))
    new_params
}
