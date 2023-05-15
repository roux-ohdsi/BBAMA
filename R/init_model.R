#' @title init_model
#' @description
#' @param data
#' @param id
#' @param ...
#' @return
#' @examples
#' @export

init_model <- function(data, id, ...){
    id <- substitute(id)
    id_col <- eval(id, data)
    if (anyNA(id_col)) stop("Some ids are missing.")
    if (anyDuplicated(id_col)) stop("The ids for the main effects of interest should be unique. Please pass the negative controls....")

    data$.id <- id_col
    data$.site <- as.numeric(factor(data$.id))
    data <- data[order(data$.site),]

    structure(list(
        data = data,
        id = id,
        stan_data = list(dist_main = 0, dist_bias = 0, fixed_main = 0, fixed_bias = 0,
                         M = nrow(data), N = 0, site = integer(0),
                         y_int = integer(0), x_int = integer(0),
                         y_star = numeric(0), x_star = numeric(0),
                         y = numeric(0), x = numeric(0), s_0 = numeric(0), s_j = numeric(0),
                         n_vals = 0, ll = array(numeric(0), dim = c(0,0)), vals_evaled = numeric(0), min_val = 0, max_val = 0,
                         K = 0, p = array(numeric(0), dim = c(0,0)), L = 0, q = array(numeric(0), dim = c(0,0)),
                         mu_mean = 0, mu_sd = 10, tau_mean = 0, tau_sd = 10, lambda_mean = 0,
                         lambda_sd = 10, eta_mean = 0, eta_sd = 10, gamma_mean = 0, gamma_sd = 10),
        call = deparse(match.call(), width.cutoff = 500),
        likelihood_call = NULL,
        likelihood = NULL,
        has_likelihood = 0,
        has_negative_control_likelihood = 0
    ),
        class = "bbama_model")
}
