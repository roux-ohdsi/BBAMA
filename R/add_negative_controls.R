#' @title add_negative_controls
#' @description
#' @param bbama_model
#' @param data
#' @param id
#' @param likelihood
#' @param fixed
#' @param y1
#' @param y0
#' @param pt1
#' @param pt0
#' @param est
#' @param se
#' @param same_args
#' @return
#' @examples
#' @export

add_negative_controls <- function(bbama_model, data, id = bbama_model$id, likelihood = bbama_model$likelihood,
                                  fixed = NULL,
                                  y1 = NULL, y0 = NULL, pt1 = NULL, pt0 = NULL,
                                  est = NULL, se = NULL, same_args = TRUE) {

    ### DEAL WITH IDS
    id <- substitute(id)

    if (same_args) {
        if (!identical(eval(id), bbama_model$id)) stop("Which id do you want?")
        if (!identical(eval(likelihood), bbama_model$likelihood)) stop("Which likelihood do you want?")
        likelihood <- bbama_model$likelihood
        likelihood_call <- bbama_model$likelihood_call
    } else {
        likelihood_call <- match.call()
    }

    id_col <- eval(eval(id), data)
    if (anyNA(id_col)) stop("Some ids are missing.")

    main_ids <- eval(bbama_model$id, bbama_model$data)
    ids_not_in_main <- setdiff(id_col, main_ids)
    ids_not_in_ncs <- setdiff(main_ids, id_col)
    if (length(ids_not_in_main) > 0) warning("There are negative controls for id(s) {",
                                             paste(ids_not_in_main, collapse = ", "),
                                             "} but no corresponding main effect(s). ",
                                             "These negative controls will be ignored.")
    if (length(ids_not_in_ncs) > 0) warning("There are no negative controls for id(s) {",
                                            paste(ids_not_in_ncs, collapse = ", "),
                                            "}. Check the id values if this is unexpected.")

    data$.site <- id_col
    data <- subset(data, !.site %in% ids_not_in_main)
    data$.site <- as.numeric(factor(data$.site))
    data <- data[order(data$.site),]

    ### DEAL WITH LIKELIHOOD
    if (likelihood == "grid") stop("Grid not supported for negative controls at the moment")
    bbama_model$stan_data$dist_bias <- switch(likelihood, poisson = 1, normal = 2)

    need_args <- switch(likelihood, poisson = c("y1", "y0", "pt1", "pt0"),
                        normal = c("est", "se"))
    got_args <- as.list(likelihood_call)[need_args]

    if (any(sapply(got_args, is.null))) stop("Missing required argument for this likelihood")

    need_data <- lapply(substitute(got_args), eval, data)

    if (likelihood == "poisson") {
        need_data$x_star <- need_data$pt1*(need_data$y0/need_data$pt0)
        need_data$x_int <- need_data$y1
        need_data <- need_data[c("x_int", "x_star")]
    } else {
        names(need_data) <- c("x", "s_j")
    }

    bbama_model$stan_data <- modifyList(bbama_model$stan_data, need_data)
    bbama_model$stan_data$N <- nrow(data)
    bbama_model$stan_data$site <- data$.site
    bbama_model$has_negative_control_likelihood <- 1

    ### DEAL WITH FIXED EFFECTS
    if (!is.null(fixed)) {
        bbama_model$stan_data$fixed_bias <- 1
        bbama_model$stan_data$q <- as.matrix(model.matrix(fixed, bbama_model$data)[,-1])
        bbama_model$stan_data$L <- ncol(bbama_model$stan_data$q)
    }

    ca <- match.call()
    ca$bbama_model <- NULL
    bbama_model$call <- c(bbama_model$call, deparse(ca, width.cutoff = 500))

    return(bbama_model)

}
