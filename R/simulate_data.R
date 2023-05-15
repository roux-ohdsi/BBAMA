#' @title simulate_data
#' @description
#' @param n_outcomes
#' @param databases
#' @param x1_databases
#' @param x2_databases
#' @param n_per_database
#' @param min_rate
#' @param max_rate
#' @param obs_days
#' @param prop_included
#' @param mu
#' @param tau
#' @param lambda
#' @param eta
#' @param min_gamma
#' @param max_gamma
#' @param x1_fixed
#' @param x2_fixed
#' @param include_ll
#' @return
#' @examples
#' @importFrom Cyclops createCyclopsData fitCyclopsModel getCyclopsProfileLogLikelihood
#' @export

simulate_data <- function(n_outcomes = 100, databases = c("A", "B", "C", "D", "E"),
                          x1_databases = c(), x2_databases = c(),
                          n_per_database = round(100000*runif(length(databases))),
                          min_rate = 0, max_rate = 0.1, obs_days = 30, prop_included = 1,
                          mu = log(1.5), tau = 0.2, lambda = log(1.25), eta = 0.2,
                          min_gamma = 0.1, max_gamma = 0.3, x1_fixed = 0, x2_fixed = 0,
                          include_ll = TRUE) {

    outcomes <- 1:n_outcomes
    base_rates <- runif(n_outcomes, min_rate, max_rate)
    negativeControls <- c(FALSE, rep(TRUE, length(outcomes) - 1))
    n_databases <- length(databases)

    delta_i <- rnorm(n_databases, lambda + (databases %in% x2_databases)*x2_fixed, eta)
    gamma_i <- runif(n_databases, min_gamma, max_gamma)

    combos <- expand.grid(outcomeId = outcomes, databaseId = databases)
    dat <- merge(
        merge(combos,
              data.frame(databaseId = databases, n = n_per_database, delta_i, gamma_i),
              by = "databaseId"),
        data.frame(outcomeId = outcomes, rate_unexposed = base_rates, isNegativeControl = negativeControls),
        by = "outcomeId"
    )

    dat <- dat[sample(nrow(dat), prop_included*nrow(dat)),]
    dat$x1 <- ifelse(dat$databaseId %in% x1_databases, 1, 0)
    dat$x2 <- ifelse(dat$databaseId %in% x2_databases, 1, 0)
    dat$beta_ij <- rnorm(nrow(dat), dat$delta_i, dat$gamma_i)
    dat$theta_ij <- ifelse(dat$isNegativeControl, 0, rnorm(nrow(dat), mu + dat$x1*x1_fixed, tau))
    dat$theta_tilde_ij <- dat$beta_ij + dat$theta_ij

    seq_vec <- Vectorize(seq.default, vectorize.args = c("to"))

    dat$personTime_exposed <- sapply(seq_vec(1, dat$n), sample, 1)*obs_days
    dat$personTime_unexposed <- dat$n*obs_days - dat$personTime_exposed
    dat$outcomes_unexposed <- rpois(nrow(dat), dat$rate_unexposed*dat$personTime_unexposed)
    dat$outcomes_exposed <- rpois(nrow(dat), dat$rate_unexposed*dat$personTime_exposed*exp(dat$theta_tilde_ij))

    if (include_ll) {
        cyclops_data <- apply(dat[,c("outcomes_exposed", "outcomes_unexposed",
                                     "personTime_exposed", "personTime_unexposed")], 1, function(.x) {
                                         newdat <- data.frame(outcomes = c(.x["outcomes_exposed"],
                                                                           .x["outcomes_unexposed"]),
                                                              exposure = c(1, 0),
                                                              logpersonTime = c(log(.x["personTime_exposed"]),
                                                                                log(.x["personTime_unexposed"])))
                                         Cyclops::createCyclopsData(outcomes ~ exposure, data = newdat, time = newdat$logpersonTime,
                                                                    modelType = "pr")
                                     })

        cyclops_fits <- lapply(cyclops_data, Cyclops::fitCyclopsModel)

        dat$logRr <- sapply(cyclops_fits, function(.x) tryCatch(coef(.x)[["exposure"]], error = function(e) NA))
        dat$seLogRr <- sapply(cyclops_fits, function(.x) tryCatch(diff(confint(.x, "exposure")[c(2,3)]) /
                                                                      (qnorm(0.975) * 2), error = function(e) NA))

        dat$ll <- lapply(cyclops_fits, function(.x) tryCatch(Cyclops::getCyclopsProfileLogLikelihood(.x, parm = "exposure",
                         bounds = c(log(.1), log(10))), error = function(e) NULL))
        dat <- subset(dat, !(is.na(logRr) | is.na(seLogRr) | is.null(ll)))
        # necessary right now due to bug in Cyclops
        dat <- subset(dat, seLogRr > 0)
    } else {
        glm_data <- apply(dat[,c("outcomes_exposed", "outcomes_unexposed",
                                     "personTime_exposed", "personTime_unexposed")], 1, function(.x) {
                                        data.frame(outcomes = c(.x["outcomes_exposed"],
                                                                           .x["outcomes_unexposed"]),
                                                              exposure = c(1, 0),
                                                              logpersonTime = c(log(.x["personTime_exposed"]),
                                                                                log(.x["personTime_unexposed"])))
                                     })
        glm_fits <- lapply(glm_data, function(.x) tryCatch(glm(outcomes ~ exposure, offset = logpersonTime, family = poisson(),
                                                      data = .x), error = function(e) NULL))
        dat$logRr <- sapply(glm_fits, function(.x) tryCatch(.x$coefficients[2], error = function(e) NA))
        dat$seLogRr <- sapply(glm_fits, function(.x) tryCatch(sqrt(chol2inv(.x$qr$qr)[2,2]), error = function(e) NA))
        dat <- subset(dat, !(is.na(logRr) | is.na(seLogRr)))
    }

    dat <- dat[order(dat$outcomeId),]
    dat <- dat[order(dat$databaseId),]
    dat
}
