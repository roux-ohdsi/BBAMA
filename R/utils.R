#' Nicely wrap a message
#' @noRd
wrapmessage <- function(mess, width = 0.9 * getOption("width")) {
    message(paste(strwrap(paste(mess, collapse = " "), width = width), collapse = "\n"))
}
