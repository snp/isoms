#' rSquared
#'
#' calculates R^2 for the fitting model
#'
#' @param model object of class nls or lm, or similar providing \code{predict}
#' function.
#' @param y numeric vector of dependent variable used to make \code{model}.
#'
#' @return single numeric R^2 value, or NA if \code{class(model) == "try-error"}.
#' @export
#'
#' @examples
rSquared <- function(model, y) {
    ## Determine R2 of a fitted model
    if (class(model) != "try-error") {
        ssTot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
        ssRes <- sum((y - predict(model))^2, na.rm = TRUE)
        r2 <- 1 - ssRes/ssTot
    } else {
        r2 <- NA_real_
    }
    return(r2)
}
