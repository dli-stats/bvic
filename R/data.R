#' Test data generated based on clayton copula 
#'
#' Longer description of the data.
#'
#' @format A data frame with 2000 rows and 6 columns:
#' \describe{
#'   \item{U1}{minimum of event time T1 and censoring time C}
#'   \item{U2}{minimum of event time T2 and censoring time C}
#'   \item{C}{censoring time C, minimum of death time D and right-censoring time Ca}
#'   \item{delta_1}{indicator for observing T1}
#'   \item{delta_2}{indicator for observing T2}
#'   \item{delta_D}{indicator for observing D}
#' }
#' @source \url{http://www.diamondse.info/}
"test_data"
