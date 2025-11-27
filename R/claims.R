#' Personal Injury Insurance Claims
#'
#' A dataset containing personal injury insurance claims made in Australia
#' from January 1998 to June 1999.
#'
#' @format A data frame with 769 rows and 4 variables:
#' \describe{
#' \item{amount}{amount of paid money by an insurance policy or claim size, in Australian dollars}
#' \item{legrep}{with legal representation or not}
#' \item{month}{month of occurrence of the accident}
#' \item{optime}{operational time in percentage with range 0 to 100}
#' }
#' @source de Jong P, Heller GZ (2008) Generalized linear models for insurance data.
#' Cambridge University Press, Cambridge
"claims"
#claims <- read.table(file = "insurance.txt", header = TRUE)
#usethis::use_data(claims,overwrite = TRUE)
