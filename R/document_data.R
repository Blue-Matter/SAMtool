#' North Atlantic Swordfish dataset
#'
#' An S4 object containing catch and index time series for North Atlantic swordfish.
#'
#' @format An object of class \linkS4class{Data}.
#' @source ASPIC Software at https://www.mhprager.com/aspic.html
#'
#' @examples
#' data(swordfish)
#'
"swordfish"


#' A two-fleet Albacore operating model
#'
#' A generic operating model for an albacore stock with two fishing fleets. The first fleet has
#' dome-shaped selectivity (similar to a baitboat fleet) while the second fleet exhibits logistic
#' selectivity (such as a longline fleet). With the \code{CatchFrac} slot, we generate a 30%-70% catch
#' ratio between the baitboat-longline fleets in the most recent historical year.
#'
#' @format An object of class \linkS4class{MOM}.
#'
#' @examples
#' \donttest{
#' ## Plot historical effort and selectivity between the 2 fleets
#' plot(Albacore_TwoFleet)
#'
#' ## Generate data (e.g., catch, length comps) from the fleets
#' Hist <- multiMSE(Albacore_TwoFleet, Hist = TRUE)
#' DataList <- Hist$Data
#' }
"Albacore_TwoFleet"
