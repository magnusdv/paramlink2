#' Dominant linkage analysis example
#'
#' A dataset with SNP genotypes for 14 members of a pedigree affected with a
#' dominant disorder.
#'
#' @format A list with 3 elements:
#'
#' * `ped`: A [pedtools::ped] object describing a pedigree with 19 individuals,
#'   including genotypes for 14 members at 248 markers on Chromosome 1.
#'
#' * `aff`: A vector indicating the affected pedigree members.
#'
#' * `map`: A data frame with 3 columns (`chrom`, `marker`, `cm`) describing the
#'   centiMorgan positions of the markers.
#'
#' @examples
#'
#' data(dominant1)
#' ped = dominant1$ped
#' aff = dominant1$aff
#' model = diseaseModel("AD")
#'
#' # Compute singlepoint LODs
#' lods = lod(ped, aff = aff, model = model)
#'
#' # LOD score graph
#' plot(lods)
#'
#'
"dominant1"
