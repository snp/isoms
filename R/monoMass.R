#' monoMass function calculates the monoisotopic mass of a molecule with a
#' composition \code{comp}
#' @param comp - named vector of molecule composition, e.g.
#' \code{c('C' = 4, 'H' = 8, 'N' = 1)}
#' @return monoisotopic mass of the molecule
#' @export
#'
#' @examples
#' monoMass(c('C' = 4, 'H' = 8, 'N' = 1)) # 70.06567
monoMass <- function(comp){
  mass <- 0
  data(isotopes)
  for(el in names(comp)){
    if(el %in% isotopes$element)
      mass <- mass + comp[el] * (isotopes %>% filter(element==el) %>% slice(1))$mass
    else
      warning(sprintf('Unknown element: %s', el))
  }
  unname(mass)
}
