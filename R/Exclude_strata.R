#' Functions to define strata to be used in the mackerel distribution work
#'
#' These strata are based on pre-defined strata or manually defined through expert knowledge:
#'
#' @param data A data frame that includes all strata names available in the study area
#' @param exclude A vector of stratum name to be excluded from the modelled area 
#' @param ... other inputs taken from the general environment
#'
#' @return
#' @export
#'
#' @examples
exclude_area <- function(data, exclude) {
  excluded <- data %>% filter(ID %in% exclude) %>% 
    st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>% 
    group_by(ID) %>% 
    summarize(geometry = st_combine(geometry)) %>% st_cast("POLYGON") 
  return(excluded)
}
