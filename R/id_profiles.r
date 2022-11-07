#' Identify Choice Profiles
#'
#' Variables that comprise a composition may have cases where not all potential components of the composition
#' are observed.   We refer to the set of observed components of the composition as a "profile".  This function
#' identifies the existing profiles and gives a count for each one.
#' @param x A matrix or data frame that represents the elements of a composition. Components with values 0 or `NA`
#' will be considered as absent from that observation's composition.
#' @param sel A `tidy-select` statement that identifies the elements of the composition.  Note, the data must already be
#' in wide format.  The default is `everything()` which chooses all columns.
#' @param ... Not currently implemented.
#' @importFrom dplyr everything select group_by as_tibble ungroup mutate across case_when tally
#' @importFrom magrittr %>%
#' @importFrom glue glue
#' @export
id_profiles <- function(x, sel = everything(), ...){
  if(!inherits(x, "data.frame"))x <- as.data.frame(x)
  x <- x %>% 
    ungroup %>% 
    select(sel) %>% 
    mutate(across(everything(), ~case_when(.x > 0 ~ 1, TRUE ~ NA_real_)))
  
  prof <- data.frame(profile = apply(x, 1, function(z)glue("c({paste(which(!is.na(z)), collapse=',')})")))
  prof <- prof %>% group_by(profile) %>% tally() %>% mutate(col = 1:n())
  
  out <- tibble(component = colnames(x))
  out <- rbind(out, data.frame(component = "N"))
  profmat <- matrix("", nrow=nrow(out)-1, ncol=nrow(prof))
  for(i in 1:nrow(prof)){
    profmat[cbind(eval(parse(text=prof$profile[i])), prof$col[i])] <- "X"
  }  
  profmat <- rbind(profmat, prof$n)
  colnames(profmat) <- glue("Profile {1:ncol(profmat)}")
  cbind(out, profmat)  
}
