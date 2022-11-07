#' Run a Compositional Regression on All Profiles
#' 
#' Runs a compositional regression model on all profiles with at least `n_min` observations
#' and saves the models in a way that is amenable to the post-model estimation routines in the 
#' package. 
#' @param form A right-hand sided formula of the independent variables in the regression model. 
#' @param composition A `tidy-select` statement that identifies the components of the composition in the data.  Note, 
#' the data must be in wide-format already. 
#' @param data A data frame that contains variables in `formula` and `composition`.  
#' @param count Logical indicating whether the data represent counts or not. 
#' @param total Logical indicating whether the components of the composition convey information about the total (`TRUE`) or 
#' their sum is irrelevant or uninformative (`FALSE`). 
#' @param relative Logical indicating whether the components of the composition are ratio level (`TRUE`) or interval level (`FALSE`).  
#' @param n_min Number of observations required to estimate a model for the profile. 
#' @param transform Transformation used for the composition.  Options are `"alr"`, `"clr"`, `"ilr"`, `"apt"`, `"cpt"` and `"ipt"`.  
#' For more information, see `help(package="compositions")`. 
#' @param ... Other arguments passed down to `lm()`. 
#' 
#' @details For more information on the `count`, `total` and `relative` arguments, see `system.file("doc/compositions_v2.html", package="compositions")`
#' @importFrom compositions rcomp acomp rplus aplus ccomp alr alrInv clr clrInv ilr ilrInv apt aptInv cpt cptInv ipt iptInv
#' @importFrom stats lm
#' @importFrom dplyr tibble
#' @importFrom DAMisc unformulate
#' @export
compreg <- function(form, 
                    composition, 
                    data, 
                    count, 
                    total, 
                    relative, 
                    n_min, 
                    transform = c("alr", "clr", "ilr", "apt", "cpt", "ipt"), 
                    ...){
  `%x%` <- function(X,Y)paste(ifelse(is.na(Y),"", X), collapse="")
  if(!count & !total & !relative)compfun <- rcomp
  if(!count & !total & relative)compfun <- acomp
  if(!count & total & !relative)compfun <- rplus
  if(!count & total & relative)compfun <- aplus
  if(count & total & relative)compfun <- ccomp
  if((count & !total & !relative) | (count & !total & relative) | (count & total & !relative)){
    stop("Not all combinations of count, total and relative are permitted, see system.file('doc/compositions_v2.html', package='compositions') for more information.\n")
  }  
  trans <- match.arg(transform)
  tfun <- eval(parse(text = trans))
  itfun <- eval(parse(text = glue("{trans}Inv")))
  vrs <- unformulate(form)$vars
  compvars <- data %>% ungroup %>% select(composition) %>% names()
  
  tmp <- data %>% 
    ungroup %>% 
    select(c(all_of(vrs), all_of(compvars))) 
  comp <- tmp %>% 
    select(all_of(compvars)) 
  legend <- tibble(variable = colnames(comp), 
                   letter = LETTERS[1:ncol(comp)])
  
  tmp$prof <- comp %>% 
    mutate(across(everything(), ~case_when(.x > 0 ~ 1, TRUE ~ NA_real_))) %>% 
    apply(., 1, function(comp)LETTERS[1:length(comp)] %x% comp)
  tmp <- tmp %>% filter(nchar(prof) > 1) %>% group_by(prof) %>% filter(n() > n_min) %>% ungroup
  sp <- split(tmp, tmp$prof)
  mods <- lapply(sp, function(x){
    comp <- x %>% select(all_of(compvars))
    ana <- apply(comp, 2, function(x)all(is.na(x)))
    if(any(ana)){
      comp <- comp[,-which(ana)]
    }
    form <- reformulate(vrs, 
                        response = paste0("apply(tfun(cbind(", paste(colnames(comp), collapse=","), ")), 2, as.numeric)"))
    mod <- lm(form, data=x, ...)        
    attr(mod, "transform") <- trans
    list(model = mod, data= x)
  })
  attr(mods, "legend") <- legend
  return(mods)
}
