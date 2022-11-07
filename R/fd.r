#' First Differences for Compositional Regression Models
#' 
#' Thanks to the non-linear transformation of the composition, the effects are on each component of the composition
#' are not as straightforward to assess as the linear model may suggest.  We use first differences to evaluate 
#' the effect of the variables in the model.  
#' 
#' @param models A list of models estimated by `compreg`. 
#' @param variable A string identifying the name of the variable whose first difference will be calculated. 
#' @param diff A string identifying the kind of first difference desired - `"sd"` for standard deviation change 
#' around observed values, `"range"` for a maximal change for the variable or `"unit"` for a one-unit change around
#' observed values. 
#' @param values An optional vector of values to use for the first difference.  If `values` is specified, it overrides
#' `diff`.  
#' @param R Number of bootstrap samples to use. 
#' @param add_zero Should you add in values of zero for components that are not part of the composition (but appear in other compositions)? 
#' @param ... Not implemented
#' 
#' @details This function does a parametric boostrap of the coefficient matrix to estimate the variability of the 
#' predicted compositions.  This permits the calculation of confidence intervals for first differences.  
#' 
#' @importFrom MASS mvrnorm
#' @importFrom abind abind
#' @importFrom tidyr pivot_longer
#' @export
fd <- function(models, variable, diff, values=NULL, R=2500, add_zero=TRUE, ...){
  vbl <- data %>% select(all_of(variable)) %>% pull  
  UseMethod("fd", vbl)
}

#' @method fd factor
#' @export
fd.factor <- function(models, variable,values, R, add_zero, ...){
  itfun <- eval(parse(text=paste0(attr(models[[wmx]]$model, "transform"), "Inv")))
  if(!is.null(values) & length(values) != 2)stop("The values argument must be length 2.\n")
  if(is.null(values)){
  n <- sapply(models, function(x)nrow(x$data))
  unvbl <- lapply(models, function(x)unique(x$data[[variable]]))
  levs <- unvbl[[1]]
  if(length(models) > 1){
    for(j in 2:length(models)){
      levs <- intersect(levs, unvbl[[j]])
    }
  }
  if(length(levs) < 2)stop("No 2 levels are common across all subset datasets.\n")
  ll <- sapply(unvbl, length)
  if(!all(ll == ll[1]) & length(ll) > 1)warning("Not all levels of the factor variable are present in all regression models.  Proceeding with levels common across all models.\n")
  wmx <- which.max(n)
  X <- model.matrix(models[[wmx]]$model)
  h <- diag(X %*% solve(t(X) %*% X) %*% t(X))
  d <- as.list(models[[wmx]]$data[which.min(h)[1], ])
  d[[variable]] <- levs
  d <- do.call(tibble, d)
  fit <- apply(itfun(predict(models[[wmx]]$model, newdata=d)), 2, as.numeric)
  combs <- combn(1:nrow(fit), 2)
  D <- matrix(0, nrow=nrow(fit), ncol=ncol(combs))
  D[cbind(combs[1,], 1:ncol(D))] <- -1
  D[cbind(combs[2,], 1:ncol(D))] <- 1
  wmx_pair <- which.max(colSums(abs(t(fit) %*% D)))[1]
  values <- levs[c(combs[,wmx_pair])]  
  }
  preds <- list()
  for(x in 1:length(models)){
    coefs <- coef(models[[x]]$model)
    if(!inherits(coefs, "matrix"))coefs <- matrix(coefs, ncol=1)
    B <- MASS::mvrnorm(R, c(coef(models[[x]]$model)), vcov(models[[x]]$model))
    B <- array(t(B), dim=c(nrow(coefs), 
                             ncol(coefs), 
                             R))
    D1 <- D2 <- models[[x]]$data
    D1[[variable]] <- values[1]
    D2[[variable]] <- values[2]
    X1 <- model.matrix(models[[x]]$model, data=D1)
    X2 <- model.matrix(models[[x]]$model, data=D2)
    preds1 <- lapply(1:R, function(i)apply(itfun(X1 %*% B[,,i]), 2, as.numeric))
    preds2 <- lapply(1:R, function(i)apply(itfun(X2 %*% B[,,i]), 2, as.numeric))
    preds1$along <- 3
    preds2$along <- 3
    preds1 <- do.call(abind, preds1)
    preds2 <- do.call(abind, preds2)
    preds[[x]] <- list(p1 = preds1, p2 = preds2, diff = preds2-preds1)    
  }
  names(preds) <- names(models)
  long_diffs <- list()
  for(i in 1:length(preds)){
    out <- apply(preds[[i]]$diff, c(2,3), function(z)c(mean(z)))
    rownames(out) <- strsplit(names(preds)[i], "")[[1]]
    colnames(out) <- 1:ncol(out)
    out <- t(out)
    long_diffs[[i]] <- as_tibble(out, rownames="sim") %>% 
      mutate(profile = names(preds)[i], 
             n=dim(preds[[i]]$diff)[1]) %>% 
      pivot_longer(-c("sim", "profile", "n"), names_to="letter", values_to="delta")
      
  }
  names(long_diffs) <- names(models)
  prof_diffs <- lapply(long_diffs, function(x){
    x %>% 
      group_by(profile, n, letter) %>% 
      summarise(diff = mean(delta), lwr = unname(quantile(delta, .025)), upr=unname(quantile(delta, .975)))  %>% 
      left_join(attr(models, "legend")) %>% 
      select(profile, n, variable, diff, lwr, upr)
  })
  long_diffs$.id <- "pnum"
  ag_diffs <- do.call(bind_rows, long_diffs) %>% 
    mutate(delta = delta*(n/total)) %>% 
    select(pnum, sim, letter, delta)
  if(add_zero){
    ag_diffs <- ag_diffs %>% complete(pnum, sim, letter, fill=list(delta=0))
  }
  ag_diffs <- ag_diffs %>% 
    group_by(letter) %>%
    summarise(diff = mean(delta), lwr = unname(quantile(delta, .025)), upr=unname(quantile(delta, .975)))  %>% 
    left_join(attr(models, "legend")) %>% 
    select(variable, diff, lwr, upr)
  return(list(profiles = prof_diffs, aggregate = ag_diffs))
}

fd.numeric <- function(models, variable, diff, values, R, add_zero, ...){
  itfun <- eval(parse(text=paste0(attr(models[[1]]$model, "transform"), "Inv")))
  if(!is.null(values) & length(values) != 2)stop("The values argument must be length 2.\n")
  if(is.null(values)){
    if(diff == "unit")values <- c(-.5, .5)
    if(diff == "sd"){
      values <- do.call(rbind, lapply(models, function(x)x$data)) %>% 
        select(all_of(variable)) %>% 
        pull %>% 
        sd*c(-.5,.5)
    }
    if(diff == "range"){
      values <- do.call(rbind, lapply(models, function(x)x$data)) %>% 
        select(all_of(variable)) %>% 
        pull %>% 
        range
      
    }
  }
  preds <- list()
  for(x in 1:length(models)){
    coefs <- coef(models[[x]]$model)
    if(!inherits(coefs, "matrix"))coefs <- matrix(coefs, ncol=1)
    B <- MASS::mvrnorm(R, c(coef(models[[x]]$model)), vcov(models[[x]]$model))
    B <- array(t(B), dim=c(nrow(coefs), 
                           ncol(coefs), 
                           R))
    D1 <- D2 <- models[[x]]$data
    D1[[variable]] <- D1[[variable]] + values[1]
    D2[[variable]] <- D2[[variable]] + values[2]
    if(diff == "range"){
      D1[[variable]] <- values[1]
      D2[[variable]] <- values[2]
    }
    X1 <- model.matrix(models[[x]]$model, data=D1)
    X2 <- model.matrix(models[[x]]$model, data=D2)
    preds1 <- lapply(1:R, function(i)apply(itfun(X1 %*% B[,,i]), 2, as.numeric))
    preds2 <- lapply(1:R, function(i)apply(itfun(X2 %*% B[,,i]), 2, as.numeric))
    preds1$along <- 3
    preds2$along <- 3
    preds1 <- do.call(abind, preds1)
    preds2 <- do.call(abind, preds2)
    preds[[x]] <- list(p1 = preds1, p2 = preds2, diff = preds2-preds1)    
  }
  names(preds) <- names(models)
  long_diffs <- list()
  for(i in 1:length(preds)){
    out <- apply(preds[[i]]$diff, c(2,3), function(z)c(mean(z)))
    rownames(out) <- strsplit(names(preds)[i], "")[[1]]
    colnames(out) <- 1:ncol(out)
    out <- t(out)
    long_diffs[[i]] <- as_tibble(out, rownames="sim") %>% 
      mutate(profile = names(preds)[i], 
             n=dim(preds[[i]]$diff)[1]) %>% 
      pivot_longer(-c("sim", "profile", "n"), names_to="letter", values_to="delta")
    
  }
  names(long_diffs) <- names(models)
  total <- sum(sapply(long_diffs, function(x)x$n[1]))
  prof_diffs <- lapply(long_diffs, function(x){
    x %>% 
      group_by(profile, n, letter) %>% 
      summarise(diff = mean(delta), lwr = unname(quantile(delta, .025)), upr=unname(quantile(delta, .975)))  %>% 
      left_join(attr(models, "legend")) %>% 
      select(profile, n, variable, diff, lwr, upr)
  })
  long_diffs$.id <- "pnum"
  ag_diffs <- do.call(bind_rows, long_diffs) %>% 
    mutate(delta = delta*(n/total)) %>% 
    select(pnum, sim, letter, delta)
  if(add_zero){
    ag_diffs <- ag_diffs %>% complete(pnum, sim, letter, fill=list(delta=0))
  }
  ag_diffs <- ag_diffs %>% 
    group_by(letter) %>%
    summarise(diff = mean(delta), lwr = unname(quantile(delta, .025)), upr=unname(quantile(delta, .975)))  %>% 
    left_join(attr(models, "legend")) %>% 
    select(variable, diff, lwr, upr)
  return(list(profiles = prof_diffs, aggregate = ag_diffs))
}



