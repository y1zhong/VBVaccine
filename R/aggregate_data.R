#' @import parallel
#' @import doParallel
#' @import foreach
#' @export
prep_data = function(df, vax.case = NULL, vax.control = NULL,
                     covar_cont = NULL,
                     covar_disc = NULL, aggregate = FALSE,
                     include_intercept = F){
  
  if (aggregate) {
    df <- df %>%
      filter(VAX_LABEL %in% c(vax.case, vax.control)) %>%
      mutate(VAX_LABEL = factor(as.numeric(VAX_LABEL == vax.case)))
    colnames(df)[1:3] = c("ID",  "VAX_LABEL", "AE_NAME")
    covar_cont_names = names(covar_cont)
    covar_disc_names = names(covar_disc)
    covars = c(covar_cont, covar_disc)
    p=length(covars)
    if(!is.null(covar_cont)){
      for(i in 1:length(covar_cont)){
        breaks_temp = covar_cont[[i]]
        var_temp = covar_cont_names[i]
        covar_group = df %>%
          dplyr::select(all_of(var_temp)) %>%
          unlist() %>%
          findInterval(breaks_temp) %>%
          as.factor()
        
        df = df %>%
          dplyr::select(-all_of(var_temp)) %>%
          tibble({{var_temp}} := covar_group)
      }
    }
    
    if(!is.null(covar_disc)){
      for(i in 1:length(covar_disc)){
        var_temp = covar_disc_names[i]
        df[, var_temp] = factor(df[, var_temp,drop=TRUE], levels = covar_disc[[i]], labels=c(1:length(covar_disc[[i]]))-1)
      }
    }
    
    X = list()
    for(i in 1:length(covar_cont)){
      X[[i]] = seq(from=0,to=length(covar_cont[[i]]))
    }
    
    Y = list()
    for(i in 1:length(covar_disc)){
      Y[[i]] = seq(from=0,to=length(covar_disc[[i]])-1)
    }
    V = list(V=c(0,1))
    X = c(V, X, Y)
    X = expand.grid(X)
    colnames(X) = c("VAX_LABEL", covar_cont_names, covar_disc_names)
    
    grp_cols = c("VAX_LABEL", "AE_NAME", covar_cont_names, covar_disc_names)
    data_count = df%>% dplyr::select(-ID) %>%
      group_by_at(grp_cols) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = AE_NAME, values_from = count) %>%
      replace(is.na(.), 0)
    
    data_merge = merge(X, data_count, all = T) %>%
      replace(is.na(.), 0)
    
    
    counts = df %>% count(!!sym("VAX_LABEL"), !!sym(covar_cont_names), !!sym(covar_disc_names), sort = TRUE)
    nn = merge(X, counts, all = T) %>%
      replace(is.na(.), 0)
    nn=nn$n
    
    A = data_merge
    V = A[, 1]
    X = as.matrix(A[, c(covar_cont_names, covar_disc_names)])
    A = A[, -which(colnames(A) %in%c("VAX_LABEL", covar_cont_names, covar_disc_names))]
    if(include_intercept) X = cbind(intercept=1,X)
    A = A %>% select(sort(names(.))) %>%
      as.matrix()

  } else { 
    AEs = sort(unique(df$AE_NAME))
    X <- df %>%
      distinct()%>%
      mutate(Intercept = 1) %>%
      select(Intercept, Age,Male, VAX_TYPE) %>%
      as.matrix()

    nn = rep(1, nrow(X))
    
    cores = floor(detectCores() / 1.5)
    cl = makeCluster(cores, setup_strategy = "sequential")
    registerDoParallel(cl)
    A = foreach(i = 1:length(AEs),
                .packages = c("tidyverse"),
                .combine=cbind
                ) %dopar% {
                  ifelse(df$AE_NAME == AEs[i], 1, 0)
                  }
    stopCluster(cl)
    # idx = which(colnames(X) == "VAX_LABLE")
    V = X[, 4]
    X = X[, -4]
  }
  
  return(data=list(A=A, nn=nn, X=X, V=V))
}

#' convert individual level VB data list to aggregate level
#' @export
aggregate_data = function(data) {
  A = data$A 
  X = data$X 
  V = data$V 
  nn = data$nn
  keys = sapply(1:length(V), function(i) {
    key = c(X[i,], V[i])
    return(paste0(key, collapse = '_'))
  })
  indi_sets = split(1:nrow(X), keys)
  names(indi_sets) = NULL
  new_n = length(indi_sets)
  new_nn = sapply(indi_sets, function(indi_set) {
    return(sum(nn[indi_set]))
  })
  new_A = matrix(0, new_n, ncol(A))
  new_X = matrix(0, new_n, ncol(X))
  new_V = rep(0, new_n)
  for (i in 1:new_n) {
    indi_set = indi_sets[[i]]
    new_A[i,] = colSums(A[indi_set, , drop = FALSE])
    new_X[i,] = X[indi_set[1],]
    new_V[i] = V[indi_set[1]]
  }
  colnames(new_A) = colnames(data$A)
  data$A = new_A
  data$X = new_X
  data$V = new_V
  data$nn = new_nn
  
  return(data)
}

#' @export
sequential_data = function(data, VB) {
  temp_vb = apply(cbind(VB$data$X, VB$data$V), 1, paste0, collapse="") 
  temp_data = apply(cbind(data$X, data$V), 1, paste0, collapse="")  
  
  new_A = matrix(0, nrow(VB$data$A), ncol(VB$data$A))
  new_nn = c()
  colnames(new_A) = colnames(VB$data$A)
  for(i in 1:length(temp_vb)){
    idx = which(temp_data == temp_vb[i])
    new_nn[i] = length(idx)
    new_A[i,] = colSums(data$A[idx, , drop = FALSE])
  }
  colnames(new_A) = colnames(data$A)
  return(data=list(A=new_A, nn=new_nn, X=VB$data$X, V=as.vector(VB$data$V)))
}

#' aggregate_data = function(data) {
#'   AEs = sort(unique(data$AE_NAME))
#'   counts <- count(data,Age,Male,VAX_TYPE)
#'   nn <- counts$n
#'   A = matrix(0, nrow(counts), length(AEs))
#'   colnames(A) = AEs
#'   for(i in 1:length(nn)){
#'     data_sub = subset(data, Age==counts$Age[i]& 
#'                         Male==counts$Male[i]& VAX_TYPE==counts$VAX_TYPE[i])
#'     A[i,] = sapply(1:length(AEs), function(j) sum(data_sub$AE_NAME == AEs[j]))
#'     
#'   }
#'   X = counts %>%
#'     select(-n) %>%
#'     as.matrix()
#'   X = cbind(Intercept = 1, X)
#'   
#'   V = X[, 4]
#'   X = X[, -4]
#'   return(data=list(A=A, nn=nn, X=X, V=V))
#' }
#' 
#' #' @export
#' individual_data = function(data) {
#'   AEs = sort(unique(data$AE_NAME))
#'   X <- data %>%
#'           distinct()%>%
#'           mutate(Intercept = 1) %>%
#'           select(Intercept, SEX, AGE, VAX_TYPE) %>%
#'           as.matrix()
#'   V <- X[, 4]
#'   X <- X[, -4]
#'   nn = rep(1, nrow(X))
#'   
#'   cores = floor(detectCores() / 1.5)
#'   cl = makeCluster(cores, setup_strategy = "sequential")
#'   registerDoParallel(cl)
#'   A = foreach(i = 1:length(AEs),
#'               .packages = c("tidyverse"),
#'               .combine=cbind
#'               ) %dopar% {
#'                 ifelse(data$AE_NAME == AEs[i], 1, 0)
#'                 }
#'   stopCluster(cl)
#' 
#' 
#'   return(data=list(A=A, nn=nn, X=X, V=V))
#' }

#' @export
expand_data = function(data) {
  AEs = sort(unique(data$AE_NAME))
  counts <- count(data,Age,Male,V)
  nn <- counts$n
  A = matrix(0, nrow(counts), length(AEs))
  colnames(A) = AEs
  for(i in 1:length(nn)){
    data_sub = subset(data, Age==counts$Age[i]& 
                        Male==counts$Male[i]& V==counts$V[i])
    A[i,] = sapply(1:length(AEs), function(j) sum(data_sub$AE_NAME == AEs[j]))
    
  }
  X = counts %>%
    select(-n) %>%
    as.matrix()
  X = cbind(Intercept = 1, X)
  
  V = X[, 4]
  X = X[, -4]
  return(data=list(A=A, nn=nn, X=X, V=V))
}


#' #' @export
#' prep_data <- function(data, initials="glm", 
#'                       mu_alpha=NULL,
#'                       a_alpha = 0.01, 
#'                       b_alpha  = 0.01,
#'                       include_intercept=FALSE){
#'   if (include_intercept){
#'     data$X = cbind(1, data$X)
#'   }
#'   
#'   J = ncol(data$A)
#'   p = ncol(data$X)
#'   if(is.null(mu_alpha)){
#'     mu_alpha = matrix(0, ncol = J, nrow=p)
#'   } 
#'   VB <- list()
#'   glm_coef <- get_glm_coefs(data$A, data$X, data$V,data$nn)
#'   VB[["post_mean"]] = list(beta = matrix(glm_coef$beta_glm, ncol = 1),
#'                            alpha = t(glm_coef$Alpha_glm),
#'                            inv_sigma2_beta = matrix(1, ncol=1, nrow=J),
#'                            inv_sigma2_alpha = matrix(1, nrow = p, ncol=J))
#'   
#'   VB[["data"]] = data
#'   
#'   VB[["hyper"]] = list(a_alpha = a_alpha,
#'                        b_alpha = b_alpha,
#'                        mu_alpha = mu_alpha)
#'   VB[["trace"]] = list(iters=c(), ELBO = c())
#'   
#'   return(VB)
#' }
