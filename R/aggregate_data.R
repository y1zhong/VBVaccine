#' @import parallel
#' @import doParallel
#' @import foreach
#' @export
prep_data = function(data, drug.case = NULL, drug.control = NULL,
                            covar_cont = NULL, covar_disc = NULL,
                            breaks = NULL, aggregate = FALSE){
  AEs = sort(unique(data$AE_NAME))
  if (aggregate) {
    # # comb <- expand.grid(unique(data$Age),unique(data$Male),unique(data$VAX_TYPE)) %>%
    # #   arrange(Var1, Var2, Var3)
    # # colnames(comb) <- c("Age", "Male", "VAX_TYPE")
    # nn = counts$n
    # A = matrix(0, nrow(counts), length(AEs))
    # colnames(A) = AEs
    # for(i in 1:nrow(counts)){
    #   data_sub = data[data$Age == counts$Age[i] &data$Male == counts$Male[i] &
    #            data$VAX_TYPE == counts$VAX_TYPE[i] ,]
    #   # data_sub = subset(data, Age==counts[i,1]& 
    #   #                     Male==counts[i,2]& VAX_TYPE==counts[i,3])
    #   A[i,] = sapply(1:length(AEs), function(j) sum(data_sub$AE_NAME == AEs[j]))
    #   # nn[i] = nrow(data_sub)
    # }
    # # X = cov
    # X = counts %>%
    #   select(-n) %>%
    #   as.matrix()
    # X = cbind(Intercept = 1, X)
    
    results = count_cases(data, drug.case = drug.case, drug.control = drug.control,
                covar_cont = covar_cont, covar_disc =covar_disc,
                breaks = breaks) %>%
       mutate(VAX_TYPE = as.numeric(VAX_TYPE == "DrugYes"))
    AEs = sort(unique(results$AE_NAME))
    A = matrix(results$AEYes, ncol = length(AEs))
    colnames(A) = AEs
    
    X = cbind(1, results[, c(covar_cont, covar_disc,"VAX_TYPE")])
  
  } else { 
    X <- data %>%
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
                  ifelse(data$AE_NAME == AEs[i], 1, 0)
                  }
    stopCluster(cl)
  }
  
  idx = which(colnames(X) == "VAX_TYPE")
  V = X[, idx]
  X = X[, -idx]
  
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
