
VB_Vacc_r <- function(A, X, V, nn,
                    init_beta, init_alpha, 
                    init_sigma_beta, init_sigma_alpha,
                    mu_alpha=NULL,
                    a_alpha = 0.1, b_alpha=0.1,
                    include_intercept=TRUE,
                    max_iter = 1000, paras_diff_tol = 1e-6,
                    ELBO_stop = 1, ELBO_diff_tol = 1e-6){
  
  if(include_intercept){
    X = cbind(1,X)
  }
  
  
  p = ncol(X)
  n = nrow(X)
  J = ncol(A)
  ELBO = -Inf
  
  if(is.null(mu_alpha)) {
    mu_alpha = matrix(0, nrow = p, ncol = J)
  }
  
  
  ###initialize
  a_alpha = a_alpha
  b_alpha = b_alpha
  
  
  E_beta = init_beta
  E_inv_sigma2_beta = rep(1, J)
  Var_beta = E_inv_sigma2_beta
  E_beta_sq = Var_beta + E_beta^2
  E_inv_a_beta = 0.5
  
  
  E_alpha = init_alpha
  E_inv_sigma2_alpha = matrix(1, nrow = p, ncol = J)
  Var_alpha = matrix(E_inv_sigma2_alpha, nrow=p, ncol=J)
  E_alpha_sq = Var_alpha + E_alpha^2
  
  
  
  E_phi = X%*%E_alpha + tcrossprod(V, E_beta)
  E_omega = sweep(tanh(E_phi/2)/(2*E_phi), 1,nn, "*")
  
  
  ELBO_all = c()
  iter = 1
  while(iter < max_iter){
    #cat(iter," ")
    #update_E_beta()
    
    #update_E_phi()
    E_phi <- X %*% E_alpha + tcrossprod(V, E_beta)
    
    #update_E_omega()
    E_omega <- sweep(tanh(E_phi/2)/(2*E_phi), 1,nn, "*")
    #E_omega[is.nan(E_omega)] <- 1000
    
    #update_E_inv_sigma2_beta()
    E_inv_sigma2_beta <- 1/(E_beta_sq/2 + E_inv_a_beta)
    
    
    #update_E_inv_sigma2_alpha()
    E_inv_sigma2_alpha <- (a_alpha + 1/2) / (b_alpha + (E_alpha_sq + mu_alpha^2 - 2*E_alpha*mu_alpha)/2)
    #E_inv_sigma2_alpha <- (a_alpha + 1/2) / (b_alpha + Var_alpha/2)
    
    #update_E_inv_a_beta()
    E_inv_a_beta <- ((J+1)/2)/(sum(E_inv_sigma2_beta)+1)
    
    Var_beta <- sapply(1:J, function(j) 1/(sum(V * E_omega[, j]) + E_inv_sigma2_beta[j]))
    E_beta <- sapply(1:J, function(j) Var_beta[j] * V %*% (A[, j] - nn/2 - E_omega[, j] * (X %*% E_alpha[, j])))
    E_beta_sq <- Var_beta + E_beta^2
    
    #cat("E_beta: ",E_beta, " ")
    
    #update_E_alpha()
    
    for(j in 1:J){
      cov_mat_j = crossprod(sweep(X, 1, E_omega[,j], "*"), X)
      diag(cov_mat_j) = diag(cov_mat_j) + E_inv_sigma2_alpha[,j]
      cov_mat_j = solve(cov_mat_j);
      Var_alpha[, j] <- diag(cov_mat_j)
      
      E_alpha[,j] <- tcrossprod(cov_mat_j,((A[,j] - nn/2 - E_omega[,j]*V*E_beta[j]) %*% X + 
                                                      E_inv_sigma2_alpha[,j]*mu_alpha[,j]))
    }
    E_alpha_sq <- sapply(1:J, function(j) Var_alpha[,j] + E_alpha[,j]^2)
    
    
    
    ELBO_old = ELBO
    
    #update_ELBO()
    ELBO <- sum(sweep(A, 1, nn/2, `-`) * E_phi)
    ELBO <- ELBO + 0.5 * sum(log(Var_beta))
    ELBO <- ELBO - sum(log(0.5 * E_beta_sq + E_inv_a_beta))
    ELBO <- ELBO + E_inv_a_beta * sum(E_inv_sigma2_beta)
    
    temp_var_alpha <- sapply(1:J, function(j){
      temp_var_j = crossprod(sweep(X, 1, E_omega[,j], "*"), X)
      diag(temp_var_j) = diag(temp_var_j) + E_inv_sigma2_alpha[,j]
      ldet(temp_var_j)
    })
    
    # temp_cov_alpha <- sapply(1:J, function(j){
    #   temp_var_j = crossprod(sweep(X, 1, E_omega[,j], "*"), X)
    #   diag(temp_var_j) = diag(temp_var_j) + E_inv_sigma2_alpha[,j]
    #   cov_mat_j = solve(temp_var_j)
    #   trace_j=sum(diag(matrix(diag(E_inv_sigma2_alpha[,j]),ncol=p) %*% cov_mat_j))
    #   tcrossprod(E_alpha[,j]-mu_alpha[,j], matrix(diag(E_inv_sigma2_alpha[,j]),ncol=p)) %*% (E_alpha[,j]-mu_alpha[,j])+trace_j
    # })
    
    ELBO <- ELBO - 0.5*sum(temp_var_alpha)
    # ELBO <- ELBO - 0.5*sum(temp_cov_alpha)
    ELBO <- ELBO - sum((a_alpha+1/2)*log(b_alpha+0.5*(E_alpha_sq + mu_alpha^2 - 2*E_alpha*mu_alpha)))
    #ELBO <- ELBO - sum((a_alpha+1/2)*log(b_alpha+0.5*Var_alpha))
    ELBO <- ELBO + sum(Var_alpha * E_inv_sigma2_alpha)/2
    ELBO <- ELBO - (J+1)/2 * log(sum(E_inv_sigma2_beta)+1)
    ELBO <- ELBO - sum(log(cosh(E_phi/2)))
    if(iter %% 100 == 0) cat("iter: ", iter, ", ELBO: ", ELBO, "\n")
    #save profile
    if(abs(ELBO - ELBO_old) < ELBO_diff_tol) break
    ELBO_all = c(ELBO_all, ELBO)
    iter = iter + 1
  }
  
  output = list("post_mean" = list("beta" = E_beta,
                                   "alpha" = E_alpha,
                                   "inv_sigma2_beta" = E_inv_sigma2_beta,
                                   "inv_sigma2_alpha" = E_inv_sigma2_alpha),
                "iter" = 1:iter,
                "trace" = ELBO_all)
  return(output)
}



get_rv_coefs = function(A,X,V,
                        cutoff = 10) {
  # run glm models
  Alpha_glm = matrix(0, ncol(A), ncol(X) + 1)
  #suppressWarnings({
  for (j in 1:ncol(A)) {
    data <- data.frame(Aj=as.factor(A[,j]),X[,-1],V) 
    model = bayesreg(Aj~.,data, model = "logistic", prior = "horseshoe", n.samples = 1e3)
    alpha_glm = rowMeans(rbind(model$beta0, model$beta))
    # alpha_glm = pmin(alpha_glm,cutoff)
    # alpha_glm = pmax(alpha_glm,-cutoff)
    Alpha_glm[j,] =  alpha_glm
  }
  #})
  
  rownames(Alpha_glm) = NULL
  beta_glm = Alpha_glm[, ncol(Alpha_glm)]
  Alpha_glm = Alpha_glm[,-ncol(Alpha_glm)]
  colnames(Alpha_glm) = colnames(X)
  
  return(list(beta_glm = beta_glm,
              Alpha_glm = Alpha_glm))
}

# Compute the log-determinant of a matrix
ldet <- function(X) {
  if(!is.matrix(X)) return(log(X))
  determinant(X,logarithm = TRUE)$modulus
}


get_glm_coefs = function(A,X,V,
                         cutoff = 10) {
  # run glm models
  Alpha_glm = matrix(0, ncol(A), ncol(X) + 1)
  suppressWarnings({
    for (j in 1:ncol(A)) {
      Aj = A[, j]
      model = glm(Aj  ~ 0 + X + V,
                  family = binomial())
      alpha_glm = coef(model)
      alpha_glm = pmin(alpha_glm,cutoff)
      alpha_glm = pmax(alpha_glm,-cutoff)
      Alpha_glm[j,] =  alpha_glm
    }
  })
  
  rownames(Alpha_glm) = NULL
  beta_glm = Alpha_glm[, ncol(Alpha_glm)]
  Alpha_glm = Alpha_glm[,-ncol(Alpha_glm)]
  colnames(Alpha_glm) = colnames(X)
  
  return(list(beta_glm = beta_glm,
              Alpha_glm = Alpha_glm))
}
# VB_Vacc <- function(A, X, V, 
#                     init_beta, init_alpha, 
#                     init_sigma_beta, init_sigma_alpha,
#                     a_alpha = 0.1, b_alpha=0.1,
#                     include_intercept=TRUE,
#                     max_iter = 1000, paras_diff_tol = 1e-6,
#                     ELBO_stop = 1, ELBO_diff_tol = 1e-6){
#   
#   if(include_intercept){
#     X = cbind(1,X)
#   }
#   
#   p = ncol(X)
#   n = nrow(X)
#   J = ncol(A)
#   ELBO = -Inf
#   
#   ###initialize
#   a_alpha = a_alpha
#   b_alpha = b_alpha
#   
#   
#   E_beta = init_beta
#   E_inv_sigma2_beta = rep(1, J)
#   Var_beta = E_inv_sigma2_beta
#   E_beta_sq = Var_beta + E_beta^2
#   E_inv_a_beta = 0.5
#   
#   
#   E_alpha = init_alpha
#   E_inv_sigma2_alpha = rep(1,p)
#   Var_alpha = matrix(E_inv_sigma2_alpha, nrow=p, ncol=J)
#   E_alpha_sq = Var_alpha + E_alpha^2
#  
#   
# 
#   E_phi = X%*%E_alpha + tcrossprod(V, E_beta)
#   E_omega = matrix(0.1, nrow=n, ncol=J)
#   
#   
#   ELBO_all = c()
#   iter = 1
#   while(iter < max_iter){
#     #cat(iter," ")
#     #update_E_beta()
#     
# 
#     Var_beta <- sapply(1:J, function(j) 1/(sum(V * E_omega[, j]) + E_inv_sigma2_beta[j]))
#     E_beta <- sapply(1:J, function(j) Var_beta[j] * V %*% (A[, j] - 0.5 - E_omega[, j] * (X %*% E_alpha[, j])))
#     E_beta_sq <- Var_beta + E_beta^2
#     
#     #cat("E_beta: ",E_beta, " ")
#     
#     #update_E_alpha()
#    
#     for(j in 1:J){
#       cov_mat_j = crossprod(sweep(X, 1, E_omega[,j], "*"), X)
#       diag(cov_mat_j) = diag(cov_mat_j) + E_inv_sigma2_alpha
#       cov_mat_j = solve(cov_mat_j);
#       Var_alpha[, j] <- diag(cov_mat_j)
# 
#       E_alpha[,j] <-  tcrossprod(cov_mat_j, X) %*% (A[,j] - 0.5 - E_omega[,j]*V*E_beta[j])
#     }
#     E_alpha_sq <- sapply(1:J, function(j) Var_alpha[,j] + E_alpha[,j]^2)
#     
#     
#     #update_E_inv_sigma2_beta()
#     E_inv_sigma2_beta <- 1/(E_beta_sq/2 + E_inv_a_beta)
#     
#     
#     #update_E_inv_sigma2_alpha()
#     E_inv_sigma2_alpha <- sapply(1:p,  function(l) (a_alpha + J/2) / (b_alpha + sum(E_alpha_sq[l,])/2))
#     
#     #update_E_inv_a_beta()
#     E_inv_a_beta <- ((J+1)/2)/(sum(E_inv_sigma2_beta)+1)
#     
#     #update_E_phi()
#     E_phi <- X %*% E_alpha + tcrossprod(V, E_beta)
#     
#     #update_E_omega()
#     E_omega <- tanh(E_phi/2)/(2*E_phi)
#     #E_omega[is.nan(E_omega)] <- 0.25
#     
#     ELBO_old = ELBO
#     
#     #update_ELBO()
#     # ELBO <- sum(log(Var_beta))/2;
#     # ELBO <- ELBO + ldet(diag(E_inv_sigma2_alpha))*J/2;
#     # ELBO <- ELBO - sum(log(E_beta_sq/2 + E_inv_a_beta));
#     # ELBO <- ELBO + E_inv_a_beta*(sum(E_inv_sigma2_beta));
#     # ELBO <- ELBO - (p*b_alpha +  sum(Var_alpha)/2)*log(a_alpha+J/2);
#     # ELBO <- ELBO - (J+1)*log(sum(E_inv_sigma2_beta)+1)/2;
#     # ELBO <- ELBO + sum((A - 0.5)*E_phi);
#     ELBO <- sum((A - 0.5) * E_phi)
#     ELBO <- ELBO + 0.5 * sum(log(Var_beta))
#     ELBO <- ELBO - sum(log(0.5 * E_beta_sq + E_inv_a_beta))
#     ELBO <- ELBO + E_inv_a_beta * sum(E_inv_sigma2_beta)
#     
#     temp_var_alpha <- sapply(1:J, function(j){
#       temp_var_j = crossprod(sweep(X, 1, E_omega[,j], "*"), X)
#       diag(temp_var_j) = diag(temp_var_j) + E_inv_sigma2_alpha
#       ldet(temp_var_j)
#     })
#     ELBO <- ELBO - 0.5*sum(temp_var_alpha)
#     ELBO <- ELBO - sum((a_alpha+J/2)*log(b_alpha+0.5*rowSums(E_alpha_sq)))
#     ELBO <- ELBO - (J+1)/2 * log(sum(E_inv_sigma2_beta)+1)
#     ELBO <- ELBO - sum(log(cosh(E_phi/2)))
#     if(iter %% 100 == 0) cat("iter: ", iter, ", ELBO: ", ELBO, "\n")
#     #save profile
#     if(abs(ELBO - ELBO_old) < ELBO_diff_tol) break
#     ELBO_all = c(ELBO_all, ELBO)
#     iter = iter + 1
#   }
#   
#   output = list("post_mean" = list("beta" = E_beta,
#                                    "alpha" = E_alpha,
#                                    "inv_sigma2_beta" = E_inv_sigma2_beta,
#                                    "inv_sigma2_alpha" = E_inv_sigma2_alpha),
#                  "iter" = 1:iter,
#                  "trace" = ELBO_all)
#   return(output)
# }
# 
# 
# 
# get_rv_coefs = function(A,X,V,
#                          cutoff = 10) {
#   # run glm models
#   Alpha_glm = matrix(0, ncol(A), ncol(X) + 1)
#   #suppressWarnings({
#     for (j in 1:ncol(A)) {
#       data <- data.frame(Aj=as.factor(A[,j]),X[,-1],V) 
#       model = bayesreg(Aj~.,data, model = "logistic", prior = "horseshoe", n.samples = 1e3)
#       alpha_glm = rowMeans(rbind(model$beta0, model$beta))
#       # alpha_glm = pmin(alpha_glm,cutoff)
#       # alpha_glm = pmax(alpha_glm,-cutoff)
#       Alpha_glm[j,] =  alpha_glm
#     }
#   #})
#   
#   rownames(Alpha_glm) = NULL
#   beta_glm = Alpha_glm[, ncol(Alpha_glm)]
#   Alpha_glm = Alpha_glm[,-ncol(Alpha_glm)]
#   colnames(Alpha_glm) = colnames(X)
#   
#   return(list(beta_glm = beta_glm,
#               Alpha_glm = Alpha_glm))
# }
# # Compute the log-determinant of a matrix
# ldet <- function(X) {
#   if(!is.matrix(X)) return(log(X))
#   determinant(X,logarithm = TRUE)$modulus
# }
# 
# 
# get_glm_coefs = function(A,X,V,
#                          cutoff = 10) {
#   # run glm models
#   Alpha_glm = matrix(0, ncol(A), ncol(X) + 1)
#   suppressWarnings({
#     for (j in 1:ncol(A)) {
#       Aj = A[, j]
#       model = glm(Aj  ~ 0 + X + V,
#                   family = binomial())
#       alpha_glm = coef(model)
#       alpha_glm = pmin(alpha_glm,cutoff)
#       alpha_glm = pmax(alpha_glm,-cutoff)
#       Alpha_glm[j,] =  alpha_glm
#     }
#   })
#   
#   rownames(Alpha_glm) = NULL
#   beta_glm = Alpha_glm[, ncol(Alpha_glm)]
#   Alpha_glm = Alpha_glm[,-ncol(Alpha_glm)]
#   colnames(Alpha_glm) = colnames(X)
#   
#   return(list(beta_glm = beta_glm,
#               Alpha_glm = Alpha_glm))
# }