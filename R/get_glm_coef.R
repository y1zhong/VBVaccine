get_glm_coefs = function(data,
                         cutoff = 10) {
  # run glm models
  A=data$A
  X=data$X
  nn=data$nn
  V=data$V
  Alpha_glm = matrix(0, ncol(A), ncol(X) + 1)
  suppressWarnings({
    for (j in 1:ncol(A)) {
      Aj = A[, j]
      model = glm(Aj/nn  ~ 0 + X + V,
                  family = binomial())
      alpha_glm = coef(model)
      alpha_glm = pmin(alpha_glm,cutoff)
      alpha_glm = pmax(alpha_glm,-cutoff)
      
      pos_inf = which(alpha_glm == cutoff)
      alpha_glm[pos_inf] <- rnorm(length(pos_inf), cutoff, 0.5)

      neg_inf = which(alpha_glm == -cutoff)
      alpha_glm[neg_inf] <- rnorm(length(neg_inf), cutoff, 0.5)
      
      Alpha_glm[j,] =  alpha_glm
    }
  })
  
  rownames(Alpha_glm) = NULL
  beta_glm = Alpha_glm[, ncol(Alpha_glm)]
  Alpha_glm = Alpha_glm[,-ncol(Alpha_glm), drop=F]
  colnames(Alpha_glm) = colnames(X)
  
  return(list(beta_glm = beta_glm,
              Alpha_glm = Alpha_glm))
}
