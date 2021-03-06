#' @export
VB_Vacc = function(data,  VB=NULL,
                   initial = NULL,
                   mu_alpha=NULL,
                   a_alpha = 0.01,
                   b_alpha  = 0.01,
                   include_intercept=FALSE,
                   ELBO_stop = 1,
                   max_iter=5000,
                   ELBO_diff_tol=1e-6,
                   verbose=0, save_profile = 1){
  
  if (is.null(VB)) {
    AEs = colnames(data$A)
    if (include_intercept){
      data$X = cbind(1, data$X)
    }

    J = ncol(data$A)
    p = ncol(data$X)
    if(is.null(mu_alpha)){
      mu_alpha = matrix(0, ncol = J, nrow=p)
    }

    if(is.null(initial)) {
      glm_coef <- get_glm_coefs(data,cutoff = 10)
      init_beta = glm_coef$beta_glm
      init_alpha = t(glm_coef$Alpha_glm)
    } else {
      #need change beta_glm name
      init_beta = initial$beta_glm
      init_alpha = t(initial$Alpha_glm)
    }

    vb_fit <- Bayes_Vacc_cpp(data$A, data$X, data$V,data$nn,
                         initial_beta = init_beta,
                         initial_alpha = init_alpha,
                         initial_inv_sigma2_beta = rep(1, J),
                         initial_inv_sigma2_alpha = matrix(1, ncol = J, nrow = p),
                         initial_inv_a_beta = 0.5,
                         mu_alpha = mu_alpha,
                         a_alpha = a_alpha,
                         b_alpha  = b_alpha,
                         ELBO_stop = ELBO_stop,
                         ELBO_diff_tol = ELBO_diff_tol,
                         max_iter = max_iter,
                         verbose = verbose,
                         save_profile = save_profile)
    colnames(vb_fit$data$A) = AEs
  } else {
    AEs = colnames(VB$data$A)
    AE_same = intersect(colnames(VB$data$A), colnames(data$A))
    AE_add = colnames(data$A)[!colnames(data$A) %in% colnames(VB$data$A)]
    if(length(AE_same) != 0){
      VB$data$A[, AE_same] = VB$data$A[, AE_same] + data$A[, AE_same]
    }
    if(length(AE_add) != 0){
      VB$data$A = cbind(VB$data$A, data$A[, AE_add])
      AEs = c(AEs, colnames(data$A[, AE_add,drop=FALSE]))
      data$A =  data$A[, AE_add,drop=FALSE]
      glm_coef <- get_glm_coefs(data,cutoff = 10)
      VB$post_mean$beta = c(VB$post_mean$beta, glm_coef$beta_glm)
      VB$post_mean$inv_sigma2_beta = c(VB$post_mean$inv_sigma2_beta, rep(1, length(AE_add)))
      
      VB$post_mean$alpha = cbind(VB$post_mean$alpha, t(glm_coef$Alpha_glm))
      VB$post_mean$inv_sigma2_alpha = cbind(VB$post_mean$inv_sigma2_alpha, matrix(1, ncol = length(AE_add), 
                                                                                  nrow = nrow(VB$post_mean$alpha)))
      VB$hyper$mu_alpha = cbind(VB$hyper$mu_alpha, matrix(0, ncol = length(AE_add), 
                                                          nrow = nrow(VB$post_mean$alpha)))
    }

    VB$data$nn = VB$data$nn + data$nn
    vb_fit <- Bayes_Vacc_cpp(VB$data$A, VB$data$X, VB$data$V,VB$data$nn,
                         initial_beta = VB$post_mean$beta,
                         initial_alpha = VB$post_mean$alpha,
                         initial_inv_sigma2_beta = VB$post_mean$inv_sigma2_beta,
                         initial_inv_sigma2_alpha = VB$post_mean$inv_sigma2_alpha,
                         initial_inv_a_beta = VB$post_mean$inv_a_beta,
                         mu_alpha = VB$hyper$mu_alpha,
                         a_alpha = VB$hyper$a_alpha,
                         b_alpha  = VB$hyper$b_alpha,
                         ELBO_stop = ELBO_stop,
                         ELBO_diff_tol = ELBO_diff_tol,
                         max_iter=max_iter,
                         verbose=verbose, save_profile = save_profile)
    colnames(vb_fit$data$A) = AEs
  }
 return(vb_fit)
}
