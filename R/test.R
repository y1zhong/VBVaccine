# library(mvtnorm)
# load("../simulation/data/true_values.RData")
# n=1000
# p=2
# J=100
# set.seed(123)
# 
# n_sigs = 20
# beta = c(rep(1,n_sigs), rep(-1,n_sigs), rep(0, J-2*n_sigs))
# 
# sigma2_alpha = matrix(runif((p+1)*J, 0, 1), ncol = J, nrow=(p+1))
# mu_alpha =  matrix(c(rnorm(J, -3, 1), rnorm(p*J,0,1)), ncol = J, nrow=(p+1), byrow = T)
# Alpha = sapply(1:J, function(j) rmvnorm(1, mu_alpha[,j], diag(sigma2_alpha[,j])))
# 
# x = cbind(round(runif(n, 16,80)), rbinom(n,1,0.5)) #generate age within range and 0.5 sex proportion
# x = matrix(rnorm(p*n), nrow = n, ncol = p)
# X = cbind(1,x)    # add intercept
# 
# Alpha = Alpha[1,,drop=F]
# X=X[,1,drop=FALSE]
# 
# p_vacc = 0.7 #probability of getting target vaccine
# V = rbinom(n, 1, p_vacc)
# sum(V)
# 
# 
# # n_new = 50 #observe new data
# # x_new = matrix(rnorm(p*n_new), nrow = n_new, ncol = p)
# # X_new = cbind(1,x_new)
# # V_new = rbinom(n_new, 1, p_vacc)
# 
# Phi <- X%*%Alpha + tcrossprod(V, beta)
# A <- sapply(1:J, function(j) rbinom(n,1,prob = plogis(Phi[,j])))
# colMeans(A)
# 
# 
# glm_coef <- get_glm_coefs(A,X,V,rep(1,nrow(A)))
# #a_alpha=0.1; b_alpha=0.1
# init_beta = glm_coef$beta_glm
# init_alpha=t(glm_coef$Alpha_glm)
# a_alpha=0.01;b_alpha=0.01
# include_intercept=F;max_iter = 100
# mu_alpha =  matrix(c(rep(-5,J), rep(0,p*J)),ncol = J, nrow=(p+1), byrow = T)
# 
# devtools::load_all()
# mu_alpha =  matrix(0,ncol = ncol(t(glm_coef$Alpha_glm)), nrow=nrow(t(glm_coef$Alpha_glm)), byrow = T)
# VB_fit <- VB_Vacc_cpp(A,X,V,rep(1,nrow(A)),
#             initial_beta=glm_coef$beta_glm,
#             initial_alpha=t(glm_coef$Alpha_glm),
#             initial_mu_alpha = mu_alpha,
#             include_intercept = F,
#             initial_a_alpha = 0.01, initial_b_alpha  = 0.01,
#             ELBO_stop = 1,
#             max_iter=1000,
#             verbose=1, save_profile = 1)
# 
# VB_fit_r <- VB_Vacc_r(A,X,V,nn=rep(1,nrow(A)),
#             init_beta=glm_coef$beta_glm,
#             init_alpha=t(glm_coef$Alpha_glm),
#             include_intercept = F,
#             mu_alpha= mu_alpha,
#             a_alpha = 0.01, b_alpha  = 0.01,
#             max_iter=2000)
# plot(VB_fit_r$trace,type = 'l')
# plot(VB_fit$trace$ELBO,type = 'l')
# mean(abs(VB_fit$post_mean$beta - beta)^2)
# mean(abs(glm_coef$beta_glm - beta)^2)
# 
# mean(abs(VB_fit$post_mean$alpha - Alpha)^2)
# mean(abs(t(glm_coef$Alpha_glm) - Alpha)^2)
# 
# 
# mean(abs(VB_fit_r$post_mean$beta - beta)^2)
# mean(abs(glm_coef$beta_glm - beta)^2)
# 
# mean(abs(VB_fit_r$post_mean$alpha - Alpha)^2)
# mean(abs(t(glm_coef$Alpha_glm) - Alpha)^2)
