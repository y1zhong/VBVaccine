# VB_prep <- function(data, case, control, aggregate = F, age_cutoff = 18, AE_cutoff = 50) {
#   data <- data %>%
#     filter(AGE >= age_cutoff,
#            VAX_LABEL %in% c(case, control))  %>%
#     mutate(SEX = as.numeric(SEX == 'M'),
#            VAX = ifelse(VAX_LABEL == case, 1,0)) %>%
#   arrange(-desc(VAERS_ID))
# 
#   data <- data %>%
#     group_by(AE_NAME) %>%
#     mutate(n = length(unique(VAERS_ID))) %>%
#     filter(n >= AE_cutoff) %>%
#     select(-n) %>%
#     ungroup()
# 
#   data %>%
#     group_by(VAERS_ID) %>%
#     mutate(n = length(unique(VAX_LABEL))) %>%
#     filter(n == 1) %>%
#     ungroup() %>%
#     select(-n) ->
#     data
# 
#   AEs = unique(data$AE_NAME)
#   if (aggregate == T) {
#     data= data %>%
#       mutate(SEX = as.factor(SEX)) %>%
#       select(VAERS_ID, VAX_LABEL, AE_NAME, AGE, SEX)%>%
#       distinct()
#     count_all=count_cases(data, drug.case = case, drug.control = control,
#                   covar_cont = c("AGE"), covar_disc = c("SEX"),
#                   breaks = list(c(30,50,65)))
#                   #breaks = list(c(30,40,50,60,70,80,90)))
# 
#     nn <-  count_all %>%
#       filter(AE_NAME == AEs[1])%>%
#       mutate(nn = AEYes + AENo) %>%
#       pull(nn)
# 
#     X <- count_all %>%
#       distinct(AGE, SEX, DRUG_TYPE) %>%
#       mutate(V = as.numeric(DRUG_TYPE == 'DrugYes'),
#              Intercept = 1) %>%
#       select(Intercept, SEX, AGE, V) %>%
#       mutate(across(where(is.factor), as.character))
#     X <- as.matrix(sapply(X, as.numeric))
# 
#     V <- X[, 4]
#     X <- X[, -4]
# 
# 
#     cores = floor(detectCores() / 1.5)
#     cl = makeCluster(cores, setup_strategy = "sequential")
#     registerDoParallel(cl)
#     A = foreach(i = 1:length(AEs),
#                       .packages = c("tidyverse"),
#                       .combine=cbind
#     ) %dopar% {
#       count_all %>%
#         filter(AE_NAME == AEs[i]) %>%
#         pull(AEYes)
#     }
#     stopCluster(cl)
# 
#     colnames(A) = AEs
#   } else {
#     X <- data %>%
#       distinct()%>%
#       mutate(Intercept = 1) %>%
#       select(Intercept, SEX, AGE, VAX) %>%
#       as.matrix()
# 
# 
# 
#     V <- X[, 4]
#     X <- X[, -4]
# 
#     cores = floor(detectCores() / 1.5)
#     cl = makeCluster(cores, setup_strategy = "sequential")
#     registerDoParallel(cl)
#     A = foreach(i = 1:length(AEs),
#                 .packages = c("tidyverse"),
#                 .combine=cbind
#     ) %dopar% {
#       ifelse(data$AE_NAME == AEs[i], 1, 0)
#     }
#     stopCluster(cl)
# 
#     # A = sapply(AEs, function(AE) {
#     #
#     # })
#     colnames(A) = AEs
#     #A = A[, sort(colnames(A))]
# 
#     nn = rep(1, nrow(A))
#   }
# 
# 
# 
# 
#   # ID_AE_lst = split(data$VAERS_ID, data$AE_NAME)
#   # IDs = sort(unique(data$VAERS_ID))
#   # A = sapply(ID_AE_lst, function(ID_AE) {
#   #   IDs %in% ID_AE + 0
#   # })
# 
# 
#   return(data=list(A=A,
#                    X=X,
#                    V=V,
#                    nn=nn))
# 
# }
# 
# 
# source("~/Dropbox (University of Michigan)/BayesianVaccine/Rcode/yuan/code/count_cases.R")
# load("data/MedDRA.RData")
# load("data/neg_control.RData")
# library(readr)
# library(dplyr)
# library(tidyr)
# library(parallel)
# library(doParallel)
# data <- read_csv("../data/full_dds_covid_flu_2016_20220603.csv",
#                  col_types = cols(DATEDIED = col_skip()))
# data %>%
#   mutate(AE_NAME = replace(AE_NAME,
#                            AE_NAME %in% MedDRA[MedDRA$soc == "Neoplasms benign, malignant and unspecified (incl cysts and polyps)",]$PT,
#                            "Cancer"),
#          MEDDRA_ID = replace(MEDDRA_ID,
#                              grepl("Magnetic resonance imaging", AE_NAME, fixed = TRUE),
#                              unique(MEDDRA_ID[AE_NAME =="Nuclear magnetic resonance imaging"])),
#          AE_NAME = replace(AE_NAME,
#                            grepl("Magnetic resonance imaging", AE_NAME, fixed = TRUE) |
#                              grepl("magnetic resonance imaging", AE_NAME, fixed = TRUE),
#                            "Nuclear magnetic resonance imaging")) ->data
# 
# 
# data <- data %>%
#   mutate(AE_NAME = replace(AE_NAME,
#                            (AE_NAME %in% MedDRA[MedDRA$soc %in% c("Product issues","Investigations","Others"), ]$PT)
#                            & (!AE_NAME %in% neg_control$PT), "Others"))
# 
# data %>%
#   filter(!VAX_NAME == "COVID19 (COVID19 (UNKNOWN))",
#          VAX_TYPE %in% c("FLU3", "FLU4", "COVID19"))%>%
#   filter(AE_NAME %in% MedDRA$PT)%>%
#   filter(SEX != 'U') %>%
#   distinct(VAERS_ID, AGE, SEX, AE_NAME, VAX_LABEL) %>%
#   drop_na() ->
#   data
# 
# 
# 
# age_cutoff = 18
# AE_cutoff = 25
# case="COVID"
# control="FLU"
# data <- data[sample(1:nrow(data), 500000), ]
# data_p <- VB_prep(data, "COVID", "FLU", age_cutoff = 18, AE_cutoff = 50,aggregate = T)
# 
# 
# ##sum nn not equal to total #AEs (one person reported several AEs)
# 
# #data_p <- VB_prep(data, "COVID", "FLU", age_cutoff = 18, AE_cutoff = 50,aggregate = T)
# 
# A=data_p$A
# X=data_p$X
# V=data_p$V
# nn=data_p$nn
# init_beta=glm_coef$beta_glm
# init_alpha=t(glm_coef$Alpha_glm)
# devtools::load_all()
# glm_coef <- get_glm_coefs(data_p$A, data_p$X, data_p$V, data_p$nn)
# mu_alpha =  matrix(0,ncol = ncol(t(glm_coef$Alpha_glm)), nrow=nrow(t(glm_coef$Alpha_glm)), byrow = T)
# VB_fit <- VB_Vacc_cpp(data_p$A, data_p$X, data_p$V,data_p$nn,
#             initial_beta=glm_coef$beta_glm,
#             initial_alpha=t(glm_coef$Alpha_glm),
#             initial_mu_alpha = mu_alpha,
#             include_intercept = F,
#             initial_a_alpha = 0.01, initial_b_alpha  = 0.01,
#             ELBO_stop = 1,
#             max_iter=10000,
#             verbose=1, save_profile = 1)
# plot(VB_fit$trace$ELBO,type = 'l' )
# 
# colnames(A)[order(VB_fit$post_mean$beta,decreasing = T)[1:50]]
# 
# VB_fit_r <- VB_Vacc_r(data_p$A, data_p$X, data_p$V,data_p$nn,
#             init_beta=glm_coef$beta_glm,
#             init_alpha=t(glm_coef$Alpha_glm),
#             include_intercept = F,
#             a_alpha = 0.01, b_alpha  = 0.01,
#             max_iter=2000)
# plot(VB_fit_r$trace,type = 'l' )
