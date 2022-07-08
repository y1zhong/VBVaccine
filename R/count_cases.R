## Function count_cases: ------------------------------------------------------
## Convert Type 1 data to Type 2 data, grouping by covariates.
###
###  Input: 
###     1. Type I data: 
###                      ID   VAX_TYPE   AE_NAME    AGE   SEX
###                     201     FLUN      Insomnia   69     F
###                     ...     ...       ...        ...   ...
###                     299     FLU       Chills     68     M
###       
###     2. drug.case: Target vaccine type
###     3. drug.control: Reference vaccine type
###     4. covar_cont: A vector of continuous covariates
###     5. covar_disc: A vector of categorical covariates
###     6. breaks: A list of vectors. Should have the same length as covar_cont.
###  Output:
###     1. Type II data:
###                  VAX_TYPE   AE_NAME   AEyes    AEno   AGE   SEX
###                    FLUN      Insomnia   640     6544    1     F
###                    FLU       Insomnia   200     3291    4     M
###                    ...        ...       ...     ...    ...   ...
###                    FLU       Chills     586     3720    3     M                     
###
###     AEyes: Number of observations that have this AE.
###     AEno: Number of observations that do not have this AE.
###     AE_NAME: The name of the adverse event.
# 79: -------------------------------------------------------------------------
count_cases = function(data, drug.case = drug.case, drug.control = NULL,
                       covar_disc = NULL, covar_cont = NULL, breaks = NULL){
  data = as_tibble(data) 
  ## check the breaks-covar_cont pairs
  if(!length(breaks) == length(covar_cont)){
    stop("The length of breaks does not match that of continuous covariates")
  }
  ## check the basic first three columns
  if(length(names(data)) < 3){
    stop("Unexpected data type")
  }
  ## check the existence of continuous covariates
  if(!is.null(covar_cont)){
    if(!all(covar_cont %in% names(data))){
      stop("nonexistent column") } else {
        if(!all(sapply(data, is.numeric)[covar_cont])){
          stop("Discrete covariates misclassified as continuous")
        }
      }
  }
  ## check the existence of discrete covariates
  if(!is.null(covar_disc)){
    if(!all(covar_disc %in% names(data))){
      stop("nonexistent column") }else {
        if(any(sapply(data, is.numeric)[covar_disc])){
          stop("Continuous covariates misclassified as discrete")
        }
      }
  }
  ## rename the first three columns
  names(data)[1:3] = c('ID', 'VAX_TYPE', 'AE_NAME')
  ## filter by drug.case and drug.control
  if(!is.null(drug.control)){
    drug_list = c(drug.case, drug.control)
    data = data[data$VAX_TYPE %in% drug_list, ]
    }
  ## remove NAs
  data = data[complete.cases(data), ]
  data = data %>%
    mutate(VAX_TYPE = ifelse(VAX_TYPE %in% drug.case,
                                "DrugYes",
                                "DrugNo") )
  ## filter out reports with both case and control vaccines
  ID_No = data %>%
    filter(VAX_TYPE == "DrugNo") %>%
    select(ID) %>%
    unique()
  ID_Yes = data %>%
    filter(VAX_TYPE == "DrugYes") %>%
    select(ID) %>%
    unique()
  Confused_ID = intersect(ID_No[[1]], ID_Yes[[1]])
  
  ## convert all the character variable to factor
  data_temp = data %>%
    filter(!ID %in% Confused_ID) %>%
    mutate_if(sapply(data, is.character), as.factor)
  
  ## for continuous covariates, classifying each obs by argument 'breaks'
  if(!is.null(covar_cont)){
    for(i in 1:length(breaks)){
      breaks_temp = breaks[[i]]
      var_temp = covar_cont[i]
      covar_group = data_temp %>%
        dplyr::select(all_of(var_temp)) %>%
        unlist() %>%
        findInterval(breaks_temp) %>%
        as.factor()
      data_temp = data_temp %>%
        dplyr::select(-all_of(var_temp)) %>%
        tibble({{var_temp}} := covar_group)
      }
    }
    
  data_comp = data_temp
  AE_SET = data_comp %>%
    dplyr::select(AE_NAME) %>%
    unique() %>%
    unlist() %>%
    unname() %>%
    as.character()

  cores = floor(detectCores() / 1.5) 
  cl = makeCluster(cores, setup_strategy = "sequential")
  registerDoParallel(cl)
  results = foreach(i = 1:length(AE_SET),
                    .packages = c("tidyverse"),
                    .combine = bind_rows 
                    ) %dopar% {
    AE = AE_SET[i]
    AE_yes = data_comp %>%
      filter(AE_NAME == AE) %>%
      mutate(AE_NAME = "AEYes") %>%
      distinct(ID, .keep_all = TRUE)
    
    ## a list of target ID corresponding to AE yes 
    ID_list = AE_yes %>%
      dplyr::select(ID) %>%
      unlist() %>%
      unname()
    
    ## filter by ID
    AE_no = data_comp %>%
      filter(! ID %in% ID_list) %>%
      mutate(AE_NAME = "AENo") %>%
      distinct(ID, .keep_all = TRUE)
    
    ## combine together
    data_AE = AE_yes %>%
      bind_rows(AE_no) %>%
      mutate(AE_NAME = as.factor(AE_NAME))
    ## grouping variables
    grp_cols = c("VAX_TYPE", "AE_NAME", covar_cont, covar_disc)
    ## convert into count
    data_count = data_AE %>%
      dplyr::select(-ID) %>%
      group_by_at(grp_cols) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = AE_NAME, values_from = count) %>%
      mutate(AEYes = replace_na(AEYes, 0),
             AENo = replace_na(AENo, 0)) %>%
      mutate(AE_NAME = {{AE}}) %>%
      relocate(VAX_TYPE, AE_NAME, AEYes, AENo)
    data_count
    }
  stopCluster(cl)
  return(results)
}
