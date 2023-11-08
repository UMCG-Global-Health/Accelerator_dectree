################################################################################ 
# This script runs the decision tree model for:                                #
# "Cost-effectiveness of test-and-treat strategies to reduce the antibiotic    #
#  prescription rate for acute febrile illness in primary healthcare clinics   #
#  in Africa"                                                                  #
#                                                                              #
#                                                                              # 
# Authors:                                                                     #
#     - Pim W.M. van Dorst                                                     # 
#     - Simon van der Pol                                                      #
#     - Piero Olliaro                                                          #
#     - Sabine Dittrich                                                        #
#     - Juvenal Nkeramahame                                                    # 
#     - Maarten J. Postma                                                      #
#     - Cornelis Boersma                                                       #
#     - Antoinette D.I. van Asselt                                             #
#                                                                              #
################################################################################ 

## Clean the workspace
rm(list=ls())

## Load packages
library(dplyr)
library(tidyverse)
library(here)
library(ggplot2)
library(gridExtra)
library(gt)
library(xlsx)

## Appointing to the working directory
#  - Make sure the script + data_files are in the same directory
setwd(here())


### FUNCTION DECISION TREE
dec_tree <- function(p_table, df_c_personnel_i, c_equipment, c_tests_exc, country_i, df_conv, n_cohort, prob = FALSE, runs = 1000){

  ## p_table            : aggregated input dataframe, specific to the country and scenario. Based on patient level data. 
  ## df_c_personnel_i   : dataframe of cost of personnel per country
  ## c_equipment        : overview of cost of equipment
  ## c_tests_exc        : dataframe of cost per test per country
  ## country_i          : country to run the model for
  ## df_conv            : dataframe with currency conversion
  ## n_cohort           : number of patient in cohort
  ## prob               : to run deterministic or probabilistic analysis. Default is deterministic (FALSE)
  ## runs               : number of runs for the probabilistic analysis. Default is 1000 runs
  
  
  ## Check if p_table is empty
  if(length(p_table) == 1){
    
    df_total_sum <- 0
  
  } else {
    
  n_train_n <- 5     # number of nurses trained
  n_train_d <- 2     # number of doctors trained
  n_train_l <- 2     # number of lab-technicians trained
  n_sites   <- 1     # number of sites modelled
  
  ## Calculation of cost of hospitalization
  ppp <- data.frame(
    country = c("uganda", "burkina-faso", "ghana"),
    c_hosp = c(124.83, 51.14, 38.03)) %>%
    filter(country == country_i) %>%
    dplyr::select(c_hosp) %>%
    unlist()
  
  ## Calculate cost of training
  df_c_training <- df_c_personnel_i %>%
    dplyr::mutate(c_training = h_rate * time_t,
                  c_training = base::ifelse(role == "doctor", c_training * n_train_d,
                                            base::ifelse(role == "nurse", c_training * n_train_n,
                                                         c_training * n_train_l))) %>%
    group_by(rarm) %>%
    dplyr::summarise(c_training = sum(c_training)/2)
  
  ## For deterministic analysis
  if (prob == FALSE) {
    
    ## Calculate personnel cost per patient 
    df_c_personnel <- df_c_personnel_i %>% 
      dplyr::mutate(cost_pp = h_rate * time_p) %>%
      group_by(rarm) %>%
      dplyr::summarise(cost_pp = sum(cost_pp))
    
  ## Collect the list of parameters from the input file
  df_consequences <- p_table %>%
    
    # remove the probabilistic elements
    dplyr::select(-contains("alpha"), -contains("beta")) %>%
    
    # Create a longer table to join the table with the Test cost table
    pivot_longer(cols = c("malaria_rapid_test":"mean_hb_test"), names_to = "test", values_to = "mean") %>%
    mutate(test = str_remove(test, "mean_")) %>%
    left_join(., c_tests_exc[c_tests_exc$country_t == country_i, ], by = c("test" = 'Tests in scope of trial')) %>%
    
    # Calculate the mean cost of tests based on the average number of tests done per decision tree branch
    dplyr::mutate(c_dx_mean = mean * cost) %>%
    dplyr::select(-c("country_t", "mean", "cost")) %>%
    pivot_wider(names_from = "test", values_from = "c_dx_mean") %>%
    
    # Per branch, calculate the total cost of tests
    rowwise() %>%  
    dplyr::mutate(c_dx_mean = sum(c_across("malaria_rapid_test":"hb_test"), na.rm = TRUE)) %>%
    dplyr::select(-c("malaria_rapid_test":"hb_test")) %>%
    
    # Determine total probability per branch (end-to-end)
    dplyr::mutate(p_tot = p_resp_mean * p_mal_mean * p_abx_mean * p_adh_mean * p_f_mean * p_sae_unfav_mean) %>%
    
    # Add the cost of personnel
    left_join(., df_c_personnel, by = c("group" = "rarm")) %>%
    
    # Calculate all the outcome parameters
    dplyr::mutate(c_dx_total = n_cohort * p_tot * c_dx_mean,   # Cost of diagnostic tests
                  c_tx_total = n_cohort * p_tot * c_tx_mean,   # cost of therapuetics
                  c_labour_total = n_cohort * p_tot * cost_pp, # Cost of labour
                  c_oop_add_abx = n_cohort * p_tot * mean_add_abx * total_cost_mean_add,  # Cost of out-of-pocket expenses 
                  n_oop_add_abx = n_cohort * p_tot * mean_add_abx, # Number of antibiotics prescribed, assuming 1 course of abx is given per time 
                  ddd_total = n_cohort * p_tot * cons_ddd + n_cohort * p_tot * mean_add_abx * ddd_prescribed_add, # Number of DDD
                  dot_total = n_cohort * p_tot * cons_dot + n_cohort * p_tot * mean_add_abx * DOT_add, # Number of days on therapy
                  n_fav = base::ifelse(fav == TRUE, n_cohort * p_tot, 0),  # Number of favourable outcomes
                  n_nonfav = base::ifelse(fav == FALSE, n_cohort * p_tot, 0), # Number of unfavourable outcomes
                  daly = base::ifelse(fav == FALSE & sae == FALSE, n_nonfav * 0.051/52, # Number of DALYs
                                     base::ifelse(fav == FALSE & sae == TRUE, n_nonfav * 0.133/52, 0)),
                  c_nonfav = base::ifelse(fav == TRUE, 0, # Cost in case of unfavorable outcome
                                          base::ifelse(fav == FALSE & sae == FALSE, n_nonfav * (c_tx_mean + cost_pp), 
                                                       base::ifelse(fav == FALSE & sae == TRUE, n_nonfav * (ppp[[1]] + c_tx_mean + cost_pp),0))),
                  n_abx = n_cohort * p_tot * n_abx + n_oop_add_abx, # Total number of antibiotics prescribed
                  abx_presc_rate = n_cohort * p_tot * abx_prescribed, # Antibiotic prescription rate
                  c_total = n_cohort * p_tot * (c_dx_mean + c_tx_mean + cost_pp) + c_oop_add_abx + c_nonfav) %>% # Total cost
    
    # Enable calculation per group (intervention or control)
    group_by(group) %>%
    
    # Calculate the sum of all outcome parameters
    dplyr::summarise(across(.cols = c("c_dx_total":"c_total", "n_abx"), 
                            .fns = ~sum(.x, na.rm = TRUE),
                            .names = "{.col}")) %>%
    
    # Add one-off costs to the table:
    #  - Cost of training
    #  - Cost of equipment
    
    left_join(., df_c_training, by = c("group" = "rarm")) %>%
    left_join(., c_equipment, by = c("group" = "group")) %>%
    dplyr::mutate(c_equipment = c * n_sites,
                  c_total = base::ifelse(group == "int",
                                           c_training + c_equipment + c_total,
                                           c_total)) %>%
    dplyr::select(-c)
  
  ## Define the output data-frame
  df_total_sum <- df_consequences
  
  } else {
    
    # For probabilistic calculation

    ## Calculate cost per patient 
    df_c_personnel <- list(cost_pp_doc_con = rnorm(runs, mean = df_c_personnel_i$time_p[[1]], sd = 0.2 * df_c_personnel_i$time_p[[1]]) * df_c_personnel_i$h_rate[[1]],
                           cost_pp_doc_int = rnorm(runs, mean = df_c_personnel_i$time_p[[2]], sd = 0.2 * df_c_personnel_i$time_p[[2]]) * df_c_personnel_i$h_rate[[2]], 
                           cost_pp_nur_con = rnorm(runs, mean = df_c_personnel_i$time_p[[3]], sd = 0.2 * df_c_personnel_i$time_p[[3]]) * df_c_personnel_i$h_rate[[3]], 
                           cost_pp_nur_int = rnorm(runs, mean = df_c_personnel_i$time_p[[4]], sd = 0.2 * df_c_personnel_i$time_p[[4]]) * df_c_personnel_i$h_rate[[4]], 
                           cost_pp_lt_con = rnorm(runs, mean = df_c_personnel_i$time_p[[5]], sd = 0.2 * df_c_personnel_i$time_p[[5]]) * df_c_personnel_i$h_rate[[5]], 
                           cost_pp_lt_int = rnorm(runs, mean = df_c_personnel_i$time_p[[6]], sd = 0.2 * df_c_personnel_i$time_p[[6]]) * df_c_personnel_i$h_rate[[6]])
                         
    ## Calculate the sum cost of all personnel per patient                          
    df_c_personnel_int <- df_c_personnel[[2]] + df_c_personnel[[4]] + df_c_personnel[[6]] 
    df_c_personnel_con <- df_c_personnel[[1]] + df_c_personnel[[3]] + df_c_personnel[[5]] 
    
    ## Calculate DALY losses in case of adverse events
    daly_prob <- list(
      no_sae = rgamma(runs, shape = (0.051/52)^2 / (0.2 * (0.051/52))^2, rate = (0.051/52) / (0.2 * (0.051/52))^2),
      sae = rgamma(runs, shape = (0.133/52)^2 / (0.2 * (0.133/52))^2, rate = (0.133/52) / (0.2 * (0.133/52))^2)
    ) 
    
    ## Calculate cost of hospitalisation
    hosp_prob <- list(
      hosp_p = rgamma(runs, shape =  ppp[[1]]^2 / (0.2 * ppp[[1]])^2, rate = ppp[[1]] / (0.2 * ppp[[1]])^2)
    )
    
    ## Calculate cost of additional antibiotics taken
    add_abx <- p_table %>%
      dplyr::select(c("total_cost_mean_add":"total_cost_rate_add")) %>%
      dplyr::filter(total_cost_rate_add > 0) %>%
      unique() %>%
      dplyr::mutate(total_cost_prob_add = list((rgamma(runs, total_cost_shape_add, total_cost_rate_add)))) %>%
      dplyr::select(total_cost_prob_add) %>%
      unlist() %>%
      as.numeric() %>%
      as.list()
    
    ## Calculate the cost of tests
    df_dx <- p_table %>%
      
      # remove the deterministic elements
      dplyr::select(-contains("mean")) %>%
      
      dplyr::select(group, resp, mal, c("alpha_typhoid_test":"beta_hb_test")) %>%
    
      # Create only unique data.frame
      unique() %>%

      group_by(group, resp, mal) %>%
      
      # Determine the probability of using a test per decision tree branch
      dplyr::mutate(prob_crp_test = list(rbeta(runs, alpha_crp_test, beta_crp_test)), 
                    prob_malaria_rapid_test = list(rbeta(runs, malaria_rapid_test_alpha, malaria_rapid_test_beta)), 
                    prob_typhoid_test = list(rbeta(runs, alpha_typhoid_test, beta_typhoid_test)),
                    prob_a_strept_test = list(rbeta(runs, alpha_a_strept_test, beta_a_strept_test)),
                    prob_influenza_test = list(rbeta(runs, alpha_influenza_test, beta_influenza_test)),
                    prob_rsv_test = list(rbeta(runs, alpha_rsv_test, beta_rsv_test)),
                    prob_strept_pneu_test = list(rbeta(runs, alpha_strept_pneu_test, beta_strept_pneu_test)),
                    prob_scrub_test = list(rbeta(runs, alpha_scrub_test, beta_scrub_test)),
                    prob_hiv_test = list(rbeta(runs, alpha_hiv_test, beta_hiv_test)),
                    prob_blood_cul_test = list(rbeta(runs, alpha_blood_cul_test, beta_blood_cul_test)),
                    prob_blood_film_test = list(rbeta(runs, alpha_blood_film_test, beta_blood_film_test)),
                    prob_urine_mcs_test = list(rbeta(runs, alpha_urine_mcs_test, beta_urine_mcs_test)),
                    prob_xray_test = list(rbeta(runs, alpha_xray_test, beta_xray_test)),
                    prob_sput_afb_test = list(rbeta(runs, alpha_sput_afb_test, beta_sput_afb_test)),
                    prob_stool_mic_test = list(rbeta(runs, alpha_stool_mic_test, beta_stool_mic_test)),
                    prob_sput_gram_test = list(rbeta(runs, alpha_sput_gram_test, beta_sput_gram_test)),
                    prob_widal_test = list(rbeta(runs, alpha_widal_test, beta_widal_test)),
                    prob_dengue_test = list(rbeta(runs, alpha_dengue_test, beta_dengue_test)),
                    prob_urine_test = list(rbeta(runs, alpha_urine_test, beta_urine_test)),
                    prob_cbc_test = list(rbeta(runs, alpha_cbc_test, beta_cbc_test)),
                    prob_hb_test = list(rbeta(runs, alpha_hb_test, beta_hb_test))) %>%
      
      # Drop fields not used anymore
      dplyr::select(-c("alpha_typhoid_test":"beta_hb_test")) %>%
      
      # Create a longer table to join the table with the test cost table
      pivot_longer(cols = c("prob_crp_test":"prob_hb_test"), names_to = "test", values_to = "prob") %>%
      dplyr::mutate(test = str_remove(test, "prob_")) %>%
      left_join(., c_tests_exc[c_tests_exc$country_t == country_i,], by = c("test" = 'Tests in scope of trial')) %>%
      dplyr::mutate(df_cost_test = list(0))
    
    for(i in 1:168){
      
      ## Calculate the test cost for each of the PSA runs
      df_dx$df_cost_test[[i]] <- unlist(lapply(df_dx$prob[[i]], function(.x) df_dx$cost[[i]] * .x))
      
    }
    
    ## Create a wide table
    df_dx_cont <- df_dx %>%
      dplyr::select(-c("country_t", "prob", "cost")) %>%
      pivot_wider(names_from = "test", values_from = "df_cost_test") %>%
      ungroup() %>%
      dplyr::mutate(c_dx_prob = list(0))
    
    for(i in 1:8){
      
      ## Calculate the total test cost per branch
      df_dx_cont$c_dx_prob[[i]] <-  unlist(df_dx_cont$crp_test[[i]] + df_dx_cont$malaria_rapid_test[[i]] + df_dx_cont$typhoid_test[[i]] + df_dx_cont$a_strept_test[[i]] 
      + df_dx_cont$influenza_test[[i]] + df_dx_cont$rsv_test[[i]] + df_dx_cont$strept_pneu_test[[i]] + df_dx_cont$scrub_test[[i]] + df_dx_cont$hiv_test[[i]]  
      +  df_dx_cont$blood_cul_test[[i]] + df_dx_cont$blood_film_test[[i]] + df_dx_cont$urine_mcs_test[[i]] + df_dx_cont$xray_test[[i]] + df_dx_cont$sput_afb_test[[i]]  
      +   df_dx_cont$stool_mic_test[[i]] + df_dx_cont$sput_gram_test[[i]] + df_dx_cont$widal_test[[i]] + df_dx_cont$dengue_test[[i]] + df_dx_cont$urine_test[[i]]  
      +   df_dx_cont$cbc_test[[i]] + df_dx_cont$hb_test[[i]])
      
    }
    
    ## Calculate final list of dx cost
    df_dx_f <- df_dx_cont %>% 
      dplyr::select(-c("crp_test":"hb_test")) 
      
    df_consequences <- p_table %>%
      
      # remove the deterministic elements and columns 
      dplyr::select(-contains("mean"), -c("malaria_rapid_test":"beta_hb_test"), -c("total_cost_mean_add":"total_cost_rate_add")) %>%
      
      # Join the calculation of Test cost per patient
      left_join(., df_dx_f, by = c("group" = "group", "resp" = "resp", "mal" = "mal")) %>%
    
      # Calculate the probability for respiratory (and non-respiratory)
      dplyr::mutate(p_resp_prob = list(rbeta(runs, alpha_resp, beta_resp)),
                    p_resp_prob = base::ifelse(resp == TRUE, p_resp_prob, lapply(p_resp_prob, function(.x) 1 - .x))) %>%
      
      dplyr::select(-c("alpha_resp","beta_resp")) %>%
      
      # calculate the probability for malaria (and non-malaria)
      group_by(group, resp) %>%
      dplyr::mutate(p_mal_prob = list(rbeta(runs, alpha_mal, beta_mal)),
                    p_mal_prob = base::ifelse(mal == TRUE, p_mal_prob, lapply(p_mal_prob, function(.x) 1-.x))) %>%
      ungroup() %>%
      
      dplyr::select(-c("alpha_mal","beta_mal")) %>%
      
      # Calculate the probability for prescribing antibiotics (and no antibiotics) 
      group_by(group, resp, mal) %>%
      dplyr::mutate(p_abx_prob = list(rbeta(runs, alpha_abx, beta_abx)),
                    p_abx_prob = base::ifelse(abx == FALSE, lapply(p_abx_prob, function(.x) 1 - .x), p_abx_prob)) %>%
                      
      ungroup() %>%
      dplyr::select(-c("alpha_abx","beta_abx")) %>%
      
      # Calculate the  probabilities for adherence + cost of therapeutics
      group_by(group, resp, mal, abx) %>%
      dplyr::mutate(p_adh_prob = list(rbeta(runs, alpha_adh, beta_adh)),
                    c_tx_prob = list(rgamma(runs, shape_tx, rate_tx)),
                    p_adh_prob = base::ifelse(adh == TRUE, p_adh_prob, lapply(p_adh_prob, function(.x) 1 - .x))) %>%
      ungroup() %>%
      dplyr::select(-c("alpha_adh","beta_adh", "shape_tx", "rate_tx")) %>%
      
      # Calculate the probabilities for favorable outcome (and unfavourable)
      group_by(group, resp, mal, abx, adh) %>%
      dplyr::mutate(p_f_prob = list(rbeta(runs, alpha_f, beta_f)),
                    p_f_prob = base::ifelse(fav == TRUE, p_f_prob, lapply(p_f_prob, function(.x) 1 - .x)),
                    p_f_prob = base::ifelse(is.na(p_f_prob) == TRUE, list(rep(0,runs)), p_f_prob)) %>%
      
      ungroup() %>%
      dplyr::select(-c("alpha_f","beta_f")) %>%
      
      # Calculate the probability of a sae 
      group_by(group, resp, mal, abx, adh, fav) %>%
      dplyr::mutate(p_sae_prob = list(rbeta(runs, alpha_sae_unfav, beta_sae_unfav))) %>%
      
      # In case sae is false, take 1 - prob
      dplyr::mutate(p_sae_prob = base::ifelse(fav == FALSE & sae == TRUE, p_sae_prob, 
                                              base::ifelse(fav == FALSE & sae == FALSE, lapply(p_sae_prob, function(.x) 1 - .x), list(0))),
                    
                    # Correct for favorable outcome which is always without SAE 
                    p_sae_prob = base::ifelse(fav == TRUE & sae == TRUE, list(rep(0,runs)), 
                                              base::ifelse(fav == TRUE & sae == FALSE, list(rep(1,runs)), p_sae_prob)),
                    p_tot_prob = list(0)) %>%
      ungroup() %>%
      dplyr::select(-c("alpha_sae_unfav","beta_sae_unfav")) 
    
    ## Calculate the total probability for each branch for the PSA
    for (i in 1:nrow(df_consequences)){
      
      df_consequences$p_tot_prob[[i]] <- df_consequences$p_resp_prob[[i]] * df_consequences$p_mal_prob[[i]] * df_consequences$p_abx_prob[[i]] * df_consequences$p_adh_prob[[i]] * df_consequences$p_f_prob[[i]] * df_consequences$p_sae_prob[[i]] 
  
    }
    
    df_consequences_f <- df_consequences %>%
        
      # For non-adherence, calculate the total cost and additional DOT and DDD
        group_by(group, resp, mal, abx, adh) %>%
        dplyr::mutate(p_add_abx = list(rbeta(runs, alpha_add_abx, beta_add_abx)),
                      p_add_abx = base::ifelse(adh == TRUE, list(rep(0,runs)), p_add_abx)) %>%
        
        ungroup() %>%
        dplyr::select(-c("alpha_add_abx","beta_add_abx")) %>%
        
      # Add the cost of personnel
      dplyr::mutate(df_c_personnel = base::ifelse(group == "int", list(df_c_personnel_int), list(df_c_personnel_con))) 
      
      ## Calculate the separate outcome parameters 
      df_total <-  df_consequences_f %>%
        dplyr::select(c("group":"sae")) %>%
        dplyr::mutate(c_dx_total = list(rep(0,runs)),
                        c_tx_total = list(rep(0,runs)),
                        c_labour_total = list(rep(0,runs)),
                        c_oop_add_abx = list(rep(0,runs)),
                        n_oop_add_abx = list(rep(0,runs)),
                        ddd_total = list(rep(0,runs)),
                        dot_total = list(rep(0,runs)),
                        n_fav = list(rep(0,runs)),
                        n_nonfav = list(rep(0,runs)),
                        daly = list(rep(0,runs)),
                        c_nonfav = list(rep(0,runs)),
                        n_abx = list(rep(0,runs)),
                        abx_presc_rate = list(rep(0,runs)),
                        c_total = list(rep(0,runs)))
        
      ## Transform list into data frame list for additional antibiotics taken
      l_add_abx <- rep(list(unlist(add_abx)), nrow(df_total))
      
      ## Calculate for each run the different outcome parameters
      df_total$c_dx_total <- lapply(Map('*', df_consequences_f$p_tot_prob, df_consequences_f$c_dx_prob), "*", n_cohort)
      df_total$c_tx_total <- lapply(Map('*', df_consequences_f$p_tot_prob, df_consequences_f$c_tx_prob), "*", n_cohort)
      df_total$c_labour_total <- lapply(Map('*', df_consequences_f$p_tot_prob, df_consequences_f$df_c_personnel), "*", n_cohort)
      df_total$c_oop_add_abx <- lapply(Map('*', l_add_abx , Map('*', df_consequences_f$p_tot_prob, df_consequences_f$p_add_abx)), "*", n_cohort)
      df_total$n_oop_add_abx <- lapply(Map('*', df_consequences_f$p_tot_prob, df_consequences_f$p_add_abx), "*", n_cohort) # assuming 1 abx is given
      df_total$ddd_total <- Map('+', lapply(Map('*', df_consequences_f$p_tot_prob, df_consequences_f$cons_ddd), "*", n_cohort), lapply(Map('*', Map('*', df_consequences_f$p_tot_prob, df_consequences_f$p_add_abx), df_consequences_f$ddd_prescribed_add), "*", n_cohort))
      df_total$dot_total <- Map('+', lapply(Map('*', df_consequences_f$p_tot_prob, df_consequences_f$cons_dot), "*", n_cohort), lapply(Map('*', Map('*', df_consequences_f$p_tot_prob, df_consequences_f$p_add_abx), df_consequences_f$DOT_add), "*", n_cohort))
      df_total$n_fav <- base::ifelse(df_consequences_f$fav == TRUE, lapply(df_consequences_f$p_tot_prob, "*", n_cohort), list(rep(0,runs)))
      df_total$n_nonfav <- base::ifelse(df_consequences_f$fav == FALSE, lapply(df_consequences_f$p_tot_prob, "*", n_cohort), list(rep(0,runs)))
      df_total$daly <- base::ifelse(df_consequences_f$fav == FALSE & df_consequences_f$sae == FALSE, lapply(df_total$n_nonfav, "*", daly_prob[[1]]),
                                              base::ifelse(df_consequences_f$fav == FALSE & df_consequences_f$sae == TRUE, lapply(df_total$n_nonfav, "*", daly_prob[[2]]), list(rep(0,runs))))
      df_total$c_nonfav <-  base::ifelse(df_consequences_f$fav == TRUE, list(rep(0,runs)),
                                                   base::ifelse(df_consequences_f$fav == FALSE & df_consequences_f$sae == FALSE, Map("*", df_total$n_nonfav , Map('+', df_consequences_f$c_tx_prob, df_consequences_f$df_c_personnel)), 
                                                                base::ifelse(df_consequences_f$fav == FALSE & df_consequences_f$sae == TRUE, Map("*", df_total$n_nonfav , lapply(Map('+', df_consequences_f$c_tx_prob, df_consequences_f$df_c_personnel), "+", hosp_prob[[1]])),
                                                                             list(rep(0,runs)))))
      df_total$n_abx <- Map("+", lapply(Map('*', df_consequences_f$p_tot_prob, df_consequences_f$n_abx), "*", n_cohort), df_total$n_oop_add_abx)
      df_total$abx_presc_rate <- lapply(Map('*', df_consequences_f$p_tot_prob, df_consequences_f$abx_prescribed), "*", n_cohort)
      df_total$c_total <- Map("+", df_total$c_nonfav, 
                              Map("+", df_total$c_oop_add_abx, 
                                  Map("+", df_total$c_labour_total, 
                                      Map("+", df_total$c_dx_total, df_total$c_tx_total))))
      
     df_list <- c(
       c_dx_total = list(rep(0,runs)),
        c_tx_total = list(rep(0,runs)),
        c_labour_total = list(rep(0,runs)),
        c_oop_add_abx = list(rep(0,runs)),
        n_oop_add_abx = list(rep(0,runs)),
        ddd_total = list(rep(0,runs)),
        dot_total = list(rep(0,runs)),
        n_fav = list(rep(0,runs)),
        n_nonfav = list(rep(0,runs)),
        daly = list(rep(0,runs)),
        c_nonfav = list(rep(0,runs)),
        n_abx = list(rep(0,runs)),
        abx_presc_rate = list(rep(0, runs)),
        c_total = list(rep(0,runs)))
      
      df_total_sum <- c(
        int = list(df_list),
        con = list(df_list)
      )
      
      df_total_int <- df_total %>%
        dplyr::filter(group == "int")
      
      df_total_con <- df_total %>%
        dplyr::filter(group == "con")
    
     ## Add all consequences up to come to a total value per run per outcome parameters
      for (i in (1:nrow(df_total_int))){

        df_total_sum$int$c_dx_total <- df_total_sum$int$c_dx_total + df_total_int$c_dx_total[[i]]
        df_total_sum$int$c_tx_total <- df_total_sum$int$c_tx_total + df_total_int$c_tx_total[[i]]
        df_total_sum$int$c_labour_total <- df_total_sum$int$c_labour_total + df_total_int$c_labour_total[[i]]
        df_total_sum$int$c_oop_add_abx <- df_total_sum$int$c_oop_add_abx + df_total_int$c_oop_add_abx[[i]]
        df_total_sum$int$n_oop_add_abx <- df_total_sum$int$n_oop_add_abx + df_total_int$n_oop_add_abx[[i]]
        df_total_sum$int$ddd_total <- df_total_sum$int$ddd_total + df_total_int$ddd_total[[i]]
        df_total_sum$int$dot_total <- df_total_sum$int$dot_total + df_total_int$dot_total[[i]]
        df_total_sum$int$n_fav <- df_total_sum$int$n_fav + df_total_int$n_fav[[i]]
        df_total_sum$int$n_nonfav <- df_total_sum$int$n_nonfav + df_total_int$n_nonfav[[i]]
        df_total_sum$int$daly <- df_total_sum$int$daly + df_total_int$daly[[i]]
        df_total_sum$int$c_nonfav <- df_total_sum$int$c_nonfav + df_total_int$c_nonfav[[i]]
        df_total_sum$int$n_abx <- df_total_sum$int$n_abx + df_total_int$n_abx[[i]]
        df_total_sum$int$c_total <- df_total_sum$int$c_total + df_total_int$c_total[[i]]
        df_total_sum$int$abx_presc_rate <- df_total_sum$int$abx_presc_rate + df_total_int$abx_presc_rate[[i]]
        
      }
      
      for (i in (1:nrow(df_total_con))){
        
        df_total_sum$con$c_dx_total <- df_total_sum$con$c_dx_total + df_total_con$c_dx_total[[i]]
        df_total_sum$con$c_tx_total <- df_total_sum$con$c_tx_total + df_total_con$c_tx_total[[i]]
        df_total_sum$con$c_labour_total <- df_total_sum$con$c_labour_total + df_total_con$c_labour_total[[i]]
        df_total_sum$con$c_oop_add_abx <- df_total_sum$con$c_oop_add_abx + df_total_con$c_oop_add_abx[[i]]
        df_total_sum$con$n_oop_add_abx <- df_total_sum$con$n_oop_add_abx + df_total_con$n_oop_add_abx[[i]]
        df_total_sum$con$ddd_total <- df_total_sum$con$ddd_total + df_total_con$ddd_total[[i]]
        df_total_sum$con$dot_total <- df_total_sum$con$dot_total + df_total_con$dot_total[[i]]
        df_total_sum$con$n_fav <- df_total_sum$con$n_fav + df_total_con$n_fav[[i]]
        df_total_sum$con$n_nonfav <- df_total_sum$con$n_nonfav + df_total_con$n_nonfav[[i]]
        df_total_sum$con$daly <- df_total_sum$con$daly + df_total_con$daly[[i]]
        df_total_sum$con$c_nonfav <- df_total_sum$con$c_nonfav + df_total_con$c_nonfav[[i]]
        df_total_sum$con$n_abx <- df_total_sum$con$n_abx + df_total_con$n_abx[[i]]
        df_total_sum$con$c_total <- df_total_sum$con$c_total + df_total_con$c_total[[i]]
        df_total_sum$con$abx_presc_rate <- df_total_sum$con$abx_presc_rate + df_total_con$abx_presc_rate[[i]]
        
      }
      
      for (i in (1:runs)){
        
      ## Add one-off equipment and training costs
      df_total_sum$int$c_total[[i]] <- df_total_sum$int$c_total[[i]] + c_equipment$c[[2]] + df_c_training$c_training[[2]]
      
      }
      
  }
  }
  
  ## dataframe with results
  return(df_total_sum)
  
}



## Function to run the decision tree function
run_dec_tree <- function(country_i, data_list, scenario, prob_i, c_tests_exc_i){
  
  ## country_i      : country to run the model for
  ## data_list      : list of 6 dataframes (for 6 scenarios) per country
  ## scenario       : ability to run tree for different scenarios
  ##   - scenario 1 : base case
  ##   - scenario 2 : mal
  ##   - scenario 3 : mal + crp
  ##   - scenario 4 : mal + cbc
  ##   - scenario 5 : mal + crp + cbc 
  ##   - scenario 6 : mal + all as per guidelines
  ## prob_i         : indicate whether it is probabilistic (TRUE) or deterministic (FALSE)
  ## c_tests_exc_i  : dataframe with cost of tests
  
  
  ## Conversion rates local currency to dollars using PPP 2021 -------
  df_conv_global <- data.frame(
    country = c("uganda","ghana","burkina-faso"),
    conv_f = c(0.00076, 0.43478, 0.00484))
  
  ## Filter the conversion table based on the country
  df_conv <- df_conv_global %>%
    dplyr::filter(country == country_i) %>%
    dplyr::select(conv_f) %>%
    unlist()
  
  ## Cost of consultation and related research for intervention group and control group + cost of training
  df_c_pers <- data.frame(
    country = c("uganda","ghana","burkina-faso"),
    
    # Monthly wages for personnel per country in local currency
    doctor =    c(2300000, 5527.55, 611000),
    nurse =     c(1030000, 1296.97, 295000),
    lab_staff = c(1030000, 1296.97, 300000)) %>%
    
    # Create a longer data-frame per country
    pivot_longer(c("doctor","nurse","lab_staff"), names_to = "role", values_to = "m_rate") %>%
    filter(country == country_i) %>%
    
    # Transform the time spend per patient from local currency to International $ and calculate hourly rate
    dplyr::mutate(h_rate = m_rate * df_conv * 12 / 52 / 40)  
  
  ## Cost of personnel
  if (country_i == "burkina-faso") {
    
    if (scenario == 1 | scenario == 6){
      t_base <- 45
      t_base_l <- 20
      
    }else if (scenario == 2){
      t_base <- 20
      t_base_l <- 0
      
    }else{
      t_base <- 32.5
      t_base_l <- 0
      
    }
    
    df_time <- data.frame(
      #           Nur_c, Nur_i, Doc_c, Doc_i, Lab_c, lab_i                 
      time_p =   c(20/60 , t_base/60,  0, 0, 0, t_base_l/60), # Time per patient
      time_t =   c(0     , 8,     0, 0, 0, 12))
    
  } else {
    
    if (scenario == 1 | scenario == 6){
      t_base <- 45
      t_base_l <- 20
      
    }else if (scenario == 2){
      t_base <- 20
      t_base_l <- 0
      
    }else{
      t_base <- 32.5
      t_base_l <- 0
      
    }
    df_time <- data.frame(
      #           Nur_c, Nur_i, Doc_c, Doc_i, Lab_c, lab_i                 
      time_p =   c(0 ,   0,  20/60, t_base/60, 0, t_base_l/60),   # Time per patient
      time_t =   c(0 ,   0,      0,     8,     0,  12))
    
  } 
  
  ## Create a grid for time per patient
  df_h_p_patient <- expand.grid(
    rarm = c("con", "int"),
    role = c("nurse", "doctor", "lab_staff")) %>%
    cbind(., df_time)
  
  df_c_personnel <- left_join(df_c_pers, df_h_p_patient, by = c("role" = "role"))
  
  ## Determine the cost of CRP and CBC analyzers, depending on the senario whether it is applicable or not
  c_crp <- base::ifelse(scenario == 2 | scenario == 4, 0, 1122.75)
  c_cbc <- base::ifelse(scenario == 2 | scenario == 3, 0, 2544.91)
  
  ## CBC is standard of care in Ghana, no new equipment needed
  c_cbc <- base::ifelse(country_i == "ghana", 0, c_cbc)
  
  ## Get the overview of the costs of test equipment
  c_equipment <- expand_grid(
    group = c("int","con"),
    c_cat = c("c_equipment_crp","c_equipment_cbc")) %>%
    dplyr::mutate(c = c(c_crp/10, c_cbc/10, 0, 0)) %>%
    group_by(group) %>%
    dplyr::summarise(across(.cols = "c",
                            .fns = ~round(sum(.),2),
                            .names = "{.col}"),
                     .groups = "drop")
  
  data_input_i <- data_list[[scenario]]
  
  ## Input all data in decision tree
  dec_tree_outcome <- dec_tree(p_table = data_input_i, df_c_personnel_i = df_c_personnel, c_equipment = c_equipment, c_tests_exc = c_tests_exc_i,
                               country_i = country_i, df_conv = df_conv, n_cohort = 1000, prob = prob_i, runs = 1000)
  
  return(dec_tree_outcome)
  
}

## Load the country and scenario specific aggregated data in the folder
df_burkina <- df_ghana <- df_uganda <- list(scenario_1 = NA, 
                   scenario_2 = NA,
                   scenario_3 = NA,
                   scenario_4 = NA,
                   scenario_5 = NA,
                   scenario_6 = NA)
for (i in 1:6){
 
  df_burkina[[i]] <- read.xlsx("Data_Burkina_Faso.xlsx", sheetIndex = i, sheetName= str_c("scenario_", i))
  df_ghana[[i]] <- read.xlsx("Data_Ghana.xlsx", sheetIndex = i, sheetName= str_c("scenario_", i))
  df_uganda[[i]] <- read.xlsx("Data_Uganda.xlsx", sheetIndex = i, sheetName= str_c("scenario_", i))
  
}

## retrieve the cost of the test
df_c_tests <- read.xlsx("Cost tests.xlsx", sheetIndex = 1, header = TRUE) %>% 
  dplyr::mutate('burkina-faso' = burkina.faso,
                'Tests in scope of trial' = Tests.in.scope.of.trial) %>%
  dplyr::select(-c(burkina.faso, Tests.in.scope.of.trial)) %>%
  pivot_longer(
    cols = c("uganda":"burkina-faso"),
    names_to = "country_t",
    values_to = "cost"
  )



## Run decision tree for 6 scenarios
scenarios <- c(1, 2, 3, 4, 5, 6)

df_dectree_uganda <- map(.x = scenarios, 
                            ~run_dec_tree("uganda", df_uganda, scenario = .x, prob_i = FALSE, c_tests_exc_i = df_c_tests))

df_dectree_ghana <- map(.x = scenarios, 
                         ~run_dec_tree("ghana", df_ghana, scenario = .x, prob_i = FALSE, c_tests_exc_i = df_c_tests))

df_dectree_burkina_faso <- map(.x = scenarios, 
                         ~run_dec_tree("burkina-faso", df_burkina, scenario = .x, prob_i = FALSE, c_tests_exc_i = df_c_tests))


## Calculate different outcomes for each of the scenarios
df_results <- function(country_data, variant, country_loc){
  
  ## country_data   : country specific results of the decision tree
  ## variant        : outcome variable to calculate (cost, abx or DALY)
  ## country_loc    : country to run the model for
  
  if (variant == "cost"){
    
    ## Create an empty table
    df_out <- tibble(
      Strategy = c("Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All", "Standard-of-Care"),
      Cost_test = rep(0,7),
      Cost_tx = rep(0,7),
      Cost_labour = rep(0,7),
      Cost_tx_oop = rep(0,7),
      Cost_nonfav = rep(0,7),
      Cost_training = rep(0,7),
      Cost_equipment = rep(0,7)
    )
    
    for (j in 1:7){
      
      k <- base::ifelse(j == 7, 6, j)
      n <- base::ifelse(j == 7, 1, 2)
      
      ## Fill the cost dataframe
      df_out$Cost_test[[j]] <- country_data[[k]][["c_dx_total"]][[n]]
      df_out$Cost_tx[[j]] <- country_data[[k]][["c_tx_total"]][[n]]
      df_out$Cost_labour[[j]] <- country_data[[k]][["c_labour_total"]][[n]]
      df_out$Cost_tx_oop[[j]] <- country_data[[k]][["c_oop_add_abx"]][[n]]
      df_out$Cost_nonfav[[j]] <- country_data[[k]][["c_nonfav"]][[n]]
      df_out$Cost_training[[j]] <- country_data[[k]][["c_training"]][[n]] 
      df_out$Cost_equipment[[j]] <- country_data[[k]][["c_equipment"]][[n]]
      
    }
    
  ## Add up all costs and calculate cost per patient  
  df_out <- df_out %>%
    pivot_longer(cols = c("Cost_test":"Cost_equipment"), values_to = "Cost", names_to = "Cost_variable") %>%
    mutate(country = str_c(country_loc),
           Cost = Cost / 1000)
  
  }else if(variant == "abx"){
    
    ## Create an empty table
    df_out <- tibble(
      Strategy = c("Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All", "Standard-of-Care"),
      Cost = rep(0,7),
      Abx_presc_rate = rep(0,7),
      country = str_c(country_loc))
    
    for (i in 1:7){
      
      k <- base::ifelse(i == 7, 6, i)
      n <- base::ifelse(i == 7, 1, 2)
      
      ## Calculate the cost and abx prescription rate per patient
      df_out$Cost[[i]] <- country_data[[k]][["c_total"]][[n]]/1000
      df_out$Abx_presc_rate[[i]] <- country_data[[k]][["abx_presc_rate"]][[n]]/1000
      
    }
    
    
    }else if(variant == "DALY"){
      
      ## Create an empty table
      df_out <- tibble(
        Strategy = c( "Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All", "Standard-of-Care"),
        Cost = rep(0,7),
        DALY = rep(0,7),
        country = str_c(country_loc))
      
      for (i in 1:7){
        
        k <- base::ifelse(i == 7, 6, i)
        n <- base::ifelse(i == 7, 1, 2)
        
        ## Calculate the cost and DALY per patient
        df_out$Cost[[i]] <- country_data[[k]][["c_total"]][[n]]/1000
        df_out$DALY[[i]] <- country_data[[k]][["daly"]][[n]]/1000
        
      }
 
  }
  
  return(df_out)
}


## Calculate the cost outcomes
df_result_u <- df_results(df_dectree_uganda, "cost", "uganda") 
df_result_g <- df_results(df_dectree_ghana, "cost", "ghana") 
df_result_b <- df_results(df_dectree_burkina_faso, "cost", "burkina-faso")

## Combine the cost outcomes and prepare for figure 3
df_result_tot <- rbind(df_result_u, df_result_g, df_result_b) %>%
  group_by(country, Strategy) %>%
  dplyr::mutate(Strategy = factor(Strategy, levels = c( "Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All", "Standard-of-Care")),
                  Cost_variable = factor(Cost_variable, levels = c("Cost_labour", "Cost_tx_oop", "Cost_nonfav", 
                                                  "Cost_training", "Cost_equipment", 
                                                  "Cost_tx", "Cost_test")),
                Cost_variable = recode(Cost_variable, Cost_labour = 'Labour', Cost_tx_oop = "Therapeutics OOP", Cost_nonfav = "Unfavourable outcome", 
                                       Cost_training = "Training", Cost_equipment = "Equipment", 
                                       Cost_tx = "Therapeutics", Cost_test = "Tests"),
                country = recode(country, 'burkina-faso' = "Burkina Faso", ghana = "Ghana", uganda = "Uganda")) %>%
  dplyr::rename(., "Cost category" = Cost_variable)

## Color blind friendly palet
cbPal <- c("#009E73","#56B4E9", "#E69F00" ,"#999999",  "#CC79A7","#D55E00", "#0072B2", "#F0E442")

## Figure 3 manuscript 
ggplot(data = df_result_tot, aes(x = Strategy, y = Cost, fill = `Cost category`)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Cost in International $") +
  xlab("Strategy") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_text(vjust = 3.5),
        axis.title.x = element_text(vjust = -1.1)) +
  scale_fill_manual(values=cbPal)+
  facet_grid(. ~ country) +
  scale_y_continuous(labels = scales::dollar_format())


## Prepare the cost data
c_table <- df_result_tot %>%
  group_by(country, Strategy) %>%
  dplyr::mutate(total = sum(Cost),
                Cost = round(Cost, 2),
                Percentage = round(Cost/total * 100,1)) %>%
  dplyr::select(-total) %>%
  pivot_wider(names_from = 'Cost category', values_from = c(Cost, Percentage)) %>%
  rename(Country = "country") 

# Table A6 appendix
c_table_out <- gt(c_table,
              rowname_col = "sub_group",
              groupname_col = "Country") %>%
  cols_merge_n_pct(., col_n = c("Cost_Tests"), col_pct = c("Percentage_Tests")) %>%
  cols_merge_n_pct(., col_n = c("Cost_Therapeutics"), col_pct = c("Percentage_Therapeutics")) %>%
  cols_merge_n_pct(., col_n = c("Cost_Labour"), col_pct = c("Percentage_Labour")) %>%
  cols_merge_n_pct(., col_n = c("Cost_Therapeutics OOP"), col_pct = c("Percentage_Therapeutics OOP")) %>%
  cols_merge_n_pct(., col_n = c("Cost_Unfavourable outcome"), col_pct = c("Percentage_Unfavourable outcome")) %>%
  cols_merge_n_pct(., col_n = c("Cost_Training"), col_pct = c("Percentage_Training")) %>%
  cols_merge_n_pct(., col_n = c("Cost_Equipment"), col_pct = c("Percentage_Equipment")) %>%
  cols_label(
    Cost_Tests = ("Test"),
    Cost_Therapeutics = ("Therapeutics"),
    Cost_Labour = ("Labour"),
    "Cost_Therapeutics OOP" = ("Therapeutics OOP"),
    "Cost_Unfavourable outcome" = ("Non-favorable outcome"),
    Cost_Training = ("Training"),
    Cost_Equipment = ("Equipment")) %>%
  tab_spanner(
    label = "Costs, I$ (%)",
    columns = c("Cost_Tests":"Cost_Equipment")
  ) 


## Function to calculate the different outcome measures
df_sub_strat <- function(country_data, country_loc, PSA = FALSE){
  
  ## country_data   : country specific results of the decision tree
  ## country_loc    : country to run the model for
  ## PSA            : using probabilistic results or not
  
  ## Create empty table
  df_scen_group <- as_tibble(c(scenarios, 7)) %>%
    dplyr::mutate(abx = 0,
                  DDD = 0,
                  DOT = 0,
                  DALY = 0,
                  Cost = 0,
                  Abx_presc_rate = 0,
                  Country = country_loc)
  
  if(PSA == FALSE){
    
  for (i in 1:7){
    
    j <- base::ifelse(i == 7, 1, 2)
    h <- base::ifelse(i == 7, 6, i)
    
    if(length(country_data[[h]]) > 1){
      
      ## Fill dataframe with country specific data and calculate per patient values
      df_scen_group$abx[i] <- country_data[[h]][["n_abx"]][[j]]/1000
      df_scen_group$DDD[i] <- country_data[[h]][["ddd_total"]][[j]]/1000
      df_scen_group$DOT[i] <- country_data[[h]][["dot_total"]][[j]]/1000
      df_scen_group$DALY[i] <- country_data[[h]][["daly"]][[j]]/1000
      df_scen_group$Cost[i] <- country_data[[h]][["c_total"]][[j]]/1000
      df_scen_group$Abx_presc_rate[i] <- country_data[[h]][["abx_presc_rate"]][[j]]/1000
      
    }else{
      
      df_scen_group$abx[i] <- 0
      df_scen_group$DDD[i] <- 0
      df_scen_group$DOT[i] <- 0
      df_scen_group$DALY[i] <- 0
      df_scen_group$Cost[i] <- 0
      df_scen_group$Abx_presc_rate[i] <- 0
      
    }  
  
  }
    
  }else{
    
    for (i in 1:7){
      
      j <- base::ifelse(i == 7, "con", "int")
      h <- base::ifelse(i == 7, 6, i)
      
      if(length(country_data[[h]]) > 1){
        
        ## Fill dataframe with country specific data and calculate per patient values based on PSA results
        df_scen_group$abx[i] <- list(lapply(country_data[[h]][[j]][["n_abx"]], "/", 1000))
        df_scen_group$DDD[i] <- list(lapply(country_data[[h]][[j]][["ddd_total"]], "/", 1000))
        df_scen_group$DOT[i] <- list(lapply(country_data[[h]][[j]][["dot_total"]], "/", 1000))
        df_scen_group$DALY[i] <- list(lapply(country_data[[h]][[j]][["daly"]], "/", 1000))
        df_scen_group$Cost[i] <- list(lapply(country_data[[h]][[j]][["c_total"]], "/", 1000))
        df_scen_group$Abx_presc_rate[i] <- list(lapply(country_data[[h]][[j]][["abx_presc_rate"]], "/", 1000))
        
      }else{
        
        df_scen_group$abx[i] <- list(rep(0, 10000))
        df_scen_group$DDD[i] <- list(rep(0, 10000))
        df_scen_group$DOT[i] <- list(rep(0, 10000))
        df_scen_group$DALY[i] <- list(rep(0, 10000))
        df_scen_group$Cost[i] <- list(rep(0, 10000))
        df_scen_group$Abx_presc_rate[i] <- list(rep(0, 10000))
      }  
      
    }
    
  }
  
  colnames(df_scen_group) <- c("scenario", "abx", "DDD", "DOT", "DALY", "Cost", "Abx_presc_rate", "Country")
  
  return(df_scen_group)
    
}

## Create overview of all outcome variables per scenario - Deterministic
df_strategy_burkina <- df_sub_strat(df_dectree_burkina_faso, "burkina-faso", PSA = FALSE)
df_strategy_uganda <- df_sub_strat(df_dectree_uganda, "uganda", PSA = FALSE)
df_strategy_ghana <- df_sub_strat(df_dectree_ghana, "ghana", PSA = FALSE)

## Prepare for figure 4 in the manuscript
df_abx_fig_f <- rbind(df_strategy_burkina, df_strategy_uganda, df_strategy_ghana) %>%
  dplyr::mutate(Scenario = factor(scenario, levels = c(1,2,3,4,5,6,7), labels = c( "Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All", "Standard-of-Care")),
                Country = base::ifelse(Country == "ghana", "Ghana",
                                       base::ifelse(Country == "uganda", "Uganda",
                                                    "Burkina Faso"))) %>%
  select(-scenario, -DDD, -DOT, -DALY, -Cost, -abx)

## Color blind friendly palet
cbPal2 <- c("#009E73","#56B4E9", "#E69F00" ,"#999999",  "#D55E00", "#0072B2", "#CC79A7","#F0E442")

## Figure 4 Plot without error bars
ggplot(df_abx_fig_f, aes(fill = Scenario, x = 1, y = Abx_presc_rate)) +
  geom_bar(position="dodge", stat="identity") +
  facet_grid(. ~ Country) +
  scale_fill_manual(values=cbPal2) + 
  theme_bw() 

## Calculate Credible intervals based on PSA
df_dectree_uganda_psa <- map(.x = scenarios, 
                         ~run_dec_tree("uganda", df_uganda, scenario = .x, prob_i = TRUE, c_tests_exc_i = df_c_tests))

df_dectree_ghana_psa <- map(.x = scenarios, 
                        ~run_dec_tree("ghana", df_ghana, scenario = .x, prob_i = TRUE, c_tests_exc_i = df_c_tests))

df_dectree_burkina_faso_psa <- map(.x = scenarios, 
                               ~run_dec_tree("burkina-faso", df_burkina, scenario = .x, prob_i = TRUE, c_tests_exc_i = df_c_tests))

## Create overview of all outcome variables per scenario - Probabilistic
df_strategy_burkina_psa <- df_sub_strat(df_dectree_burkina_faso_psa, "burkina-faso", PSA = TRUE)
df_strategy_uganda_psa <- df_sub_strat(df_dectree_uganda_psa, "uganda", PSA = TRUE)
df_strategy_ghana_psa <- df_sub_strat(df_dectree_ghana_psa, "ghana", PSA = TRUE)

## Gather the data
df_abx_fig_f_p <- rbind(df_strategy_burkina_psa, df_strategy_uganda_psa, df_strategy_ghana_psa)  %>%
  dplyr::mutate(Scenario = factor(scenario, levels = c(1,2,3,4,5,6,7), labels = c( "Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All", "Standard-of-Care")),
                Country = base::ifelse(Country == "ghana", "Ghana",
                                       base::ifelse(Country == "uganda", "Uganda",
                                                    "Burkina Faso"))) %>%
  select(-scenario) %>%
  rowwise() %>%
  dplyr::mutate(low = quantile(unlist(Abx_presc_rate), 0.025),
                high = quantile(unlist(Abx_presc_rate), 0.975)) %>%
  dplyr::select(-abx, -DDD, -DOT, -DALY, -Cost, -Abx_presc_rate) 

## Prepare for combined figure of deterministic and probabilistic data
df_abx_fig_combined <- left_join(df_abx_fig_f_p, df_abx_fig_f, by = c("Scenario" = "Scenario", "Country" = "Country")) %>%
  dplyr::mutate(across(.cols = c("Abx_presc_rate","low","high"),
                       .fns = ~base::ifelse(.x == 0, NA, .x))) %>%
  dplyr::rename(., "Strategy" = Scenario)

## Figure 4 with error bars
ggplot(subset(df_abx_fig_combined, !is.na(Abx_presc_rate)), aes(fill = Strategy, x = 1, y = Abx_presc_rate)) +
  geom_bar(position="dodge", stat="identity") +
  facet_grid(. ~ Country) +
  scale_fill_manual(values=cbPal2) + 
  theme_bw() +
  geom_errorbar(stat = "identity", position =  position_dodge(0.9), aes(ymin = low, ymax = high), width=0.4, alpha=0.5, linewidth = 0.8) +
  ylab("Antibiotics prescription rate") +
  xlab("Age group") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.y = element_text(vjust = 3.5),
        axis.title.x = element_text(vjust = -1.2)) +
  scale_y_continuous(labels = scales::percent)


## Calculate the Credible interval for all outcome measures
df_psa_incremental_vtot <- rbind(df_strategy_burkina_psa, df_strategy_uganda_psa, df_strategy_ghana_psa) %>%
  dplyr::mutate(Scenario = factor(scenario, levels = c(1,2,3,4,5,6,7), labels = c( "Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All", "Standard-of-Care")),
                Country = base::ifelse(Country == "ghana", "Ghana",
                                       base::ifelse(Country == "uganda", "Uganda",
                                                    "Burkina-Faso"))) %>%
  rowwise() %>%
  dplyr::mutate(across(.cols = c("Abx_presc_rate", "Cost"),
                       .fns = list(low = ~quantile(unlist(.x), 0.025),
                                   high = ~quantile(unlist(.x), 0.975)),
                       .names = "{.col}_{.fn}")) %>%
  dplyr::select(-abx, -DDD, -DOT, -DALY, -Cost, -Abx_presc_rate, -scenario) 

## Calculate the deterministic values for all outcome measures
df_det_incremental_vtot <- rbind(df_strategy_burkina, df_strategy_uganda, df_strategy_ghana) %>%
  dplyr::mutate(Scenario = factor(scenario, levels = c(1,2,3,4,5,6,7), labels = c( "Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All", "Standard-of-Care")),
                Country = base::ifelse(Country == "ghana", "Ghana",
                                       base::ifelse(Country == "uganda", "Uganda",
                                                    "Burkina-Faso")))

## Table 3 of manuscript
icer_tab <- left_join(df_det_incremental_vtot,  df_psa_incremental_vtot, by = c("Country" = "Country", "Scenario" = "Scenario")) 

## CEAC
## Function for calculation of incremental cost per abx percentage point reduction
incremental_calc <- function(country_input, country_loc){
    
  ## country_data : country specific results of the decision tree
  ## country_loc  : country to run the model for
  
    scenarios <- c(1,2,3,4,5,6)
    
    df_scen_group <- as_tibble(scenarios) %>%
      dplyr::mutate(Cost = 0,
                    Abx_presc_rate = 0,
                    Country = country_loc)
    
      for (i in 1:6){
        
        if(length(country_input) > 1){
          
          # Calculate the incremental cost and abx prescription rate
          df_scen_group$Cost[i] <- list(lapply(Map('-', country_input[[i]][["int"]][["c_total"]], country_input[[i]][["con"]][["c_total"]]), "/", 1000))
          df_scen_group$Abx_presc_rate[i] <- list(lapply(Map('-', country_input[[i]][["int"]][["abx_presc_rate"]], country_input[[i]][["con"]][["abx_presc_rate"]]), "/", 1000))
          
        }else{
          
          df_scen_group$Cost[i] <- list(rep(0, 10000))
          df_scen_group$Abx_presc_rate[i] <- list(rep(0, 10000))
        }  
        
      }
      
    colnames(df_scen_group) <- c("scenario", "Cost", "Abx_presc_rate", "Country", "Group")
    
    return(df_scen_group)
    
}

## Incremental results for WTP - 1000 simulations
## Calculate the Incremental results vs. Standard-of-Care for all three countries
df_wtp_burkina <- incremental_calc(df_dectree_burkina_faso_psa, "burkina-faso")
df_wtp_uganda <- incremental_calc(df_dectree_uganda_psa, "uganda")
df_wtp_ghana <- incremental_calc(df_dectree_ghana_psa, "ghana")

## Combine the incremental values 
df_wtp_full <- rbind(df_wtp_burkina, df_wtp_uganda, df_wtp_ghana) %>%
  dplyr::mutate(Scenario = factor(scenario, levels = c(1,2,3,4,5,6), labels = c( "Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All")),
                Country = base::ifelse(Country == "ghana", "Ghana",
                                       base::ifelse(Country == "uganda", "Uganda",
                                                    "Burkina Faso"))) %>%
  select(-scenario) 

## Prepare environment for calculating the CEAC
df_psa_global <- tibble()
ceac <- seq(0, 20.5, by = 0.005)
wtp <- ceac 

## Function to calculate the CEAC based on the Net Monetary Benefit
ceac_calc <- function(wtp, df_wtp_full_fil, country_loc){
  
  ## WTP              : range of WTP over which the NMB needs to be calculated
  ## df_wtp_full_fil  : difference in cost and antibiotic prescription rate per run of PSA
  ## country_loc      : country to run the model for
  
  ## Filter the input data based on the country
  df_wtp_full_fil <- df_wtp_full_fil %>%
    dplyr::filter(Country == country_loc)
  
  ## Calculate the net monetary benefit for each strategy
  NMB <- tibble(
    Trial = unlist(Map('-', lapply(df_wtp_full_fil$Abx_presc_rate[[1]],"*", -wtp*100), lapply(df_wtp_full_fil$Cost[[1]], "*", 1))),
    Malaria = unlist(Map('-', lapply(df_wtp_full_fil$Abx_presc_rate[[2]],"*",  -wtp*100), lapply(df_wtp_full_fil$Cost[[2]], "*", 1))),
    Malaria_crp = unlist(Map('-', lapply(df_wtp_full_fil$Abx_presc_rate[[3]],"*",  -wtp*100), lapply(df_wtp_full_fil$Cost[[3]], "*", 1))),
    Malaria_cbc = unlist(Map('-', lapply(df_wtp_full_fil$Abx_presc_rate[[4]],"*",  -wtp*100), lapply(df_wtp_full_fil$Cost[[4]], "*", 1))),
    Malaria_crp_cbc = unlist(Map('-', lapply(df_wtp_full_fil$Abx_presc_rate[[5]],"*",  -wtp*100), lapply(df_wtp_full_fil$Cost[[5]], "*", 1))),
    Malaria_all = unlist(Map('-', lapply(df_wtp_full_fil$Abx_presc_rate[[6]],"*",  -wtp*100), lapply(df_wtp_full_fil$Cost[[6]], "*", 1)))) %>%
    dplyr::mutate(across(.cols = everything(),
                         .fns = ~base::ifelse(.x <= 0, 0,1),
                         .names = "{col}N"),
                  SOC = rowSums(across("TrialN":"Malaria_allN")),
                  SOC = base::ifelse(SOC > 0, 0, 1)) %>%
    dplyr::select(-c("TrialN":"Malaria_allN"))
  
  ## Calculate the probability each strategy had the highest NMB per run of the PSA
  df_dominant <- apply(NMB, 1, which.max) %>%
    as.numeric() %>%
    plyr::count() %>%
    as.data.frame() %>%
    dplyr::mutate(wtp = wtp,
                  scenario = x,
                  percentage = freq/sum(freq)*100) %>%
    dplyr::select(-c("x", "freq"))
  
  return(df_dominant)
  
}

## Fill the WTP dataframe with NMB per WTP for each country
#   - Takes ~15 minutes
df_wtp_burkina <- map(.x = wtp,
               ~ceac_calc(.x, df_wtp_full, "Burkina Faso"))

df_wtp_uganda <- map(.x = wtp,
              ~ceac_calc(.x, df_wtp_full, "Uganda"))

df_wtp_ghana <- map(.x = wtp,
              ~ceac_calc(.x, df_wtp_full, "Ghana"))

## Combine the three WTP dataframes
df_wtp_total <- list(df_wtp_burkina, df_wtp_uganda, df_wtp_ghana)

## Function to link the WTP empty table with the probability of each scenario to result in the highest NMB per WTP step
ceac_list <- function(df_ceac_in, ceac_i){
  
  df_ceac <- expand.grid(
    wtp = ceac_i,
    scenario = 1:7,
    percentages = 0)


  for (i in 1:length(df_wtp_total[[1]])){
  
    df_ceac <- left_join(df_ceac, df_ceac_in[[i]], by = c("scenario" = "scenario", "wtp" = "wtp"))  %>%
    
      dplyr::mutate(percentages = base::ifelse(is.na(percentage) == FALSE, percentage, percentages)) %>%
      dplyr::select(wtp, scenario, percentages)
    }

  df_ceac <- df_ceac %>%
    dplyr::mutate(scenario = as.factor(scenario),
                scenario = base::ifelse(scenario == 1, "Trial-based", 
                                        base::ifelse(scenario == 2, "Malaria", 
                                                     base::ifelse(scenario == 3 ,"Malaria + CRP", 
                                                                  base::ifelse(scenario == 4,"Malaria + CBC", 
                                                                               base::ifelse(scenario == 5,"Malaria + CRP + CBC", 
                                                                                            base::ifelse(scenario == 6,"Malaria + All", 
                                                                                                         base::ifelse(scenario == 7, "Standard-of-Care", NA))))))))
  
  return(df_ceac)
}

## Determine for each WTP level the probabilities per scenario 
df_ceac_total <- map(.x = df_wtp_total,
                     ~ceac_list(.x, ceac_i = ceac))

## Combine the results for each country into one dataframe
df_total_wtp <- data.frame(df_ceac_total[[1]], country= "Burkina Faso") %>%
  rbind(., data.frame(df_ceac_total[[2]], country= "Uganda"), data.frame(df_ceac_total[[3]], country = "Ghana")) %>%
  rename(Strategy = "scenario") %>%
  dplyr::mutate(Strategy = factor(Strategy, levels = c( "Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All", "Standard-of-Care"), labels = c( "Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All", "Standard-of-Care"))) %>%
  mutate(percentages = percentages/100)

# Color blind friendly palet
cbPal2 <- c("#009E73","#56B4E9", "#E69F00" ,"#999999",  "#D55E00", "#0072B2", "#CC79A7","#F0E442")

## Create the CEAC plot 
ggplot(df_total_wtp, aes(x = wtp, y = percentages, group = Strategy)) +
  ggforce::facet_col(facets = vars(country)) +
  geom_line(aes(colour = Strategy), linewidth = 1) +
  theme_bw() +
  scale_color_manual(values=cbPal2)+ 
  ylim(0, 100) +
  ylab("Probability cost-effective") + 
  xlab("Willingness-to-Pay (International $) per consultation for 1 ppt reduction in antibiotic prescription rate") +
  scale_y_continuous(labels = scales::percent) + 
  scale_x_continuous(labels = scales::dollar_format(), limits = c(0,20)) 


## Scatter plot - Create empty tibble
df_scat <- tibble(
  Cost = NA,
  Abx_presc_rate = NA,
  Country = NA,
  Strategy = NA) 

for (i in 1:nrow(df_wtp_full)){
  
  ## Create a long list with cost and abx prescribtion rate to create a scatter plot
  df_scatter <- tibble(
    Cost = unlist(df_wtp_full$Cost[i]),
    Abx_presc_rate = unlist(df_wtp_full$Abx_presc_rate[i]),
    Country = df_wtp_full$Country[i],
    Strategy = df_wtp_full$Scenario[i]) 
  
  df_scat <- rbind(df_scat, df_scatter)
  
}

## Rename the scatter plot variables
df_scat <- df_scat[-1,] %>%
  dplyr::mutate(Strategy = factor(Strategy, levels = c( "Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All"), 
                                  labels = c( "Trial-based", "Malaria", "Malaria + CRP", "Malaria + CBC", "Malaria + CRP + CBC", "Malaria + All")))

## Create the scatter plot
ggplot(df_scat, aes(x = -Abx_presc_rate*100, y = Cost, col = Strategy)) +
  facet_grid(. ~Country) +
  geom_point(size = 1, alpha = 0.2) +
  theme_bw() +
  scale_color_manual(values=cbPal2)+ 
  xlim(-20, 65) + 
  ylab("Incremental Cost (in International $)") + 
  xlab("Incremental percentage point difference in antibiotic prescription rate") + 
  theme(axis.title.y = element_text(vjust = 3.5),
        axis.title.x = element_text(vjust = -1.2)) +
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  scale_y_continuous(labels = scales::dollar_format(), limits = c(-21, 61)) +
  guides(color = guide_legend(override.aes = list(alpha = 1,
                                                  size = 3) ) )
