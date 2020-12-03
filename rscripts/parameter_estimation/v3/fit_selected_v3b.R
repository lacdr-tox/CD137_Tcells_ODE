setwd("/home/richard/Documents/03_CTLsGoesNative/01_draft_2/02_scripts/parameter_estimation/v3")
set.seed(170842342)
source(file = "fit_model_v3b.R")

Ngens <- 500
Npop <- 200
Ncore <- 16



for(callmewhateveryoulikebaby in 1:5){

  fit_model(treatment="control",CI = FALSE,QSS = FALSE,set_growth = 0.4,save_dir = "fitted_output",save_name = "ivd7_CI-1_QSS-0_g-4_NPAR-7",fitted_pars = c("s1","ki","kr","ke","kq","dq","di"),gens = Ngens,pop=Npop,ncore=Ncore)
  fit_model(treatment="mAB",CI = FALSE,QSS = FALSE,set_growth = 0.4,save_dir = "fitted_output",save_name = "ivd7_CI-1_QSS-0_g-4_NPAR-7",fitted_pars = c("s1","ki","kr","ke","kq","dq","di"),gens = Ngens,pop=Npop,ncore=Ncore)  
  
  fit_model(treatment="control",CI = FALSE,QSS = FALSE,set_growth = 0.5,save_dir = "fitted_output",save_name = "ivd7_CI-1_QSS-0_g-5_NPAR-7",fitted_pars = c("s1","ki","kr","ke","kq","dq","di"),gens = Ngens,pop=Npop,ncore=Ncore)
  fit_model(treatment="mAB",CI = FALSE,QSS = FALSE,set_growth = 0.5,save_dir = "fitted_output",save_name = "ivd7_CI-1_QSS-0_g-5_NPAR-7",fitted_pars = c("s1","ki","kr","ke","kq","dq","di"),gens = Ngens,pop=Npop,ncore=Ncore)  
  
  fit_model(treatment="control",CI = FALSE,QSS = FALSE,set_growth = 0.6,save_dir = "fitted_output",save_name = "ivd7_CI-1_QSS-0_g-6_NPAR-7",fitted_pars = c("s1","ki","kr","ke","kq","dq","di"),gens = Ngens,pop=Npop,ncore=Ncore)
  fit_model(treatment="mAB",CI = FALSE,QSS = FALSE,set_growth = 0.6,save_dir = "fitted_output",save_name = "ivd7_CI-1_QSS-0_g-6_NPAR-7",fitted_pars = c("s1","ki","kr","ke","kq","dq","di"),gens = Ngens,pop=Npop,ncore=Ncore)  
  
  fit_model(treatment="control",CI = FALSE,QSS = FALSE,set_growth = 0.7,save_dir = "fitted_output",save_name = "ivd7_CI-1_QSS-0_g-7_NPAR-7",fitted_pars = c("s1","ki","kr","ke","kq","dq","di"),gens = Ngens,pop=Npop,ncore=Ncore)
  fit_model(treatment="mAB",CI = FALSE,QSS = FALSE,set_growth = 0.7,save_dir = "fitted_output",save_name = "ivd7_CI-1_QSS-0_g-7_NPAR-7",fitted_pars = c("s1","ki","kr","ke","kq","dq","di"),gens = Ngens,pop=Npop,ncore=Ncore)  
  
  fit_model(treatment="control",CI = FALSE,QSS = FALSE,set_growth = 0.8,save_dir = "fitted_output",save_name = "ivd7_CI-1_QSS-0_g-8_NPAR-7",fitted_pars = c("s1","ki","kr","ke","kq","dq","di"),gens = Ngens,pop=Npop,ncore=Ncore)
  fit_model(treatment="mAB",CI = FALSE,QSS = FALSE,set_growth = 0.8,save_dir = "fitted_output",save_name = "ivd7_CI-1_QSS-0_g-8_NPAR-7",fitted_pars = c("s1","ki","kr","ke","kq","dq","di"),gens = Ngens,pop=Npop,ncore=Ncore)  
  

}