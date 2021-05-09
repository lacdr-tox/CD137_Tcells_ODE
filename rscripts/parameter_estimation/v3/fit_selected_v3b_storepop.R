setwd("/home/richard/Documents/03_CTLsGoesNative/000_archive/rscripts/parameter_estimation/v3")
set.seed(170842342)
source(file = "fit_model_v3b.R")

Ngens <- 100
Npop <- 100
Ncore <- 16



fit_model_storepop(treatment="control",CI = FALSE,QSS = FALSE,set_growth = 0.5,save_dir = "stored_populations",save_name = "ivd7_CI-1_QSS-0_g-5_NPAR-7",fitted_pars = c("s1","ki","kr","ke","kq","dq","di"),gens = Ngens,pop=Npop,ncore=Ncore)
fit_model_storepop(treatment="mAB",CI = FALSE,QSS = FALSE,set_growth = 0.5,save_dir = "stored_populations",save_name = "ivd7_CI-1_QSS-0_g-5_NPAR-7",fitted_pars = c("s1","ki","kr","ke","kq","dq","di"),gens = Ngens,pop=Npop,ncore=Ncore)  
  
