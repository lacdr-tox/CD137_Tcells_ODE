# CD137_Tcells_ODE

Contains scripts to analyse data from: Weigelin, Bettina, et al. "Focusing and sustaining the antitumor CTL effector killer response by agonist anti-CD137 mAb." Proceedings of the National Academy of Sciences 112.24 (2015): 7551-7556. Statistics taken from imaging data containing information about CTL behaviours are analysed in the presence or absence of stimulating CD137 antibody. ODE models are developed and applied to understand differences in CTL function in these two different settings. 

Two master R markdown scripts can be used to reproduce figures. These are:

/rscripts/data_analysis.Rmd
(reproduces all figures that do not rely on ODE model) 

/rscripts/parameter_estimation/v3/ODE_MODEL_v3.rmd
(reproduces all figures that do rely on ODE model)

additional metadata and other notes about useage are contained within each script (can also be viewed in the corresponding .html files).



