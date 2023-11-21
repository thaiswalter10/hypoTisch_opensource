library(targets)
library(here)
source(here("code", "packages.R"))
tar_option_set(packages = 
                 c("tidyverse",
                   "readxl",
                   "writexl",
                   "cowplot",
                   "hms",
                   "rlang",
                   "ggsignif",
                   "lubridate"))

source(here("code", "functions_analyze.R"))
source(here("code", "functions_visualization.R"))

list(
################
### Analyses ###
################
  tar_target(flow_chart,
             do_flow_chart(awake_mandonnet_raw, gliomes_nomoni)
             ),
  tar_target(table_missing_demog,
             do_missing_demog_table(demog)
             ),
  tar_target(table_missing_pec,
            do_table_missing_pec(gliomes_awake_1121, analyze_df, demog)
            ),
  tar_target(table1,
             do_table1(demog, rounding_precision = 1)
             ),
  tar_target(table_raw_pec,
           do_table_raw_pec(gliomes_awake_1121, analyze_df, demog)
            ),
  tar_target(table_pec,
           do_table_pec(table_raw_pec, rounding_precision = 1)
            ),
  tar_target(analyze_df,
             do_analyze_df(gliomes_awake_1121, moni)
             ),
  tar_target(gradeIV,
             gliomes_awake_1121 %>% filter(grade == "GBM")
             ),
  tar_target(analyze_dfIV,
             do_analyze_df(gradeIV, moni)
             ),
  tar_target(univariate_glm,
             do_univariate_glm(analyze_df,
                                vars = c("PAS_under_90p_mins", "PAM_mean", "PAD_mean", 
                                         "bilan_entree_sortie_normalized", "volume_exp_normalized", "mannitol",
                                         "diuresis_normalized", "blood_loss", "catheco", "moni_dur_dbl", "perc_resection"
                                         ),
                                alpha = 0.975, 
                                p_precision=4, 
                                ci_precision=2)
            ),
  tar_target(univariate_glmIV,
             do_univariate_glm(analyze_dfIV,
                               vars = c("PAS_under_90p_mins", "PAM_mean", "PAD_mean", 
                                        "bilan_entree_sortie_normalized", "volume_exp_normalized", "mannitol",
                                        "diuresis_normalized", "blood_loss", "catheco", "moni_dur_dbl", "perc_resection"
                               ),
                               alpha = 0.975, 
                               p_precision=4, 
                               ci_precision=2)
             ),
  tar_target(multivariate_glm,
           do_multivariate_glm(analyze_df, 
                            univariate_df = univariate_glm, 
                            alpha = 0.975, 
                            p_precision = 3,
                            ci_precision = 4)
            ),
  tar_target(multivariate_glmIV,
             do_univariate_glm(analyze_dfIV,
                               vars = c("PAS_under_90p_mins", "PAM_mean", "PAD_mean", 
                                        "bilan_entree_sortie_normalized", "volume_exp_normalized", "mannitol",
                                        "diuresis_normalized", "blood_loss", "catheco", "moni_dur_dbl", "perc_resection"
                               ),
                               alpha = 0.975, 
                               p_precision=4, 
                               ci_precision=2)
             ),
  tar_target(forest_plot,
           do_forest_plot(univariate_glm, 
                          multivariate_glm)
           ),
  tar_target(univariate_kendall,
             do_univariate_kendall(analyze_df,
                                  vars = c("preop_vol", "PAS_under_90p_mins", "PAM_mean", "PAD_mean", "moni_dur_dbl", 
                                           "diuresis_normalized", "volume_exp_normalized", "mannitol",
                                           "blood_loss", "perc_resection",
                                           "bilan_entree_sortie_normalized"))
          ),
  tar_target(univariate_kendallIV,
             do_univariate_kendall(analyze_dfIV,
                                   vars = c("preop_vol", "PAS_under_90p_mins", "PAM_mean", "PAD_mean", "moni_dur_dbl", 
                                            "diuresis_normalized", "volume_exp_normalized", "mannitol",
                                            "blood_loss", "perc_resection",
                                            "bilan_entree_sortie_normalized"))
             ),
  tar_target(tau_per_bootstrap, 
             do_tau_per_bootstrap(analyze_df, 
                                  univariate_kendall)
             ),
  tar_target(kendall_table, 
             do_kendall_table(analyze_df,
                              univariate_kendall)
           ),
  tar_target(kendall_tableIV,
             do_kendall_table(analyze_dfIV,
                              univariate_kendallIV)
             ),
  tar_target(hist_tau,
             draw_tau_per_bootstrap(tau_per_bootstrap)
             )
)