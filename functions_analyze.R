
under_fixed_cutoff <- function(moni, cut_off){
  moni %>% 
    group_by(surgery_date) %>%
    mutate("PAM_{{cut_off}}" := PAM < cut_off) %>%
    ungroup()
}

under_cutoff <- function(moni, cut_off){
  # cut_off is a number between 0 and 100, eg: 90 / 80 ... 
  ref_values <- moni %>% 
    group_by(surgery_date) %>%
    filter(row_number()==1) %>%
    select(surgery_date, PAS, PAM, PAD, ppulsee) %>%
    rename(ref_PAS = PAS, 
           ref_PAM = PAM, 
           ref_PAD = PAD,
           ref_ppulsee = ppulsee)
  
  moni %>%
    left_join(ref_values, by = "surgery_date") %>%
    mutate("PAS_under_{{cut_off}}" := PAS < cut_off * ref_PAS/100,
           "PAM_under_{{cut_off}}" := PAM < cut_off * ref_PAM/100,
           "PAD_under_{{cut_off}}" := PAD < cut_off * ref_PAD/100,
           "ppulsee_under_{{cut_off}}" := ppulsee < cut_off * ref_ppulsee/100) %>%
    group_by(name, nipp, surgery_date) %>%
    summarise(across(starts_with(c("PAS_under", "PAM_under", "PAD_under", "ppulsee_under")), ~sum(.x, na.rm = TRUE)*5, .names = "{.col}p_mins")) %>%
    ungroup()
}

variability <- function(x){
  sd(x, na.rm = TRUE) /mean(x, na.rm = TRUE)
}

summarize_vars <- function(data, vars){
  data %>% 
    group_by(name, nipp, surgery_date) %>%
    summarise(
      across({{vars}}, variability, .names = "{.col}_var"),
      across({{vars}}, ~min(.x, na.rm = TRUE), .names = "{.col}_min"),
      across({{vars}}, ~max(.x, na.rm = TRUE), .names = "{.col}_max"),
      across({{vars}}, ~median(.x, na.rm = TRUE), .names = "{.col}_median"),
      across({{vars}}, ~mean(.x, na.rm = TRUE), .names = "{.col}_mean")) %>%
    ungroup()
}

vect_wilcox <- function(df, grouping_variable, x, y){
  # x & y are characters
  # y is the outcome
  # grouping_variable is an object
  df %>% 
    group_by({{grouping_variable}}) %>%  
    group_modify(~ broom::tidy(wilcox.test(pull(.x[x]) ~ pull(.x[y])))) %>%
    ungroup %>%
    mutate(p.value = round(p.value, 3)) %>%
    #mutate(p.value = ifelse(p.value <= 0.01, "< 0.01", as.character(p.value))) %>%
    select({{grouping_variable}}, p.value)
}

vect_kendall <- function(df, grouping_variable, x, y){
  # x & y are characters
  # y is the outcome
  # grouping_variable is an object name
  df %>% 
    group_by({{grouping_variable}}) %>%  
    group_modify(~ broom::tidy(cor.test(pull(.x[x]), pull(.x[y]), method = "kendall"))) %>%
    ungroup()
}

kendall_for_bootstrap <- function(split, x, y) {
  dat <- analysis(split)
  tidy(cor.test(pull(dat[x]), pull(dat[y]), method = "kendall")) %>%
    mutate(term = "cor")
}

vect_lme <- function (df_time_point_as_col, var){
  # give var as a character
  df_time_point_as_col %>% 
    group_by(name_var) %>%
    group_modify(
      ~broom::tidy(
        anova(lme(eval(substitute(j ~ 1, list(j = as.name({{var}})))), random = ~1 |
                    gupi, na.action = na.omit, data = .))
      )
    ) %>%
    ungroup() %>%
    filter(column == "p-value") %>% 
    select(name_var, column, mean) %>%
    pivot_wider(names_from = column, values_from = mean)
}

desc_continuous_var <- function(df, list_continuous_var, rounding_precision){
  # df has to contain the variable "postop_ischemia"
  # list_continuous_var is a list of strings of the names of the variable which
  # descriptive statistics have to be computed. 
  cont <- df %>%
    summarise(across(all_of(list_continuous_var), 
                     list(median = ~round(median(.x, na.rm = TRUE), rounding_precision),
                          q1 = ~round(quantile(.x, 0.25, na.rm = TRUE), rounding_precision),
                          q3 = ~round(quantile(.x, 0.75, na.rm = TRUE), rounding_precision)))) %>%
    pivot_longer(cols = -postop_ischemia,
                 names_to = c("variable", "computation"),
                 names_pattern = "(.*)_(.*)$",
                 values_to = c("value")) %>%
    mutate(postop_ischemia = ifelse(postop_ischemia == TRUE, "isch", "no_isch")) %>%
    pivot_wider(names_from = c(computation, postop_ischemia), 
                values_from = value) %>%
    mutate(isch =  str_glue("{median_isch} [{q1_isch}-{q3_isch}]"),
           no_isch = str_glue("{median_no_isch} [{q1_no_isch}-{q3_no_isch}]")) %>%
    select(-starts_with(c("median", "q")))
}

desc_continuous_var_whole_cohort <- function(df, list_continuous_var, rounding_precision){
  # df has to contain the variable "postop_ischemia"
  # list_continuous_var is a list of strings of the names of the variable which
  # descriptive statistics have to be computed. 
  cont <- df %>%
    summarise(across(all_of(list_continuous_var), 
                     list(median = ~round(median(.x, na.rm = TRUE), rounding_precision),
                          q1 = ~round(quantile(.x, 0.25, na.rm = TRUE), rounding_precision),
                          q3 = ~round(quantile(.x, 0.75, na.rm = TRUE), rounding_precision)))) %>%
    pivot_longer(cols = everything(),
                 names_to = c("variable", "computation"),
                 names_pattern = "(.*)_(.*)$",
                 values_to = c("value")) %>%
    pivot_wider(names_from = c(computation), 
                values_from = value) %>%
    mutate(whole_cohort =  str_glue("{median} [{q1}-{q3}]")) %>%
    select(-starts_with(c("median", "q")))
}

desc_binary_variable <- function(df, list_binary_var, rounding_precision){
  # df has to contain the postop_ischemia variable
  # list_binary_var is a list of strings of the names of the variable which
  # descriptive statistics have to be computed. 
  tot <- dim(df)[1]
  bin <- df %>%
    summarise(across(
      all_of(list_binary_var), 
      list(n = ~round(sum(.x, na.rm= TRUE), rounding_precision), 
           perc = ~round(sum(.x, na.rm = TRUE)/tot*100, rounding_precision)))) %>%
    pivot_longer(cols = -postop_ischemia,
                 names_to = c("variable", "computation"),
                 names_pattern = "(.*)_(.*)$",
                 values_to = c("value")) %>%
    mutate(postop_ischemia = ifelse(postop_ischemia == TRUE, "isch", "no_isch")) %>%
    pivot_wider(names_from = c(computation, postop_ischemia), 
                values_from = value) %>%
    mutate(isch = str_glue("{n_isch} ({perc_isch})"),
           no_isch = str_glue("{n_no_isch} ({perc_no_isch})")) %>%
    select(-starts_with(c("n_", "perc")))
}

desc_binary_variable_whole_cohort <- function(df, list_binary_var, rounding_precision){
  # df has to contain the postop_ischemia variable
  # list_binary_var is a list of strings of the names of the variable which
  # descriptive statistics have to be computed. 
  tot <- dim(df)[1]
  bin <- df %>%
    summarise(across(
      all_of(list_binary_var), 
      list(n = ~round(sum(.x, na.rm= TRUE), rounding_precision), 
           perc = ~round(sum(.x, na.rm = TRUE)/tot*100, rounding_precision)))) %>%
    pivot_longer(cols = everything(),
                 names_to = c("variable", "computation"),
                 names_pattern = "(.*)_(.*)$",
                 values_to = c("value")) %>%
    pivot_wider(names_from = c(computation), 
                values_from = value) %>%
    mutate(whole_cohort = str_glue("{n} ({perc})")) %>%
    select(-starts_with(c("n", "perc")))
}

confidence_interval <- function (estimate, sd, alpha = 0.975) {
  lower_limit = exp(estimate - qnorm(alpha) * sd)
  upper_limit = exp(estimate + qnorm(alpha) * sd)
  limit = c(lower_limit = lower_limit, upper_limit = upper_limit)
  return(limit)
}

do_flow_chart <- function(awake_mandonnet_raw, gliomes_nomoni){
  all_surgeries <- dim(awake_mandonnet_raw)[1]
  gliomes_only <- awake_mandonnet_raw %>% filter(gliomes == 1)
  chir_eveillee_only <- gliomes_only %>% filter(is.na(chir_eveil_1oui_0non))
  complete_postop_mri <- chir_eveillee_only %>% filter(is.na(IRM_with_all_sequences_found_by_nawel))
  complete_anesthesia_folder <- complete_postop_mri %>% filter(is.na(dossier_anesth_perdu))
  
  nb_no_gliomes <- all_surgeries - dim(gliomes_only)[1]
  nb_no_awake_surgery <- dim(gliomes_only)[1] - dim(chir_eveillee_only)[1]
  no_complete_postop_mri <- dim(chir_eveillee_only)[1] - dim(complete_postop_mri)[1]
  nb_lost_anesthesia_folder <- dim(complete_postop_mri)[1] - dim(complete_anesthesia_folder)[1]
  exluded = nb_no_gliomes + 
    nb_no_awake_surgery +
    no_complete_postop_mri + 
    nb_lost_anesthesia_folder
  included_surgeries <- dim(complete_anesthesia_folder)[1]
  included_patients <- gliomes_nomoni %>% distinct(nipp) %>% dim()
  
  tribble(
    ~var,               ~n,
    "all_surgeries", all_surgeries,
    "exluded",      exluded,
    "no gliomes",   nb_no_gliomes,
    "no_awake_surgery", nb_no_awake_surgery,
    "no_complete_postop_mri", no_complete_postop_mri,
    "lost_anesthesia_folder", nb_lost_anesthesia_folder,
    "included_surgeries", included_surgeries,
    "included_patients", included_patients[1]
  ) %>% unnest(cols = c(n))
}

do_analyze_df <- function(gliomes_awake_1121, moni){
  surgery_dates <- gliomes_awake_1121 %>% 
    pull(surgery_date)
  
  df_demog <- gliomes_awake_1121 %>%
    select(surgery_date, preop_vol, diuresis_normalized, bilan_entree_sortie, 
           bilan_entree_sortie_normalized, volume_exp_normalized, perc_resection, blood_loss, postop_vol, sex, hta, moni_dur, 
           moni_dur_dbl, age, high_grade, grade, catheco, postop_ischemia, volume_isch_cc, preop_creat, preop_hb,
           ends_with(c("_test", "_covar", "_w", "_sympt")), mannitol)
  
  moni_outc_raw <- moni %>% 
    filter(surgery_date %in% surgery_dates) %>%
    filter(name != "GAUCHE")
  
  df_cutoff_90 <- moni_outc_raw %>%
    mutate(ppulsee = PAS-PAD) %>%
    under_cutoff(90) %>%
    select(-nipp, -name)
  df_cutoff_80 <- moni_outc_raw %>%
    mutate(ppulsee = PAS-PAD) %>%
    under_cutoff(80) %>%
    select(-nipp, -name)
  df_cutoff_70 <- moni_outc_raw %>%
    mutate(ppulsee = PAS-PAD) %>%
    under_cutoff(70) %>%
    select(-nipp, -name)
  
  df_fixed_cutoff <- moni_outc_raw %>%
    under_fixed_cutoff(70) %>%
    under_fixed_cutoff(65) %>%
    under_fixed_cutoff(60) %>%
    select(-nipp, -name) %>% 
    group_by(surgery_date) %>% 
    summarise(PAM_70 = sum(PAM_70)*5, 
              PAM_65 = sum(PAM_65)*5, 
              PAM_60 = sum(PAM_60)*5)
  
  moni_outc <- moni_outc_raw %>%
    select(-origin) %>%
    mutate(ppulsee = PAS-PAD) %>%
    summarize_vars(c(FC, PAS, PAM, PAD, ppulsee))  %>%
    left_join(df_cutoff_90, by = ("surgery_date")) %>%
    left_join(df_cutoff_80, by = ("surgery_date")) %>%
    left_join(df_cutoff_70, by = ("surgery_date")) %>%
    left_join(df_demog, by = "surgery_date") %>%
    left_join(df_fixed_cutoff, by = "surgery_date")
}

do_missing_demog_table <- function(demog){
  demog %>% 
    summarise(across(everything(), ~sum(is.na(.x)))) %>%
    pivot_longer(cols = everything(), names_to = "variables", values_to = "number_missing_var")
}

do_table_missing_pec <- function(gliomes_awake_1121, analyze_df, demog){
  first_bp <- demog %>% select(starts_with("first"), surgery_date)
  hd <- analyze_df %>% select(PAS_under_90p_mins, PAS_mean, PAM_mean, PAD_mean, surgery_date)
  
  raw_df <- gliomes_awake_1121 %>%
    left_join(first_bp, by = "surgery_date") %>%
    left_join(hd, by = "surgery_date") %>%
    select(moni_dur_dbl, bilan_entree_sortie_normalized,  diuresis_normalized,
           volume_exp_normalized, red_blood_cell, blood_loss, catheco,
           postop_vol, perc_resection, postop_ischemia, volume_isch_cc, 
           starts_with("first"), PAS_mean, PAM_mean, PAD_mean, PAS_under_90p_mins, 
           mannitol) %>%
    select(-first_name)
  
  missing_pec <- raw_df %>% 
    summarise(across(everything(), ~sum(is.na(.x)))) %>%
    pivot_longer(cols = everything(), names_to = "variables", values_to = "number_missing_var") %>%
    arrange(desc(number_missing_var))
}

do_table1 <- function(demog, rounding_precision=1){
  
  raw_table1 <- demog %>% 
    select(-nipp, -surgery_date, -dob, -name, -volume_isch_cc, -med_hist_dex, -tbc_sevre) %>%
    rename(sex_feminine = sex) %>%
    mutate(sex_ckdepi = sex_feminine == 0,
           creat_ckdepi = preop_creat/100) %>%
    mutate(ckd_epi_dfg = CKDEpi.creat.rf(creat_ckdepi, sex_ckdepi, age)) %>%
    select(-sex_ckdepi) %>%
    group_by(postop_ischemia) 
  
  whole_cohort <- raw_table1 %>% ungroup()
  
  colname_isch <- str_glue("Patients with postoperative ischemia (n={dim(raw_table1 %>% filter(postop_ischemia == TRUE))[1]})")
  colname_no_isch <- str_glue("Patients without postoperative ischemia (n={dim(raw_table1 %>% filter(postop_ischemia == FALSE))[1]})")
  
  cont_var <- c("age", "bmi", "asa", "preop_vol", "preop_hb", "preop_creat", "ckd_epi_dfg")
  cont <- desc_continuous_var(raw_table1, cont_var, rounding_precision = 1)
  cont_whole <- desc_continuous_var_whole_cohort(whole_cohort, cont_var, rounding_precision = 1)
  
  bin_var <- c("sex_feminine", "intox_tabac", "hta", "epilepsie_2R", "htic", "high_grade", "redo")
  bin <- desc_binary_variable(raw_table1, bin_var, rounding_precision = 1)
  bin_whole <- desc_binary_variable_whole_cohort(whole_cohort, bin_var, rounding_precision = 1)
  
  var_whole <- cont_whole %>% add_row(bin_whole)
  
  ordered_raw <- c("age", "sex_feminine", "bmi", "asa", "intox_tabac", "hta", 
                   "preop_vol", "high_grade", "redo", "epilepsie_2R", "htic", 
                   "preop_hb", "preop_creat", "ckd_epi_dfg")
  
  pvalue <- raw_table1 %>%
    ungroup() %>%
    select(c(all_of(ordered_raw), postop_ischemia)) %>% 
    pivot_longer(all_of(ordered_raw), 
                 names_to = "variable",
                 values_to = "value") %>%
    vect_wilcox(variable, "value", "postop_ischemia") %>%
    mutate("p value" = round(p.value, 2)) %>%
    select(-p.value)
  
  table1 <- cont %>% 
    add_row(bin) %>%
    left_join(var_whole, by = "variable") %>%
    mutate(variable = factor(variable, levels = ordered_raw)) %>%
    arrange(variable) %>%
    rename({{colname_isch}} := isch,
           {{colname_no_isch}} := no_isch) %>%
    left_join(pvalue) %>%
    relocate(whole_cohort, .after = variable)
}

do_table_raw_pec <- function(gliomes_awake_1121, analyze_df, demog){
  first_bp <- demog %>% select(starts_with("first"), surgery_date)
  hd <- analyze_df %>% select(PAS_under_90p_mins, PAS_mean, PAM_mean, PAD_mean, surgery_date)
  
  raw_df <- gliomes_awake_1121 %>%
    left_join(first_bp, by = "surgery_date") %>%
    left_join(hd, by = "surgery_date") %>%
    select(moni_dur_dbl, bilan_entree_sortie_normalized,  diuresis_normalized,
           volume_exp_normalized, red_blood_cell, blood_loss, catheco,
           postop_vol, perc_resection, postop_ischemia, volume_isch_cc, 
           starts_with("first"), PAS_mean, PAM_mean, PAD_mean, PAS_under_90p_mins, 
           mannitol_bin) %>%
    group_by(postop_ischemia) 
}

do_table_pec <- function(table_raw_pec, rounding_precision = 1){
  whole_cohort <- table_raw_pec %>% ungroup()
  
  #Formatting
  colname_isch <- str_glue("Surgical operations with postoperative ischemia (n={dim(table_raw_pec %>% filter(postop_ischemia == TRUE))[1]})")
  colname_no_isch <- str_glue("Surgical operations without postoperative ischemia (n={dim(table_raw_pec %>% filter(postop_ischemia == FALSE))[1]})")
  
  # Descriptive statistics
  # Continous variables
  cont_var <- c("first_PAS", "first_PAM", "first_PAD", 
                "PAS_mean", "PAM_mean", "PAD_mean", "PAS_under_90p_mins", "moni_dur_dbl", "volume_exp_normalized", 
                "diuresis_normalized", "blood_loss", "bilan_entree_sortie_normalized", 
                "postop_vol", "perc_resection", "volume_isch_cc")
  cont <- desc_continuous_var(table_raw_pec, cont_var, rounding_precision = 1)
  cont_whole <- desc_continuous_var_whole_cohort(whole_cohort, cont_var, rounding_precision = 1)
  
  # Binary variables
  bin_var <- c("catheco", "mannitol_bin")
  bin <- desc_binary_variable(table_raw_pec, bin_var, rounding_precision = 1)
  bin_whole <- desc_binary_variable_whole_cohort(whole_cohort, bin_var, rounding_precision = 1)
  
  var_whole <- cont_whole %>% add_row(bin_whole)
  
  ordered_raw <- c("PAS_under_90p_mins","bilan_entree_sortie_normalized", "volume_exp_normalized", "mannitol_bin",
                   "diuresis_normalized", "blood_loss",  
                   "catheco", "moni_dur_dbl", "first_PAS", "first_PAM", "first_PAD", 
                   "PAS_mean", "PAM_mean", "PAD_mean", "postop_vol", "perc_resection", "volume_isch_cc")
  
  pvalue <- table_raw_pec %>%
    ungroup() %>%
    select(c(all_of(ordered_raw), postop_ischemia)) %>% 
    pivot_longer(all_of(ordered_raw), 
                 names_to = "variable",
                 values_to = "value") %>%
    vect_wilcox(variable, "value", "postop_ischemia") %>%
    mutate("p value" = round(p.value, 2)) %>%
    select(-p.value)
  
  # Joining all tables
  table_pec <- cont %>% 
    add_row(bin) %>%
    left_join(var_whole, by = "variable") %>%
    mutate(variable = factor(variable, levels = ordered_raw)) %>%
    arrange(variable) %>%
    left_join(pvalue, by = "variable") %>%
    rename({{colname_isch}} := isch,
           {{colname_no_isch}} := no_isch) %>%
    relocate(whole_cohort, .after = variable)
}

do_univariate_glm <- function(analyze_df, vars, alpha, p_precision, ci_precision){
  df1 <- analyze_df %>%
    pivot_longer(all_of(vars), 
                 names_to = "variable", 
                 values_to = "value") %>%
    group_by(variable) %>%
    group_modify( ~tidy(glm(postop_ischemia ~ value,
                            family="binomial", 
                            data= .))) %>%
    ungroup() %>%
    filter(term == "value") %>%
    mutate(odd_ratio = exp(estimate),
           p.value = round(p.value, p_precision),
           lower_limit_ci = round(exp(estimate - qnorm(alpha) * std.error), ci_precision),
           higher_limit_ci = round(exp(estimate + qnorm(alpha) * std.error), ci_precision)) %>%
    select(variable, odd_ratio, lower_limit_ci, higher_limit_ci, p.value) %>%
    arrange(p.value)
}

do_multivariate_glm <- function(analyze_df, univariate_df, alpha, p_precision, ci_precision){
  vars_for_multivariate <- univariate_df %>% filter(p.value < 0.05) %>% pull(variable)
  formula1 <- paste0(vars_for_multivariate, collapse=" + ")
  formula2 <- str_c(formula1,"age", "preop_vol", "high_grade", sep = " + ")
  
  tidy(glm(postop_ischemia ~ 
             diuresis_normalized + 
             PAS_under_90p_mins + 
             moni_dur_dbl + 
             preop_vol + 
             bilan_entree_sortie_normalized + 
             catheco + 
             high_grade +
             hta + 
             age,
           family="binomial", 
           data= analyze_df)) %>%
    filter(!row_number()==1) %>%
    mutate(odd_ratio = exp(estimate),
           p.value = round(p.value, p_precision),
           lower_limit_ci = round(exp(estimate - qnorm(alpha) * std.error), ci_precision),
           higher_limit_ci = round(exp(estimate + qnorm(alpha) * std.error), ci_precision)) %>%
    select(term, odd_ratio, lower_limit_ci, higher_limit_ci, p.value) %>%
    rename(variable = term) %>%
    arrange(p.value)
}

do_univariate_kendall <- function(analyze_df, vars, p_precision = 2){
  analyze_df %>%
    filter(postop_ischemia == 1) %>%
    pivot_longer(all_of(vars), 
                 names_to = "variable", 
                 values_to = "value") %>%
    vect_kendall(variable, "value", "volume_isch_cc") %>%
    mutate(p.value = round(p.value, 2)) %>%
    arrange(p.value)
}

do_tau_per_bootstrap <- function(analyze_df, univariate_kendall){
  vars <- univariate_kendall$variable
  
  estimate_per_bootstrap <- analyze_df %>% 
    select(all_of(vars), volume_isch_cc) %>%
    pivot_longer(all_of(vars), 
                 names_to = "variable", 
                 values_to = "value") %>%
    group_by(variable) %>% 
    group_modify(
      ~ bootstraps(., 5000, apparent = TRUE) %>% 
        mutate(models = map(splits, kendall_for_bootstrap, "value", "volume_isch_cc")))
}

do_kendall_table <- function(analyze_df, univariate_kendall){
  vars <- univariate_kendall$variable
  
  ci_per_var <- analyze_df %>% 
    select(all_of(vars), volume_isch_cc) %>%
    pivot_longer(all_of(vars), 
                 names_to = "variable", 
                 values_to = "value") %>%
    group_by(variable) %>% 
    group_modify(
      ~ bootstraps(., 2000, apparent = TRUE) %>% 
        mutate(models = map(splits, kendall_for_bootstrap, "value", "volume_isch_cc")) %>%
        int_pctl(models)) %>%
    select(-term, -`.estimate`, -`.alpha`, -`.method`)
  
  univariate_kendall %>% 
    select(variable, estimate, p.value) %>%
    mutate(tau = round(estimate,2)) %>%
    left_join(ci_per_var, by = "variable")  %>%
    mutate(`95% CI` = str_glue("{round(.lower,2)} - {round(.upper,2)}")) %>%
    select(-`.lower`, -`.upper`, -estimate) %>%
    relocate(p.value, .after = last_col())
}