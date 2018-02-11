
#create sample id
#clinical_data <- rownames_to_column(clinical_data, var = "s") #sample id in numeric
#obtain unique times
times <- sort(unique(clinical_data$os_months))

#create long data with observed and censored times
longdata <- survival::survSplit(Surv(time = os_months, event = deceased) ~ . , 
                      cut = times, data = (clinical_data %>%
                        filter(!is.na(os_status) & os_status != '') %>%
                        filter(os_months > 0 & !is.na(os_months))  %>%
  mutate(deceased = os_status == "DECEASED")))

#Checkin the correctness of the new dataset
longdata %>% group_by(sample) %>%
  select(sample, tstart, os_months, deceased) %>%
  top_n(3, tstart) %>% arrange(desc(sample, tstart))

longdata %>% group_by(sample) %>%
  select(sample, tstart, os_months, deceased) %>%
  dplyr::slice(1:3) 

assertthat::assert_that((longdata %>% group_by(sample) %>% dplyr::n_groups()) == (clinical_data %>% filter(!is.na(os_status) & os_status != '') %>% filter(os_months > 0 & !is.na(os_months)) %>% group_by(sample) %>% dplyr::n_groups()))


#Create timepoint id
longdata <- longdata %>%
  group_by(sample) %>%
  mutate(t = seq(n())) 

longdata %>% group_by(sample) %>%
  select(sample, tstart, t) %>%
  dplyr::slice(1:3) 

#Calculate the duration of each time point
longdata %>% 
  group_by(s) %>%
  mutate(t_dur = os_months - tstart) %>%
  select(s, t, tstart, os_months, t_dur) %>%
  View()

longdata <- longdata %>%
  group_by(sample) %>%
  mutate(t_dur = os_months - tstart) 

#Wrap in a function
#generate data

gen_stan_data <- function(data, formula = as.formula(~1)) {
  if(!inherits(formula, 'formula'))
  formula <- as.formula(formula)

  data <- tibble::rownames_to_column(data, var = "s")
  times <- sort(unique(data$os_months))
  
  longdata <- survival::survSplit(Surv(time = os_months, event = deceased) ~ . , 
                                  cut = times, data = (data %>%
                                                         filter(!is.na(os_status) & os_status != '') %>%
                                                         filter(os_months > 0 & !is.na(os_months))  %>%
                                                         mutate(deceased = os_status == "DECEASED")))
  #create time point id
  longdata <- longdata %>%
    group_by(sample) %>%
    mutate(t = seq(n())) 
  
  #calculate duration
  longdata <- longdata %>%
    group_by(sample) %>%
    mutate(t_dur = os_months - tstart) 
  
  stan_data <- list(
    N = nrow(longdata),
    S = dplyr::n_distinct(longdata$sample),
    "T" = length(times),
    s = array(as.numeric(longdata$s)),
    t = array(longdata$t),
    event = array(longdata$deceased),
    t_obs = array(longdata$os_months),
    t_dur = array(longdata$t_dur)
  )
}


#--- This function will take a formula object as input --- #
gen_stan_data <- function(data, formula = as.formula(~1)) {
  if(!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  observed_data <- data %>%
    dplyr::filter(os_status == "DECEASED")
  
  censored_data <- data %>%
    dplyr::filter(os_status != "DECEASED")
  
  Xobs_bg <- observed_data %>%
    model.matrix(formula, data = .)
  
  Xcen_bg <- censored_data %>%
    model.matrix(formula, data = .)
  
  assertthat::assert_that(ncol(Xcen_bg) == ncol(Xobs_bg))
  M_bg <- ncol(Xcen_bg)
  
  if (M_bg > 1){
    if("(Intercept)" %in% colnames(Xobs_bg))
      Xobs_bg <- array(Xobs_bg[,-1], dim = c(nrow(observed_data), M_bg - 1))
    if("(Intercept)" %in% colnames(Xcen_bg))
      Xcen_bg <- array(Xcen_bg[,-1], dim = c(nrow(censored_data), M_bg - 1))
    assertthat::assert_that(ncol(Xcen_bg) == ncol(Xobs_bg))
    M_bg <- ncol(Xcen_bg)
  }
  
  stan_data <- list(
    Nobs = nrow(observed_data),
    Ncen = nrow(censored_data),
    yobs = observed_data$os_months,
    ycen = censored_data$os_months,
    M_bg = M_bg,
    Xcen_bg = array(Xcen_bg, dim = c(nrow(censored_data), M_bg)),
    Xobs_bg = array(Xobs_bg, dim = c(nrow(observed_data), M_bg))
  )
}






stanfile <- "null_pem_survival_model.stan"
#open stan file
if (interactive())
  file.edit(stanfile)

pem_null_model <-  stan(stanfile,
                            data = gen_stan_data(clinical_data),
                            chains = 1,
                            iter = 5
)
