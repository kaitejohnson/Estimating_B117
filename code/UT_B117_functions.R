# Functions used to generate the outputs called in UT_B117_main.Rmd

fit_SGTF_data<-function(ut, beta_pars, spring_arrival_date, spring_end_semester_date, nsim){

n_future<-as.numeric(spring_end_semester_date-spring_arrival_date)
par_names<-colnames(beta_pars)
for (i in 1:length(par_names)){
  assign(par_names[i], as.double(beta_pars[i]))
}
## Select a date for time=0
mydate0 <- spring_arrival_date-1

## Data frame for fitting
ut_fit <- ut %>%
  mutate(time = as.numeric(test_date) - as.numeric(ymd(mydate0)))

## Dataframe for predicting in the future
ut_future <- rbind(
  ut_fit %>% mutate(period = "calibration"),
  ut_fit[rep(1, n_future), ] %>%
    mutate(test_date = max(ut_fit$test_date) + 1:n_future) %>%
    mutate(time = as.numeric(test_date) - as.numeric(ymd(mydate0))) %>%
    mutate(period = "projection")
)


###############################################################################
#  Fit a model for B117 growth at UT  #
###############################################################################


## Estimate a logistic growth model for B117 prevalence, while
## integrating over uncertainty in the B117/SGTF fraction. 
## 
## One Monte Carlo draw of this model works as follows -- (i) take a
## draw from the beta distribution described above for the fraction of
## B117/SGTF, (ii) impute B117 cases by multiplying SGTF cases by this
## draw from the beta, (iii) estimate the logistic growth model using
## this set of imputed B117 case numbers, and finally (iv) project
## future B117 prevalence using this model.  We combine all draws for
## projected B117 prevalence to integrate over uncertainty in the
## fraction of B117 to SGTF samples

## One Monte Carlo draw of this model works as follows -- (i) take a draw from the beta distribution described above for the fraction of B117/SGTF, (ii) impute B117 cases by multiplying SGTF cases by this draw from the beta, (iii) estimate the logistic growth model using this set of imputed B117 case numbers, and finally (iv) project future B117 prevalence using this model.  We combine all draws for projected B117 prevalence to integrate over uncertainty in the fraction of B117 to SGTF samples


mypredList <- vector("list", nsim)
mypredList_insample <- vector("list", nsim)

mypredList2 <- vector("list", nsim)
mypredList_insample2 <- vector("list", nsim)


## Matrix for imputed b117 values
bmat <- matrix(nrow = nrow(ut_fit), ncol = nsim)
bmat2 <- matrix(nrow = nrow(ut_fit), ncol = nsim)

## Store dataframes of imputed b117 values
ut_bk_list <- vector("list", nsim)
ut_bk2_list <- vector("list", nsim)

## For fitted stan_glm objects
mystan_list <- vector("list", nsim)
mystan_list2 <- vector("list", nsim)

mystan_samples_list <- vector("list", nsim)
mystan_samples_list2 <- vector("list", nsim)

## Date of first report (for fitting "original" estimates)
date_first_report <- ymd(initial_estimate_date)
for (k in 1:nsim) {
  
  ## Set the seed for reproducibility
  set.seed(k)
  
  cat(sprintf("Performing simulation %i out of %i...\n", k, nsim)) #
  
  ## Draw fraction of B117/SGTF
  mybeta <- rbeta(nrow(ut_fit), a1, b1)
  mybeta2 <- rbeta(nrow(ut_fit), a2, b2)
  
  ## Impute number of B117 cases
  bvec <- (ut_fit$n_Sdrop * mybeta)
  bmat[, k] <- bvec
  
  bvec2 <- (ut_fit$n_Sdrop * mybeta2) + ut_fit$n_Spos * runif(1, 0.01, 0.05)
  # second term accounts for the 1-5% chance of a non-SGTF being B117
  bmat2[, k] <- bvec2
  
  
  ## Create a dataframe for fitting the growth model and store it
  ut_bk <- ut_fit
  ut_bk['B117']<-round(bvec)
  ut_bk<-ut_bk%>%mutate(prop2 = (B117) / n_pos,
                        prop2_se = sqrt(prop2 * (1 - prop2) / n_pos))
  
  ut_bk_list[[k]] <- ut_bk
  
  ut_bk2 <- ut_fit
  ut_bk2['B117']<-round(bvec2)
  ut_bk2<-ut_bk2%>%mutate(prop2 = (B117) / n_pos,
                          prop2_se = sqrt(prop2 * (1 - prop2) / n_pos))
  
  ut_bk2_list[[k]] <- ut_bk2
  
  ## Build the model (original report date)
  myglm_stan2k <- rstanarm::stan_glm(
    cbind(B117, n_Sdrop + n_Spos - B117) ~ time,
    family=binomial,
    data = ut_bk %>%
      filter(n_pos >= 1) %>%
      filter(test_date <= date_first_report)
  )
  
  mystan_list[[k]] <- myglm_stan2k
  
  mystan_samples_list[[k]] <- as.matrix(myglm_stan2k)
  
  ## Predict the prevalence
  mypredList[[k]] <- posterior_epred(myglm_stan2k, newdata = ut_future)
  mypredList_insample[[k]] <- posterior_epred(myglm_stan2k, newdata = ut_bk)
  
  ## Build the model
  myglm_stan2k2 <- rstanarm::stan_glm(
    cbind(B117, n_Sdrop + n_Spos - B117) ~ time,
    family=binomial,
    data = ut_bk2 %>%
      filter(n_pos >= 1)                              
  )
  
  mystan_list2[[k]] <- myglm_stan2k2
  
  mystan_samples_list2[[k]] <- as.matrix(myglm_stan2k2)
  
  ## Predict the prevalence
  mypredList2[[k]] <- posterior_epred(myglm_stan2k2, newdata = ut_future)
  mypredList_insample2[[k]] <- posterior_epred(myglm_stan2k2, newdata = ut_bk)
  
  ## Print out an update
  cat(sprintf("\n\n***********************\n***********************\n\nSIM %i out of %i...\n\n***********************\n***********************\n\n", k, nsim))
  
}




logistic <- function(x) 1 / (1 + exp(-x))

## Concatenate samples of terms in model
mystan_samples <- mystan_samples_list %>% plyr::rbind.fill.matrix()
mystan_samples2 <- mystan_samples_list2 %>% plyr::rbind.fill.matrix()

colnames(mystan_samples)

###############################################################################
#           In sample fit...          #
###############################################################################

## Original posterior for growth rate (first report)
mystan_samples <- mystan_samples_list %>% plyr::rbind.fill.matrix()
mystan_samples[, 2] %>% mean() %>% round(3)
mystan_samples[, 2] %>% quantile(c(0.025, 0.975)) %>% round(3)

# intercept and growth rate (first report)
c1_distrib<-mystan_samples[,1]
k1_distrib<-mystan_samples[,2]

## Updated posterior for growth rate
mystan_samples2 <- mystan_samples_list2 %>% plyr::rbind.fill.matrix()
mystan_samples2[, 2] %>% mean() %>% round(3)
mystan_samples2[, 2] %>% quantile(c(0.025, 0.975)) %>% round(3)

# intercept and growth rate (second report)
c2_distrib<-mystan_samples2[,1]
k2_distrib<-mystan_samples2[,2]

## Matrix of prevalence estimates
mypredSim <- plyr::rbind.fill.matrix(mypredList)
mypredSim2 <- plyr::rbind.fill.matrix(mypredList2)

## Quantiles of prevalence estimates
mypred2_intsim <- apply(mypredSim, 2, quantile, probs = c(0.025, 0.5, 0.975, 0.25, 0.75))
mypred2_intsim2 <- apply(mypredSim2, 2, quantile, probs = c(0.025, 0.5, 0.975, 0.25, 0.75))

# Find the column that corresponds to the current date

iobs<-which.max(ut$test_date)
p_curr_distrib2<-mypredSim2[,iobs]



## Dataframe for plotting
plotdf <- ut_future 
ut_future['pred'] <- mypred2_intsim[2, ]
ut_future['pred_lo'] <- mypred2_intsim[1, ]
ut_future['pred_hi'] <- mypred2_intsim[3, ]
ut_future['pred_iqr_lo'] <- mypred2_intsim[4, ]
ut_future['pred_iqr_hi'] <- mypred2_intsim[5, ]
ut_future['pred2'] <- mypred2_intsim2[2, ]
ut_future['pred2_lo'] = mypred2_intsim2[1, ]
ut_future['pred2_hi'] = mypred2_intsim2[3, ]
ut_future['pred2_iqr_lo'] = mypred2_intsim2[4, ]
ut_future['pred2_iqr_hi'] = mypred2_intsim2[5, ]
## Deviance (goodness of fit)
ut_future['dev']<- -2 * (dbinom(round(rowMeans(bmat)), ut_future$n_pos, ut_future$pred, log=TRUE))
ut_future['dev2']<- -2 * (dbinom(round(rowMeans(bmat)), ut_future$n_pos, ut_future$pred2, log=TRUE))
## Remove in-sample proportions for dates not in dataset
ut_future<-ut_future%>%mutate(prop = ifelse(test_date > max(ut_fit$test_date), NA, prop),
                              prop_se = ifelse(test_date > max(ut_fit$test_date), NA, prop_se))
ut_future['blo']<- c(apply(bmat, 1, quantile, 0.025), rep(NA, n_future))
ut_future['bhi']<- c(apply(bmat, 1, quantile, 0.975), rep(NA, n_future))

ut_future<-ut_future%>%mutate(dataset = ifelse(test_date <= date_first_report,
                                               "Feb 12",
                                               "Apr 9")) %>%
  mutate(dataset = factor(dataset, levels = c("Feb 12", "Apr 9")))
posterior_list<-list(c1_distrib, k1_distrib, c2_distrib, k2_distrib, p_curr_distrib2, ut_future)
return(posterior_list)
}





get_daily_p_local<-function(test_data, ARRIVAL_DATE){
# Convert local infections vector into incidence vector
days_from_arrival<-as.numeric((test_data[length(test_data[,1]), "date"]) - as.Date(ARRIVAL_DATE))
days_vec<-c(0:days_from_arrival)
N_SAMPLES <-1000
time_to_symptoms <- rlnorm(N_SAMPLES, meanlog = 1.54, sdlog = 0.47)
time_to_test <- rexp(N_SAMPLES, rate = 1/2)
days_after_infection <- round(time_to_symptoms + time_to_test)

# Store whether sample is local or not
samples_local<-rep(0, nrow(test_data))


## Run samples ##
for(i in 1:N_SAMPLES){
  infection_times<-days_vec - sample(days_after_infection, size = length(days_vec), replace = TRUE)
  
  # Any negative numbers are importations, positive local
  for (m in 1:nrow(test_data)) {
    if (infection_times[m] >= 0) {
      samples_local[m] <- samples_local[m] + 1
    }
  }
}

## Calculate sample proportions for each day 
daily_p_local<-samples_local/N_SAMPLES
return(daily_p_local)
}

Rt_fxn_cases <- function(cases, case_data, arrival, last_reliable_day, MEAN_GI,STD_GI){
  daily_p_local<-get_daily_p_local(case_data, arrival)
  #days <- days_vec
  Q2p5 <- c()
  median <- c()
  Q97p5 <- c()
  BINOM_SAMPLES = 100
  RT_SAMPLES = 100
  local_sample_store <- data.frame()
  import_sample_store <- data.frame()
  ## Binomial samples for each day ##
  for (i in 1:length(daily_p_local)) {
    binom_sample <- rbinom(BINOM_SAMPLES, size = cases[i], prob = daily_p_local[i])
    import_cases <- rep(cases[i], BINOM_SAMPLES) - binom_sample
    
    Q2p5 <- append(Q2p5, quantile(binom_sample, 0.025, na.rm = TRUE))
    median <- append(median, quantile(binom_sample, 0.5, na.rm = TRUE))
    Q97p5 <- append(Q97p5, quantile(binom_sample, 0.975, na.rm = TRUE))
    
    # Store samples for Rt calculation
    colname <- toString(i-1)
    if (i == 1) {
      local_sample_store <- data.frame(binom_sample)
      colnames(local_sample_store) <- colname
      import_sample_store <- data.frame(import_cases)
      colnames(import_sample_store) <- colname
    } else {
      local_sample_store[[colname]] <- binom_sample
      import_sample_store[[colname]] <- import_cases
    }
  }
  
  # Dataframe for local vs. imported detected each day
  df <- data.frame(Q2p5, median, Q97p5, daily_p_local)
  df$daily_p_import <- 1 - daily_p_local
  
 
  SAMPLES_FROM_RT_DIST <- 1000
  
  for (i in 1:RT_SAMPLES) {
    sample_row_index <- round(runif(1, min = 1, max = BINOM_SAMPLES))
    local <- as.numeric(local_sample_store[sample_row_index,])
    import <- as.numeric(import_sample_store[sample_row_index,])
    incid <- data.frame(local, import)
    colnames(incid) <- c("local", "imported")
    
    # Calculate Rt
    config <- make_config(list(mean_si = MEAN_GI, std_si = STD_GI))
    Rt <- estimate_R(incid, method = "parametric_si", config = config)

    
    # Rename columns
    colnames(Rt$R)[11] <- "upperbound"
    colnames(Rt$R)[5]  <- "lowerbound"
    colnames(Rt$R)[8]  <- "median"
    colnames(Rt$R)[3]  <- "mean"
    colnames(Rt$R)[4]  <- "sd"
    
    # Make dataframe: each col is sample, rows are days
    if (i == 1) {
      Rt_means <- data.frame(Rt$R$t_end)
      Rt_stds <- data.frame(Rt$R$t_end)
      colnames(Rt_means) <- "time"
      colnames(Rt_stds)   <- "time"
    }
    colname <- toString(i)
    Rt_means[[colname]] <- Rt$R$mean
    Rt_stds[[colname]]   <- Rt$R$sd
    
  }
  
  Rt_times <- Rt_means$time
  Rt_medians <- rep(0, length(Rt_times))
  Rt_upperbounds  <- rep(0, length(Rt_times))
  Rt_lowerbounds <- rep(0, length(Rt_times))
  Rt_mean <- rep(0, length(Rt_times))
  Rt_var <- rep(0, length(Rt_times))
  Rt_a <- rep(0, length(Rt_times))
  Rt_b <- rep(0, length(Rt_times))
  
  # Summarize Rt estimates
  for (i in 1:nrow(Rt_means)) {
    Rt_aggregated_samples <- c()
    # For each Rt sample, use mean and std to 
    #   sample gamma, summarize all gamma draws
    #   for each day
    
    for (j in 2:RT_SAMPLES+1) {
      mean <- Rt_means[i, j]
      std  <- Rt_stds[i, j]
      var  <- std ** 2
      b <- mean / var
      a <- mean * b
      Rt_aggregated_samples <- append(
        Rt_aggregated_samples,
        rgamma(SAMPLES_FROM_RT_DIST, shape = a, rate = b)
      )
    }
    Rt_medians[i] <- quantile(Rt_aggregated_samples, .5)
    Rt_upperbounds[i] <- quantile(Rt_aggregated_samples, .975)
    Rt_lowerbounds[i] <- quantile(Rt_aggregated_samples, .025)
    Rt_mean[i]<-mean
    Rt_a[i]<-a
    Rt_b[i]<-b
    Rt_var[i]<-var
  }
  
  tvec = Rt_times - Rt_times[1]

  
  dates<-seq(as.Date(arrival), as.Date(last_reliable_day), "days")
  
  Rt_summary <- data.frame(dates)
  Rt_summary['Rt_times']<-Rt_times[1:length(dates)]
  Rt_summary['Rt_medians']<-Rt_medians[1:length(dates)]
  Rt_summary['Rt_lowerbounds']<-Rt_lowerbounds[1:length(dates)]
  Rt_summary['Rt_upperbounds']<-Rt_upperbounds[1:length(dates)]
  Rt_summary['Rt_mean']<-Rt_mean[1:length(dates)]
  Rt_summary['Rt_var']<-Rt_var[1:length(dates)]
  Rt_summary['Rt_a']<-Rt_a[1:length(dates)]
  Rt_summary['Rt_b']<-Rt_b[1:length(dates)]
  Rt_summary<-Rt_summary%>%filter(dates<=last_reliable_day-7)
  return(Rt_summary)
}



run_fall_SEIR<-function(par_table, I0bounds, Rt_fall){
# Get parameter names
  var_names<-colnames(par_table)
  for (i in 1:length(var_names)){
    assign(var_names[i], as.double(par_table[i]))
  }
  dates<-Rt_fall$dates
  tsim<-seq(from = 0, to = length(dates)-1, by=dt)
  Rt_cap<-2.5 # assume Rt never gets above this 
  # sample from uncertainty in Rts and
  Rts<-matrix(nrow=length(dates), ncol=nsamps) # rows are time, columns are samples 
  for (i in 1:length(dates)){
    mean_val<-Rt_fall$Rt_mean[i]
    if(mean_val>Rt_cap){ # Set a cap on R(t)
      mean_val <-Rt_cap
    }
    var<-Rt_fall$Rt_var[i]
    b<-mean_val/var
    a<-mean_val*b
    Rts[i,]<-rgamma(nsamps, shape = a, rate = b)
  }
  # Expand Rts
  Rts<-data.frame(Rts)
  Rts<-Rts%>%slice(rep(1:n(), each = 1/dt))
  
  
  
  
  # And from the initial number infected
  I0s<-runif(n=nsamps, min = I0bounds[1], max = I0bounds[2]) # sample uniformly from the range of I0s
  
  #Run 100 simulations randomly sampling from the Rt distribution at each time point
  # store outcomes in matrices
  cumImat<-c()
  Imat<-c()
  newImat<-c()
  betatmat<-c()
  Rtmat<-c()
  for (k in 1:nsamps){
    Rt<-Rts[,k]
    I0<-I0s[k] 
    R0 <- prop_prev_inf_Aug*N_POPULATION
    S<-rep(0,length(tsim))
    E<-rep(0,length(tsim))
    I<-rep(0,length(tsim))
    R<-rep(0,length(tsim))
    beta_t<-rep(0, length(tsim))
    cumI<-rep(0,length(tsim))
    newI<-rep(0, length(tsim))
    
    for (i in 1:length(tsim)){
      if (i==1){
        S[i]<-N_POPULATION-I0-R0 - delta/gamma*I0
        E[i]<-delta/gamma*I0
        I[i]<-I0
        R[i]<-R0
        # calculate beta from observed R(t)
        beta_t[i]<-Rt[i]*delta*N_POPULATION/S[i]
        cumI[i]<-R0 + I0
        newI[i]<-0
      }else if (i>1) {
        S[i]<-S[i-1] + dt*(-beta_t[i-1]*((S[i-1]*I[i-1])/N_POPULATION))
        E[i]<-E[i-1] + dt*(beta_t[i-1]*((S[i-1]*I[i-1])/N_POPULATION) - gamma*E[i-1])
        I[i]<-I[i-1] + dt*(gamma*E[i-1] - delta*I[i-1])
        R[i]<-R[i-1] + dt*(delta*I[i-1])
        # update beta from observed Rt
        beta_t[i]<-Rt[i]*delta*N_POPULATION/(S[i])
        cumI[i]<-cumI[i-1] + dt*(beta_t[i-1]*S[i-1]*I[i-1]/N_POPULATION)
        newI[i]<-dt*(beta_t[i-1]*S[i-1]*I[i-1]/N_POPULATION)
      }
    }
    newImat<-cbind(newImat, newI)
    cumImat<-cbind(cumImat, cumI)
    Imat<-cbind(Imat, I)
    betatmat<-cbind(betatmat, beta_t)
    Rtmat<-cbind(Rtmat, Rt)
    
    if(k==1){
      newIall<-data.frame(tsim)
      Iall<-data.frame(tsim)
    }
    colname<-toString(k)
    newIall[[colname]]<-newI
    Iall[[colname]]<-I
  }
  # Summary statistics from each of the matrices
  
  cumImedian<-apply( cumImat , 1 , quantile , probs = 0.5 , na.rm = TRUE )
  cumIlb<-apply( cumImat , 1 , quantile , probs = 0.025 , na.rm = TRUE )
  cumIub<-apply( cumImat , 1 , quantile , probs = 0.975 , na.rm = TRUE )
  Imedian<-apply( Imat , 1 , quantile , probs = 0.5 , na.rm = TRUE )
  Ilb<-apply( Imat , 1 , quantile , probs = 0.025 , na.rm = TRUE )
  Iub<-apply( Imat , 1 , quantile , probs = 0.975 , na.rm = TRUE )
  Rtmedian<-apply( Rtmat , 1 , quantile , probs = 0.5 , na.rm = TRUE )
  Rtlb<-apply( Rtmat , 1 , quantile , probs = 0.025 , na.rm = TRUE )
  Rtub<-apply( Rtmat , 1 , quantile , probs = 0.975 , na.rm = TRUE )
  betamedian<-apply( betatmat , 1 , quantile , probs = 0.5 , na.rm = TRUE )
  betalb<-apply( betatmat , 1 , quantile , probs = 0.025 , na.rm = TRUE )
  betaub<-apply( betatmat , 1 , quantile , probs = 0.975 , na.rm = TRUE )
  
  summary_df<-data.frame(tsim, Imedian, Ilb, Iub, cumImedian, cumIlb, cumIub)
  ikeep<-tsim == round(tsim)
  summary_df<-summary_df[ikeep,]
  summary_df$dates<-dates
  #return(summary_df)
  cumIJan<-as.numeric(cumImat[nrow(cumImat),])
  return(cumIJan)
  
}

run_spring_SEIR<-function(par_table, case_data_spring, Rt_spring, iRts, spring_last_date, R0s,k_tr_adv_distrib, k_p_distrib,
                          c_distrib, p_curr_distrib){
  # k_p_distrib is the growth rate use to describe the prevalence over time
  # k_tr_adv_distrib is the growth rate used to calculate transmissibility advantage
  # In one case (only using UT data) these will be the same, in the other (using literature values for k_tr_adv) these will differ

  var_names<-colnames(par_table)
  for (i in 1:length(var_names)){
    assign(var_names[i], as.double(par_table[i]))
  }
  
  # Assign dates
  springdates <- seq(Rt_spring$dates[1], spring_end_semester_date, "days")
  calibdates <- seq(Rt_spring$dates[1], spring_last_date, "days")
  projdates<-seq(spring_last_date, spring_end_semester_date, "days")
  
  # get time vectors
  tsim<-seq(from = 0, to = length(springdates)-1, by = dt)
  tobs<-seq(from = 0, to = length(calibdates)-1, by = dt)
  
  
  # Use PCT data from the first 4 days (some of which was entry testing) to get initial infection estimates
  n_pos_first_week <-sum(case_data_spring$student_PCT_pos[1:4])
  n_tests_first_week<-sum(case_data_spring$total_student_PCT_tests[1:4])
  # get the CI on the bounds of the initial proportion of student infected. 
  init_prevCI<-rbeta(n=nsamps, shape1 = 1+n_pos_first_week, shape2 = 1+n_tests_first_week-n_pos_first_week)%>%quantile(probs = c(0.5, 0.025, 0.975))
  I0lb<-init_prevCI[2]
  I0med<-init_prevCI[1]
  I0ub<-init_prevCI[3]
  
  f_ts<-matrix(nrow = length(tsim), ncol = nsamps)
  p_ts<-matrix(nrow = length(tsim), ncol = nsamps)
  
  # Assign the parameter combinations that we will apply to both the variant and no variant simulations for matched statistics.
  # R0s & cumIJan:
  R0s<-sample(R0s, nsamps)
  init_prev<-rbeta(n=nsamps, shape1 = 1+n_pos_first_week, shape2 = 1+n_tests_first_week-n_pos_first_week)
  for( i in 1:nsamps){
    p_curr<-p_curr_distrib[i] 
    k_p<-k_p_distrib[i]
    k_tr_adv<-k_tr_adv_distrib[i]
    c<-c_distrib[i]
    tr_adv<-exp(mean_GI*k_tr_adv) -1 # M-1
    p_ts[,i] <-(1/(1+exp(-(c+k_p*tsim))))
    f_ts[,i]<-(1+tr_adv*p_ts[,i])/(1+(tr_adv*p_curr)) 
    
  }
  
  
  f_t_median<-apply( f_ts , 1 , quantile , probs = 0.5 , na.rm = TRUE )
  f_t_lb<-apply( f_ts , 1 , quantile , probs = 0.025 , na.rm = TRUE )
  f_t_ub<-apply( f_ts , 1 , quantile , probs = 0.975 , na.rm = TRUE )
  df_test<-data.frame(tsim, f_t_median, f_t_lb, f_t_ub)
  
  Rts<-matrix(nrow=length(calibdates), ncol=nsamps)
  for (i in 1:length(calibdates)){
    mean_val<-Rt_spring$Rt_mean[i]
    var<-Rt_spring$Rt_var[i]
    b<-mean_val/var
    a<-mean_val*b
    Rts[i,]<-rgamma(nsamps, shape = a, rate = b)
  }
  
  # expand Rts 
  Rts<-data.frame(Rts)
  Rts<-Rts%>%slice(rep(1:n(), each = 1/dt))
  
  
  scenarios<-c('Faster spread', 'Slower spread')
  iRts_dt<-iRts*(1/dt)-9
  RthiCI<-quantile(Rts[iRts_dt[1],], probs = c(0.5, 0.025, 0.975))
  RtloCI<-quantile(Rts[iRts_dt[2],], probs = c(0.5, 0.025, 0.975))
  Rtlolb<-RtloCI[2]
  Rtlomed<-RtloCI[1]
  Rtloub<-RtloCI[3]
  Rthilb<-RthiCI[2]
  Rthimed<-RthiCI[1]
  Rthiub<-RthiCI[3]

  
  for (m in 1:length(scenarios)){
    # Set scenario
    iRt<-iRts_dt[m] 
    dates<-rep(springdates, nsamps)
    t<-rep(tsim, nsamps)
    scenario<-rep(scenarios[m], nsamps*length(tsim))
    for (j in 1:2){ # Run two times through for variant and no variant scenarions
      Imat<-matrix(nrow = length(tsim), ncol = nsamps)
      B117mat<-matrix(nrow = length(tsim), ncol = nsamps)
      cumIvec<-matrix(nrow = 1, ncol = nsamps) # total number infected in this time period
      cumIAVvec<-matrix(nrow = 1, ncol = nsamps) 
      Rvec<-matrix(nrow = 1, ncol = nsamps) # total immune at end of time period
      betavec<-matrix(nrow =1, ncol = nsamps) # Stores the constant beta for each simulation 
      cumIBVvec<-matrix(nrow = 1, ncol = nsamps)
      
      for (k in 1:nsamps){
        # get the parameter combo and initial condition
        Rt<-Rts[,k] # Get an Rt trajectory
        R0<-R0s[k]
        I0 = init_prev[k]*N_POPULATION
        variant_factor<-f_ts[,k]
        p_t<-p_ts[,k]
        
        
        # Preallocate compartments 
        S<-rep(0,length(tsim))
        E<-rep(0,length(tsim))
        I<-rep(0,length(tsim))
        R<-rep(0,length(tsim))
        cumI<-rep(0,length(tsim))
        cumIAV<-rep(0, length(tsim))
        newI<-rep(0, length(tsim))
        pV<-rep(0, length(tsim))
        beta_t<-rep(0, length(tsim))
        
        for (i in 1:length(tsim)){
          if (i==1){
            # Assign initial conditions without variant
            S[i]<-N_POPULATION-I0-R0 -delta/gamma*I0
            E[i]<-delta/gamma*I0
            I[i]<-I0
            R[i]<-R0
            cumI[i]<-I0
            beta_t[i]<-Rt[i]*delta*N_POPULATION/S[i] # calculate the first beta
            
          }
          # Before variant, infections in wildtype proceed using R(t) from cases
          else if (i>1 & i<=length(tobs)){
            S[i]<-S[i-1] + dt*(-(beta_t[i-1]*I[i-1])*(S[i-1]/N_POPULATION))
            beta_t[i]<-Rt[i]*delta*N_POPULATION/S[i]
            E[i]<-E[i-1] + dt*(beta_t[i-1]*((S[i-1]*I[i-1])/N_POPULATION) - gamma*E[i-1])
            I[i]<-I[i-1] + dt*(gamma*E[i-1] - delta*I[i-1])
            R[i]<-R[i-1] + dt*(delta*I[i-1])
            cumI[i]<-cumI[i-1] + dt*((beta_t[i-1]*I[i-1])*(S[i-1]/N_POPULATION))
            
            
          }
          # After tobs, introduce variant
          else if (i>length(tobs)) {
            
            if(j==2){
              variant_factor<- rep(1, length(tsim)) 
            }
            
            Rt_i<-Rt[iRt]
            # Hold beta constant except for the variant factor 
            beta_t[i]<-(Rt_i*delta*N_POPULATION/S[(length(tobs))])*variant_factor[i]
            S[i]<-S[i-1] + dt*(-(beta_t[i-1]*I[i-1])*(S[i-1]/N_POPULATION))
            E[i]<-E[i-1] + dt*(beta_t[i-1]*((S[i-1]*I[i-1])/N_POPULATION) - gamma*E[i-1])
            I[i]<-I[i-1] + dt*(gamma*E[i-1] - delta*I[i-1])
            R[i]<-R[i-1] + dt*(delta*I[i-1])
            cumI[i]<-cumI[i-1] + dt*((beta_t[i-1]*I[i-1])*(S[i-1]/N_POPULATION))
            cumIAV[i]<-cumIAV[i-1] + dt*((beta_t[i-1]*I[i-1])*(S[i-1]/N_POPULATION))
            
            
          }
          
          
        } # end time loop
        
        # Save total infections over time
        Imat[,k]<- I
        B117mat[,k]<-I*p_t
        cumIvec[1,k]<-cumI[length(tsim)]
        cumIAVvec[1,k]<-cumIAV[length(tsim)]
        cumIBVvec[1,k]<-cumI[length(tobs)+1] + R0
        Rvec[1,k]<-R[length(tsim)]
        betavec[1, k]<-beta_t[length(tobs)+1]
        
        
        colname<-toString(k)
        if(k==1){
          Ialls<-data.frame(tsim)
          Ialls[[colname]]<-I
        }
        else{
          Ialls[[colname]]<-I
        }
      } # end k loop for nsamps simulations
      
      # Assign Iall data frame to either be no variant or variant
      if (j==1){ # updated variant present
        Iallv<-Ialls
        # Get median and lower and upper bounds over time
        Imedianv<-apply( Imat , 1 , quantile , probs = 0.5 , na.rm = TRUE )
        Ilbv<-apply( Imat , 1 , quantile , probs = 0.025, na.rm = TRUE )
        Iubv<-apply( Imat , 1 , quantile , probs = 0.975, na.rm = TRUE )
        # Get the median and bounds of cumulative infections at end
        cumImedianv<-as.numeric(quantile(cumIvec, probs = 0.5, na.rm = TRUE))
        cumIlbv<-as.numeric(quantile(cumIvec, probs = 0.025, na.rm = TRUE))
        cumIubv<-as.numeric(quantile(cumIvec, probs = 0.975, na.rm = TRUE))
        cumIAVvmedian<-as.numeric(quantile(cumIAVvec, probs = 0.5, na.rm = TRUE))
        cumIAVvlb<-as.numeric(quantile(cumIAVvec, probs = 0.025, na.rm = TRUE))
        cumIAVvub<-as.numeric(quantile(cumIAVvec, probs = 0.975, na.rm = TRUE))
        cumIBVmedian<-as.numeric(quantile(cumIBVvec/N_POPULATION, probs = 0.5, na.rm = TRUE))
        cumIBVlb<-as.numeric(quantile(cumIBVvec/N_POPULATION, probs = 0.025, na.rm = TRUE))
        cumIBVub<-as.numeric(quantile(cumIBVvec/N_POPULATION, probs = 0.975, na.rm = TRUE))
        # Save the overall cumulative infections for the spring semester and the cumulative infections
        #after the variant, in order by the parameter combo
        cumIvec_v<-cumIvec
        cumIAVvec_v<-cumIAVvec
        
      }
      
      
      if (j==2){ # no variant
        Iallnv<-Ialls # assign I all to no variant
        # Get median and lower and upper bounds over time
        Imediannv<-apply( Imat , 1 , quantile , probs = 0.5 , na.rm = TRUE )
        Ilbnv<-apply( Imat , 1 , quantile , probs = 0.025, na.rm = TRUE )
        Iubnv<-apply( Imat , 1 , quantile , probs = 0.975, na.rm = TRUE )
        B117medianv<-apply( B117mat , 1 , quantile , probs = 0.5 , na.rm = TRUE )
        B117lbv<-apply( B117mat , 1 , quantile , probs = 0.025, na.rm = TRUE )
        B117ubv<-apply( B117mat , 1 , quantile , probs = 0.975, na.rm = TRUE )
        WTmedianv<-Imedianv-B117medianv
        WTlbv<-Ilbv-B117lbv
        WTubv<-Iubv-B117ubv
        # Get the median and bounds of cumulative infections at end
        cumImediannv<-as.numeric(quantile(cumIvec, probs = 0.5, na.rm = TRUE))
        cumIlbnv<-as.numeric(quantile(cumIvec, probs = 0.025, na.rm = TRUE))
        cumIubnv<-as.numeric(quantile(cumIvec, probs = 0.975, na.rm = TRUE))
        # Get the median and bounds of total immune at end
        cumIAVnvmedian<-as.numeric(quantile(cumIAVvec, probs = 0.5, na.rm = TRUE))
        cumIAVnvlb<-as.numeric(quantile(cumIAVvec, probs = 0.025, na.rm = TRUE))
        cumIAVnvub<-as.numeric(quantile(cumIAVvec, probs = 0.975, na.rm = TRUE))
        
        # Save the overall cumulative infections for the spring semester and the cumulative
        # infections after the variant, in order by the parameter combo
        cumIvec_nv<-cumIvec
        cumIAVvec_nv<-cumIAVvec
      }
      
    } # end j loop for variant no variant
    
    # Compute the summary statistics from the matched samples
    # April 9th to May 23rd (after variant)
    delta_inf_AV<- cumIAVvec_v - cumIAVvec_nv
    delta_inf_AV_med<-quantile(delta_inf_AV, probs = 0.5, na.rm = TRUE)
    delta_inf_AV_lb<-quantile(delta_inf_AV, probs = 0.025, na.rm = TRUE)
    delta_inf_AV_ub<-quantile(delta_inf_AV, probs = 0.975, na.rm = TRUE)
    pct_inc_AV<-100*(cumIAVvec_v - cumIAVvec_nv)/cumIAVvec_nv
    pct_inc_AV_med<-quantile(pct_inc_AV, probs = 0.5, na.rm = TRUE)
    pct_inc_AV_lb<-quantile(pct_inc_AV, probs = 0.025, na.rm = TRUE)
    pct_inc_AV_ub<-quantile(pct_inc_AV, probs = 0.975, na.rm = TRUE)
    # Overall (whole time period)
    delta_inf<-cumIvec_v- cumIvec_nv
    delta_inf_med<-quantile(delta_inf, probs = 0.5, na.rm = TRUE)
    delta_inf_lb<-quantile(delta_inf, probs = 0.025, na.rm = TRUE)
    delta_inf_ub<-quantile(delta_inf, probs = 0.975, na.rm = TRUE)
    pct_inc<-100*(cumIvec_v - cumIvec_nv)/cumIvec_nv
    pct_inc_med<-quantile(pct_inc, probs = 0.5, na.rm = TRUE)
    pct_inc_lb<-quantile(pct_inc, probs = 0.025, na.rm = TRUE)
    pct_inc_ub<-quantile(pct_inc, probs = 0.975, na.rm = TRUE)
    
    # Make the data long to stack all of the simulation runs
    I_longvariant<-Iallv%>%melt(id.vars = c("tsim"), variable.name = "samples", value.name = "I_t_v")
    I_longnv<-Iallnv%>%melt(id.vars = c("tsim"), variable.name = "samples", value.name = "I_t_nv")
    I_variant<-I_longvariant$I_t_v
    I_novariant<-I_longnv$I_t_nv
    samples<-I_longvariant$samples
    if (m==1){
      # First make a long data frame with all of the samples stacked (nrows = n_timesteps*nsamples)
      df_t_sim<-data.frame(t, scenario, samples, I_variant,  I_novariant) # Long one with all sims
      scenario<-scenario[1:length(tsim)]
      ik<-t==round(t) # only keep who timesteps (days)
      df_t_sim<-df_t_sim[ik,]
      df_t_sim$dates<-rep(springdates, nsamps)
      
      # Then make the time course summary data frame (nrows = n_timesteps)
      df_t<-data.frame(tsim, scenario, Imedianv, Ilbv, Iubv, 
                       Imediannv, B117medianv, B117lbv, B117ubv,
                       WTmedianv, WTlbv, WTubv,Ilbnv, Iubnv, f_t_median, f_t_lb,
                       f_t_ub) # time course and scenario
      scenario<-scenario[1]
      ikeep<-tsim == round(tsim) # only keep whole timesteps (days)
      df_t<-df_t[ikeep,]
      df_t$springdates<-springdates
      
      df_summary<-data.frame(scenario,delta_inf_AV_med, delta_inf_AV_lb, delta_inf_AV_ub,
                             pct_inc_AV_med, pct_inc_AV_lb, pct_inc_AV_ub,
                             delta_inf_med, delta_inf_lb, delta_inf_ub,
                             pct_inc_med, pct_inc_lb, pct_inc_ub, I0med, I0lb, I0ub,
                             cumImedianv, cumIlbv, cumIubv,
                             cumImediannv, cumIlbnv, cumIubnv,
                             cumIAVvmedian, cumIAVvlb, cumIAVvub, cumIAVnvmedian, cumIAVnvlb, cumIAVnvub,
                             cumIBVmedian, cumIBVlb, cumIBVub)
      
      
      
    }
    else{
      # all simulations
      df_t_simi<-data.frame(t,scenario, samples, I_variant, I_novariant)
      ik<-t==round(t)
      df_t_simi<-df_t_simi[ik,]
      df_t_simi$dates<-rep(springdates, nsamps)
      df_t_sim<-rbind(df_t_sim, df_t_simi) # stack the dataframes
      # just median and lb/ubs
      scenario<-scenario[1:length(tsim)]
      df_ti<-data.frame(tsim, scenario, Imedianv, Ilbv, Iubv,
                        Imediannv,B117medianv, B117lbv, B117ubv,
                        WTmedianv, WTlbv, WTubv, Ilbnv, Iubnv, f_t_median, f_t_lb,
                        f_t_ub)
      ikeep<-tsim == round(tsim)
      df_ti<-df_ti[ikeep,]
      df_ti$springdates<-springdates
      df_t<-rbind(df_t, df_ti)
      # summary statistics
      scenario<-scenario[1]
      df_summaryi<-data.frame(scenario,delta_inf_AV_med, delta_inf_AV_lb, delta_inf_AV_ub,
                              pct_inc_AV_med, pct_inc_AV_lb, pct_inc_AV_ub,
                              delta_inf_med, delta_inf_lb, delta_inf_ub,
                              pct_inc_med, pct_inc_lb, pct_inc_ub, I0med, I0lb, I0ub,
                              cumImedianv, cumIlbv, cumIubv,
                              cumImediannv, cumIlbnv, cumIubnv,
                              cumIAVvmedian, cumIAVvlb, cumIAVvub, cumIAVnvmedian, cumIAVnvlb, cumIAVnvub,
                              cumIBVmedian, cumIBVlb, cumIBVub)
      df_summary<-rbind(df_summary, df_summaryi)
    }
    
    
    
  } # end m loop for scenarios
  out_list<-list(df_t_sim, df_t, df_summary)
  return(out_list)
}
  