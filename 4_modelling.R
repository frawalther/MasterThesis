library(dplyr)
library(rstan)
library(lubridate)
library(zoo)
library(tidyr)
library(bayesplot)
library(loo)
library(shinystan)
library(purrr)
library(ggsci)
library(tictoc)
library(ggplot2)
library(gridExtra)

load("~/01Master/MasterThesis/Pius/R/sand dam/df_com_all.RData")

df_par <- df_com
# aggregate over the 4 periods
    df_y_m <- df_par %>%
      mutate (year = year(year_month)) %>%
      mutate (month= month(year_month, label=T))
    
    df_y_m$season <- factor(df_y_m$month)
    
    levels(df_y_m$season) <- list(JF = month.abb[c(1,2)],
                                       MAM = month.abb[c(3:5)],
                                       JJAS = month.abb[c(6:9)],
                                       OND = month.abb[c(10:12)])
                            #based on literature and expert knowledge
        
    df_y_m$seas_y <- paste(df_y_m$season, df_y_m$year, df_y_m$X, sep="_")
        
    df_season <- df_y_m %>%
      group_by(seas_y) %>%
      summarize(EVI_mean = mean(EVI_mean),
                Precip_mean = mean(Precip_mean),
                time_mean = mean(timeorder)) %>% 
      as.data.frame()
        
    df_s <- merge(df_season, df_y_m, by="seas_y")
        
    df_s <- df_s %>%
      mutate(evi = EVI_mean.x,
             P = Precip_mean.x) %>%
      select(-EVI_mean.y, -EVI_mean.x, -Precip_mean.y, -Precip_mean.x, -date.x, -date.y, -month, -timeorder, -X.1, -year_month) #, -year_month
    
    df_s <- distinct(df_s)
        
    #SAVE/LOAD
      # save(df_s, file = "df_seas_all.RData")
      # load("~/01Master/MasterThesis/Pius/R/sand dam/df_seas_all.RData")
    
    #in presence column are NAs 
    nrow(df_s[is.na(df_s$presence), ])  
    unique(df_s$seas_y)
    #removed:
    df_s <- df_s %>%
      drop_na("presence")
    
    # add/wrangle season column 
    df_s <- df_s %>%
      transform(seas=as.numeric(factor(season)))
    unique(df_s$seas)
    #1 = JF
    #3 = JJAS
    #2 = MAM
    #4 = OND
    
    df_s$gp_id <- paste(df_s$X, df_s$class, sep="_") #n=336
    #save(df_s, file = "df_s.RData")

### CREATE LIST OF DATA FOR STAN MODEL ###
    standat = with(df_s, list(
      evi = evi,
      P = P,
      time = time_mean,
      presence = presence, 
      lc_class = as.integer(factor(LC_proj)),
      gp_id = as.integer(factor(gp_id)),
      seas = as.integer(factor(seas))
    ))
    
    standat$gp_sampsize = table(standat$gp_id)
    standat$max_gp_sampsize = max(standat$gp_sampsize)
    standat$ngp = max(standat$gp_id)
    standat$N = length(standat$evi)
    standat$n_lc = max(standat$lc_class)
    standat$n_seas = max(standat$seas)

    #rescale precipitation data
        #standardizing by substracting the mean and dividing by 2 standard deviations
        #(Gelman 2021, p. 186)
    standat$P <- (standat$P - mean(standat$P))/(2*sd(standat$P))

#############################
### Gausian Process Model ###
#############################  

GP_unp_int = stan_model("~/01Master/MasterThesis/Pius/R/MasterThesis/MasterThesis/GP_model.stan")

#job::job({
tic()
  fit_GP_unp_int = sampling (GP_unp_int, data=standat, 
                             chains=4, 
                             cores=4,
                             iter=6000,
                             warmup=4000,
                             control= list(adapt_delta=0.92)) 
toc()
#}) 

check_hmc_diagnostics(fit_GP_unp_int)
print(fit_GP_unp_int, pars = fit_GP_unp_int@sim$pars_oi, probs = c(0.045, 0.955), digits = 2)

#   GP_int <- as.matrix(fit_GP_unp_int)
#   save(GP_int, file = "post_std.RData")
#   load("~/01Master/MasterThesis/Pius/R/sand dam/post_std.RData")
