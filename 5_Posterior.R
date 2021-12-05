library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(bayesplot)
library(matrixStats)
library(rethinking)
library(data.table)
library(reshape2)
library(ggrepel)
library(performance)
library(Metrics)

# Sampler diagnostics 
  sampler_params <- get_sampler_params(fit_GP_unp_int, inc_warmup = FALSE)
  sampler_params_chain1 <- sampler_params[[1]]
  colnames(sampler_params_chain1)
  
  #average value of acceptance probabilities of all possible samples by chain 
  mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
  print(mean_accept_stat_by_chain)
  
  # maximum value of treedepth for each chain 
  max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
  print(max_treedepth_by_chain)
  
# Model Diagnostics
  #Check Rhat 
  #Figure B3
  rhat_GP <- rhat(fit_GP_unp_int)
  mcmc_rhat_hist(rhat_GP)
  rhat <- as.data.frame(rhat_GP)
  rhat$names <- rownames(rhat)
  #check if and where Rhat is > 1.01
  rhat[( rhat_GP > 1.01),]

  #ESS: Check n_eff/N ratio (ratio of effective sample size (n_eff) to total sample size (N))
  #Figure B2
  ratios_GP <- neff_ratio(fit_GP_unp_int)
  mcmc_neff(ratios_GP)
  neff <- as.data.frame(ratios_GP)
  neff$names <- rownames(neff)
  #check if and where neff is < 0.1
  neff[( ratios_GP < 0.1),]

  fit_summary <- summary(fit_GP_unp_int)
  print(names(fit_summary))
  print(fit_summary$summary)

# Draw posterior samples 
post_evi <- as.matrix(fit_GP_unp_int, pars= "evi_pred")
#create data table with median predicted evi (from posterior distribution) for each observed evi
post_e <- data.table(evi_pred = colMedians(post_evi), sd = apply(post_evi, 2, sd), 
          lower = apply(post_evi, 2, quantile, 0.055), #ci_89
          upper = apply(post_evi, 2, quantile, 0.945), #ci_89
          evi_obs = standat$evi)
  #SAVE/LOAD
  # save(post_e, file = "post_e.RData")
  # load("~/01Master/MasterThesis/Pius/R/sand dam/post_e.RData")

# Posterior predicitve check (Figure B1): 
  #PPC distribution - comparison of empirical distribution (y) 
  #to simulated posterior predicitve distributions (yrep)

  ppc_dens_overlay(y = post_e$evi_obs, #observed
                   yrep = post_evi) #predicted
  
# Model performance: Crossvalidation (Figure 6)
  # plot observed (y) vs. predicted (x) EVI  (Pineiro et al., 2008)
  ggplot(post_e) +
    labs(x= "Predicted EVI", y="Observed EVI", title="Model performance") +
    geom_errorbarh(aes(y=evi_obs, xmin = lower, xmax=upper), color="darkgrey") + #CI = 89% 
    geom_point(aes(x=evi_pred, y=evi_obs),alpha = 1/5) + 
    geom_abline(intercept = 0, slope = 1, color="blue", size=0.7) +
    theme_bw() 
  #R-squared
  cor(post_e$evi_pred, post_e$evi_obs)^2
  #RMSE
  sqrt(mean((post_e$evi_obs - post_e$evi_pred)^2))

#SLOPES:
post_s <- as.matrix(fit_GP_unp_int, pars = c("b1", "b2", "b3[1]", "b3[2]", "b4")) 
  #SAVE/LOAD
    # save(post_s, file = "post_s.RData")
    # load("~/01Master/MasterThesis/Pius/R/sand dam/post_s.RData")
  
  #vizualise slopes:
  #Figure 7
  mcmc_areas(post_s, pars=c("b1", "b2", "b4"), 
             prob = .89, point_est = "median") + 
    scale_y_discrete(labels=c(
      "b1" = expression(paste("Precipitation (P) " ,beta[1])),
      "b2" =  expression(paste("Sand Dam Presence (SD) ",beta[2])),
      "b4" = expression(paste("Interaction P:SD ",beta[4]))
      ))
  # Figure 8            
  mcmc_areas(post_s, regex_pars = "b3", prob=0.89, point_est = "median") +
    scale_y_discrete(labels=c(
      "b3[1]" = expression(paste("Shrubs " ,beta["3(1)"])),
      "b3[2]" =  expression(paste("Cropland ",beta["3(2)"]))
      ))
        

load("~/01Master/MasterThesis/Pius/R/sand dam/df_seas_all.RData")
df_s <- df_s %>%
  filter(!is.na(presence))
post_dt <- post_e
post_dt$P <- df_s$P # unit (mm)
post_dt$presence <- standat$presence #0/1
post_dt$lc_class <-df_s$class #shrubs, cropland 
post_dt$time <- standat$time #timesteps [96]
post_dt$gp_id <- standat$gp_id
post_dt$season <- df_s$season
post_dt$area <- df_s$area_m2
post_dt$year <- df_s$year

###Visualize posterior samples 
#Overall: EVI and Precipitation
pres.labs <- c("0" = "Sand dams absent (0)",
               "1" = "Sand dams present (1)")

  #grouped by land cover class (Figure C1)
  post_dt %>%
    group_by(lc_class) %>%
    mutate(w_med = weightedMedian(evi_pred, area)) %>%
    ggplot(aes(x=P, y=evi_pred)) +
    geom_errorbar(aes(ymin=lower, ymax = upper), color = "grey") +
    geom_point(alpha = 0.4) +
    facet_grid(cols= vars(lc_class)) +
    labs(x="Precipitation (mm)", y="EVI estimates", title = "EVI and Precipitation (by Land Cover Class)") +
    theme_bw() +
    geom_hline(aes(yintercept = w_med, group = lc_class), color="blue", linetype = "dashed")

  #grouped by sand dam presence (Figure C2)
  post_dt %>%
    group_by(presence) %>%
    mutate(w_med = weightedMedian(evi_pred, area), 
           m_med = median(evi_pred)) %>% 
    ggplot(aes(x=P, y=evi_pred)) +
    geom_errorbar(aes(ymin=lower, ymax=upper), color="grey") + #, alpha=0.1
    geom_point(alpha = 0.4) +
    facet_grid(cols = vars(presence), labeller=labeller(presence = pres.labs)) +
    labs(x="Precipitation (mm)", y="EVI estimates", title="EVI and Precipitation (by Sand Dam Presence)") +
    theme_bw() +
    geom_hline(aes(yintercept=w_med, group = presence), color="blue", linetype="dashed") 


#SEASONALITY
  #Figure 9 
  #A
  post_dt %>%
    group_by(!!!syms(c("season", "presence"))) %>%
    mutate(m_med = median(evi_pred),
           w_med = weightedMedian(evi_pred, area)) %>%
    ggplot(aes(x=P, y=evi_pred, color=factor(lc_class))) +
    geom_point(alpha=0.25) +
    facet_grid(season ~ presence, 
               labeller=labeller(presence = pres.labs)) +# ,scales="free" #labeller=pres.labs
    geom_errorbar(aes(ymin=lower, ymax=upper, color=factor(lc_class)), alpha=0.15) +
    labs(x="Precipitation (mm)", y="EVI estimates", title="a") +  
    theme_bw() +
    scale_color_manual(name="Land cover",labels=c("Cropland", "Shrubs"), values=c("#D55E00", "#009E73")) +
    theme(legend.position = "none") 
  
  #B
  post_dt %>%
    group_by(!!!syms(c("season", "presence", "lc_class"))) %>%
    mutate(w_evi = weightedMedian(evi_pred, area),
           w_lower = weightedMedian(lower, area),
           w_upper = weightedMedian(upper, area)) %>%
    ggplot() +    
    facet_grid(cols=vars(lc_class), rows=vars(season)) +
    geom_errorbar(aes(x=factor(presence), ymin=w_lower, ymax=w_upper, color=factor(lc_class), width=0.3)) +
    scale_color_manual(name="",labels=c("cropland", "shrubs"), values=c("#D55E00", "#009E73")) +
    geom_point(aes(x=factor(presence), y=w_evi, color=factor(lc_class)), size=3) +
    theme_bw() +
    labs(x="Sand dams", y="EVI estimates", title="b") +
    scale_x_discrete(breaks = c(0,1), labels = c("absent","present")) 

  #C: Precipitation over seasons 
  post_dt %>%
    group_by(season) %>%
    mutate(w_P_med = weightedMedian(P, area),
           Pmed = median(P),
           lower_p = quantile(P, probs = 0.055),
           upper_p = quantile(P, probs = 0.945),
           lower_pw = weightedMedian(lower_p, area),
           upper_pw = weightedMedian(upper_p, area))%>%
    ggplot() +    
    facet_grid(rows=vars(season)) +
    geom_errorbar(aes(y=1, xmin=lower_pw, xmax=upper_pw)) + 
    geom_point(aes(x=w_P_med, y=1), size=3) +
    theme_bw() +
    labs(x="Precipitation (mm)", y="season", title="c") +
    scale_y_discrete(breaks = NULL) +
    ylab(NULL)

  
#OVER YEARS 
### Figure C3
   #A
   post_dt %>%
     group_by(!!!syms(c("year", "presence"))) %>%
     mutate(m_med = median(evi_pred),
            w_med = weightedMedian(evi_pred, area)) %>%
     ggplot(aes(x=P, y=evi_pred, color=factor(lc_class))) +
     geom_point(alpha=0.25) +
     facet_grid(year ~ presence, 
                labeller=labeller(presence = pres.labs)) +# ,scales="free" #labeller=pres.labs
     geom_errorbar(aes(ymin=lower, ymax=upper, color=factor(lc_class)), alpha=0.15) +
     labs(x="Precipitation [mm]", y="EVI estimates", title="a") +  
     theme_bw() +
     scale_color_manual(name="Land cover",labels=c("Cropland", "Shrubs"), values=c("#D55E00", "#009E73")) 
   
   #B
   post_dt %>%
    group_by(!!!syms(c("year", "presence", "lc_class"))) %>%
    mutate(w_evi = weightedMedian(evi_pred, area),
           w_lower = weightedMedian(lower, area),
           w_upper = weightedMedian(upper, area)) %>%
    ggplot() +    
    facet_grid(cols=vars(lc_class), rows=vars(year)) +
    geom_errorbar(aes(x=factor(presence), ymin=w_lower, ymax=w_upper, color=factor(lc_class), width=0.3)) +
    scale_color_manual(name="",labels=c("cropland", "shrubs"), values=c("#D55E00", "#009E73")) +
    geom_point(aes(x=factor(presence), y=w_evi, color=factor(lc_class)), size=2) +
    theme_bw() +
    labs(x="Sand dams", y="EVI estimates", title="b") +
    scale_x_discrete(breaks = c(0,1), labels = c("absent","present")) +
    theme(legend.position = "bottom") 

  #C: Precipitation histogram 
   b <-  post_dt %>%
     group_by(year) %>%
     mutate(m_P = median(P),
            w_P = weightedMedian(P, area)) %>%
     ggplot(aes(y=P)) +
     geom_histogram(alpha=0.2) + 
     facet_grid(cols=vars(year)) +
     labs(x="Count", y="Precipitation (mm)", title="b. Precipitation") + 
     theme_bw() +
     geom_hline(aes(yintercept=w_P, group=year), colour="black", linetype="dashed")
   

###Figure 10 
  e<-post_dt %>%
    group_by(!!!syms(c("year", "presence"))) %>%
    mutate(w_evi = weightedMedian(evi_pred, area),
           w_lower = weightedMedian(lower, area),
           w_upper = weightedMedian(upper, area)) %>%
    ggplot() +    
    facet_grid(cols=vars(year)) +
    geom_errorbar(aes(x=presence, ymin=w_lower, ymax=w_upper, color=factor(presence), width=0.5), size=1) +
    geom_point(aes(x=presence, y=w_evi, color=factor(presence)), size=3) +
    theme_bw() +
    labs(x=" ", y="EVI", title="a. EVI")  +
    scale_color_manual(name="Sand dam",labels=c("absent", "present"), values=c("#CC79A7", "#F0E442")) +
    scale_x_discrete(breaks = NULL) +
    xlab(NULL) +
    theme(legend.position="bottom")

grid.arrange(e,b)
