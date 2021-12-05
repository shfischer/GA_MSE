### ------------------------------------------------------------------------ ###
### analyse "multi-species" GA runs ####
### ------------------------------------------------------------------------ ###

library(doParallel)
library(doRNG)
library(GA)
library(ggplot2)
library(scales)
library(cowplot)
library(Cairo)
library(tidyr)
library(dplyr)
library(FLCore)
library(FLash)
library(FLBRP)
library(mseDL)

source("funs_GA.R")
source("funs.R")

trans_from <- function(from = 1) {
  trans <- function(x) x - from
  inv <- function(x) x + from
  trans_new("from", trans, inv, 
            domain = c(from, Inf))
}


stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)

### ------------------------------------------------------------------------ ###
### uncertainty cap ####
### ------------------------------------------------------------------------ ###
### GA runs for pollack

ga_solution <- function(object) {
  res <- tail(object@bestSol, 1)[[1]][1, ]
  names(res) <- object@names
  res[c(1:4, 8)] <- round(res[c(1:4, 8)])
  res[c(5:7)] <- round(res[c(5:7)], 1)
  res[c(9:11)] <- round(res[c(9:11)], 2)
  return(res)
}

res <- foreach(fhist = c("one-way", "random"), .combine = bind_rows) %:%
  foreach(obj = c("PA", "MSY", "MSYPA"), .combine = bind_rows) %:%
  foreach(params = c("default", "multiplier", "cap", "cap_multiplier", "full", 
                     "full_cap"), .combine = bind_rows) %:%
  foreach(stat_yrs = c("all", "last10"), .combine = bind_rows) %do% {#browser()
    #browser()
    #print(paste(fhist, obj, params, stat_yrs))
    ### load data
    par_file <- switch(params,
      "default" = "multiplier-upper_constraint-lower_constraint", 
      "multiplier" = "multiplier", 
      "cap" = "upper_constraint-lower_constraint", 
      "cap_multiplier" = "multiplier-upper_constraint-lower_constraint", 
      "full" = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                      "interval-multiplier"), 
      "full_cap" = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                      "interval-multiplier-upper_constraint-lower_constraint")
    )
    obj_file <- switch(obj,
      "PA" = "obj_ICES_PA2",
      "MSY" = "obj_SSB_C_risk_ICV",
      "MSYPA" = "obj_ICES_MSYPA"
    )
    path <- paste0("output/500_50/uncertainty_cap/", fhist, "/pol/",
                   par_file, "--", obj_file)
    path_runs <- paste0(path, "_runs", 
                        ifelse(stat_yrs == "all", "", paste0("_", stat_yrs)),
                        ".rds")
    path_res <- paste0(path, "_res", 
                        ifelse(stat_yrs == "all", "", paste0("_", stat_yrs)),
                        ".rds")
    ### use GA paper results for "full" GA
    # if (isTRUE(obj == "MSY" & params == "full")) {
    #   path <- paste0("output/500_50/ms/trial/", fhist, "/pol/",
    #                par_file, "--", obj_file)
    # }
    if (!file.exists(path_runs)) return(NULL)
    print("found something")
    ga_res <- readRDS(path_res)
    ga_runs <- readRDS(path_runs)
    
    ### optimised parameters
    if (isFALSE(params == "default")) {
      pars <- ga_solution(ga_res)
    } else {
      pars <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0)
    }
    pars[which(is.nan(pars))] <- Inf
    tmp <- as.data.frame(t(pars))
    names(tmp) <- c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
                    "exp_r", "exp_f", "exp_b", "interval", "multiplier",
                    "upper_constraint", "lower_constraint")
    #if (is.nan(tmp$upper_constraint)) tmp$upper_constraint <- Inf
    tmp$obj <- obj
    tmp$fhist <- fhist
    tmp$stat_yrs_obj <- stat_yrs
    tmp$ga_obj <- params
    
    ### stats
    par_scn <- pars
    if (isTRUE(length(which(is.na(par_scn))) > 0)) {
      par_scn <- par_scn[-which(is.na(par_scn))]
    }
    par_scn <- paste0(par_scn, collapse = "_")
    stats_tmp <- ga_runs[[par_scn]]
    stats_tmp <- as.data.frame(lapply(as.data.frame(t(stats_tmp$stats)), unlist))
    
    ### combine pars and stats
    stats_tmp <- cbind(tmp, stats_tmp)
    
    ### if different stat_yrs period used, extract also default stats
    if (isFALSE(stat_yrs == "all")) {
      stats_tmp <- rbind(stats_tmp, stats_tmp)
      stats_tmp$stat_yrs <- c("all", stat_yrs)
      stats_names <- c("risk_Blim", "risk_Bmsy", "risk_halfBmsy", 
                       "risk_collapse", "SSB", "Fbar", "Catch", "SSB_rel", 
                       "Fbar_rel", "Catch_rel", "ICV")
      stats_tmp[2, stats_names] <- stats_tmp[2, paste0(stats_names, "_", stat_yrs)]
      stats_tmp[, paste0(stats_names, "_", stat_yrs)] <- NULL
    } else {
      stats_tmp$stat_yrs <- "all"
      ### remove redundant stats
      stats_tmp[, grep(x = names(stats_tmp), pattern = "_last10")] <- NULL
    }
    
    ### recreate fitness
    yr_suffix <- ""
    if (isTRUE(obj == "PA")) {
      stats_tmp$fitness <- sapply(seq(nrow(stats_tmp)), function(x) {
        sum(stats_tmp[x, "Catch_rel"]) -
          sum(penalty(x = stats_tmp[x, "risk_Blim"], negative = FALSE, max = 5,
                      inflection = 0.06, steepness = 0.5e+3))
      })
    } else if (isTRUE(obj == "MSY")) {
      stats_tmp$fitness <- sapply(seq(nrow(stats_tmp)), function(x) {
        -sum(abs(stats_tmp[x, "SSB_rel"] - 1),
                          abs(stats_tmp[x, "Catch_rel"] - 1),
                          stats_tmp[x, "ICV"], 
                          stats_tmp[x, "risk_Blim"])
      })
    } else if (isTRUE(obj == "MSYPA")) {
      stats_tmp$fitness <- sapply(seq(nrow(stats_tmp)), function(x) {
        -sum(abs(stats_tmp[x, "SSB_rel"] - 1),
                          abs(stats_tmp[x, "Catch_rel"] - 1),
                          stats_tmp[x, "ICV"], 
                          penalty(x = stats_tmp[x, "risk_Blim"], 
                                  negative = FALSE, max = 5, 
                                  inflection = 0.06, steepness = 0.5e+3))
      })
    }
    return(stats_tmp)
}

res %>%
  filter(obj == "PA")
saveRDS(res, file = "output/500_50/uncertainty_cap/results.rds")
res <- readRDS("output/500_50/uncertainty_cap/results.rds")


### format for plotting
stats_pol <- res %>%
  select(obj, fhist, stat_yrs_obj, stat_yrs, ga_obj, risk_Blim, SSB_rel, Fbar_rel, Catch_rel,
         ICV, fitness) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV, fitness), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel", "risk_Blim", 
                                       "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value"))) %>%
  mutate(scenario = factor(ga_obj,
                           levels = c("default", "multiplier", "cap", 
                                      "cap_multiplier", "full", "full_cap"),
                           labels = c("default", "GA multiplier", "GA cap", 
                                      "GA cap+\nmultiplier", 
                                      "GA all w/o cap", "GA all")))
stats_targets <- data.frame(stat = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                     "B[lim]~risk", "ICV", "fitness~value"),
                            target = c(1, 1, 1, 0, 0, NA))


plot_stats <- function(SSB_min = 0, SSB_max = NA,
                       F_min = 0, F_max = NA,
                       C_min = 0, C_max = NA,
                       risk_min = 0, risk_max = NA,
                       ICV_min = 0, ICV_max = NA,
                       fitness_min = NA, fitness_max = NA,
                       obj = "MSY",
                       stat_yrs_obj = "all",
                       stat_yrs = "all",
                       data,
                       risk_line = FALSE
) {
  
  data <- data[data$obj == obj & data$stat_yrs_obj %in% stat_yrs_obj &
                 data$stat_yrs %in% stat_yrs, ]
  
  p_theme <- theme_bw(base_size = 8, base_family = "sans") +
    theme(panel.spacing.x = unit(0, units = "cm"),
          strip.placement.y = "outside",
          strip.background.y = element_blank(),
          strip.text.y = element_text(size = 8),
          plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  p_pol_stats_SSB <- data %>% 
    filter(stat %in% c("SSB/B[MSY]")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
             show.legend = FALSE, colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    scale_y_continuous(trans = trans_from(), limits = c(SSB_min, SSB_max))
  
  p_pol_stats_SSB <- data %>% 
    filter(stat %in% c("SSB/B[MSY]")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
             show.legend = FALSE, colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    scale_y_continuous(trans = trans_from(), limits = c(SSB_min, SSB_max))
  
  p_pol_stats_F <- data %>% 
    filter(stat %in% c("F/F[MSY]")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(), limits = c(F_min, F_max))
  p_pol_stats_C <- data %>% 
    filter(stat %in% c("Catch/MSY")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(), limits = c(C_min, C_max))
  p_pol_stats_risk <- data %>% 
    filter(stat %in% c("B[lim]~risk")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = ifelse(isTRUE(risk_line), 0.05, 0), 
               linetype = "solid", size = 0.5, 
               colour = ifelse(isTRUE(risk_line), "red", "grey")) +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(0), limits = c(risk_min, risk_max))
  p_pol_stats_ICV <- data %>% 
    filter(stat %in% c("ICV")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(0), limits = c(ICV_min, ICV_max))
  p_pol_stats_fitness <- data %>% 
    filter(stat %in% c("fitness~value")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "") +
    theme_bw(base_size = 8, base_family = "sans") +
    theme(panel.spacing.x = unit(0, units = "cm"),
          strip.text.x = element_blank(),
          strip.placement.y = "outside",
          strip.background.y = element_blank(),
          strip.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.margin = unit(x = c(0, 3, 3, 3), units = "pt")) +
    scale_y_continuous(trans = trans_from(0), 
                       limits = c(fitness_min, fitness_max)#,
                       #breaks = c(0, -0.5, -1, -1.5), 
                       #minor_breaks = c(-0.25, -0.75, -1.25)
                       )
  p_pol_stats_comb <- plot_grid(p_pol_stats_SSB, p_pol_stats_F, p_pol_stats_C,
                                p_pol_stats_risk, p_pol_stats_ICV,
                                p_pol_stats_fitness,
                                ncol = 1, align = "v",
                                rel_heights = c(1.25, 1, 1, 1, 1, 2.1))
  return(p_pol_stats_comb)
}


### plots for MSY fitness function
plot_stats(obj = "MSY", stat_yrs_obj = "all", data = stats_pol,
           SSB_max = 1.5, F_max = 1.5, C_max = 1.5, risk_max = 1.5, 
           ICV_max = 1.5, fitness_min = -1.5)

ggsave(filename = "output/plots/PA/pol_GA_params_MSY.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/pol_GA_params_MSY.pdf",
      width = 17, height = 13, units = "cm", dpi = 600)

### PA fitness function
plot_stats(obj = "PA", stat_yrs_obj = "all", data = stats_pol,
           SSB_max = NA, F_max = 1, C_max = 1, risk_max = 1, 
           ICV_max = 1, fitness_min = NA, fitness_max = 0.5, 
           risk_line = TRUE)

ggsave(filename = "output/plots/PA/pol_GA_params_PA.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/pol_GA_params_PA.pdf",
      width = 17, height = 13, units = "cm", dpi = 600)



### MSY & PA fitness function
plot_stats(obj = "MSYPA", stat_yrs_obj = "all", data = stats_pol,
           SSB_max = NA, F_max = 1, C_max = 1, risk_max = 1, 
           ICV_max = 1, fitness_min = NA, fitness_max = NA, 
           risk_line = TRUE)

ggsave(filename = "output/plots/PA/pol_GA_params_MSYPA.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/pol_GA_params_MSYPA.pdf",
      width = 17, height = 13, units = "cm", dpi = 600)

### MSY & PA fitness function & last 10 years for stats
plot_stats_2(obj = "MSYPA", 
             stat_yrs = c("all", "last10"), 
             stat_yrs_obj = c("all", "last10"),
             data = stats_pol,
             SSB_max = NA, F_max = 1.2, C_max = 1.2, risk_max = 1, 
             ICV_max = 1, fitness_min = NA, fitness_max = NA, 
             risk_line = TRUE)

ggsave(filename = "output/plots/PA/pol_GA_params_MSYPA_last10.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/pol_GA_params_MSYPA_last10.pdf",
      width = 17, height = 13, units = "cm", dpi = 600)


### plot function, now for plotting several options
plot_stats_2 <- function(SSB_min = 0, SSB_max = NA,
                       F_min = 0, F_max = NA,
                       C_min = 0, C_max = NA,
                       risk_min = 0, risk_max = NA,
                       ICV_min = 0, ICV_max = NA,
                       fitness_min = NA, fitness_max = NA,
                       obj = "MSY",
                       stat_yrs_obj = "all",
                       stat_yrs = "all",
                       data,
                       risk_line = FALSE
) {
  
  data <- data[data$obj == obj & data$stat_yrs_obj %in% stat_yrs_obj &
                 data$stat_yrs %in% stat_yrs, ]
  
  data <- data %>%
    mutate(scenario2 = paste0("GA yrs: ", stat_yrs_obj, "\n",
                             "stat yrs: ", stat_yrs))
  
  p_theme <- theme_bw(base_size = 8, base_family = "sans") +
    theme(panel.spacing.x = unit(0, units = "cm"),
          strip.placement.y = "outside",
          strip.background.y = element_blank(),
          strip.text.y = element_text(size = 8),
          plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  p_pol_stats_SSB <- data %>% 
    filter(stat %in% c("SSB/B[MSY]")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario2,
               colour = scenario2)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
             show.legend = TRUE, colour = "black", size = 0.1) +
    scale_fill_discrete("") +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(legend.key.height = unit(1, "lines"),
          legend.key.width = unit(0.5, "lines")) +
    scale_y_continuous(trans = trans_from(), limits = c(SSB_min, SSB_max))
  
  p_pol_stats_F <- data %>% 
    filter(stat %in% c("F/F[MSY]")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario2,
               colour = scenario2)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(), limits = c(F_min, F_max))
  p_pol_stats_C <- data %>% 
    filter(stat %in% c("Catch/MSY")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario2,
               colour = scenario2)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(), limits = c(C_min, C_max))
  p_pol_stats_risk <- data %>% 
    filter(stat %in% c("B[lim]~risk")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario2,
               colour = scenario2)) +
    geom_hline(yintercept = ifelse(isTRUE(risk_line), 0.05, 0), 
               linetype = "solid", size = 0.5, 
               colour = ifelse(isTRUE(risk_line), "red", "grey")) +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(0), limits = c(risk_min, risk_max))
  p_pol_stats_ICV <- data %>% 
    filter(stat %in% c("ICV")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario2,
               colour = scenario2)) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(0), limits = c(ICV_min, ICV_max))
  p_pol_stats_fitness <- data %>% 
    filter(stat %in% c("fitness~value")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario2,
               colour = scenario2)) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "") +
    theme_bw(base_size = 8, base_family = "sans") +
    theme(panel.spacing.x = unit(0, units = "cm"),
          strip.text.x = element_blank(),
          strip.placement.y = "outside",
          strip.background.y = element_blank(),
          strip.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.margin = unit(x = c(0, 3, 3, 3), units = "pt")) +
    scale_y_continuous(trans = trans_from(0), 
                       limits = c(fitness_min, fitness_max)#,
                       #breaks = c(0, -0.5, -1, -1.5), 
                       #minor_breaks = c(-0.25, -0.75, -1.25)
                       )
  p_pol_stats_comb <- plot_grid(
    plot_grid(p_pol_stats_SSB + theme(legend.position = "none"), 
              p_pol_stats_F, p_pol_stats_C, p_pol_stats_risk, p_pol_stats_ICV,
              p_pol_stats_fitness,
              ncol = 1, align = "v", rel_heights = c(1.25, 1, 1, 1, 1, 2.1)),
    get_legend(p_pol_stats_SSB), rel_widths = c(1, 0.2), ncol = 2
  )
  return(p_pol_stats_comb)
}


### ------------------------------------------------------------------------ ###
### fitness penalty visualisation ####
### ------------------------------------------------------------------------ ###

penalty <- function(x, negative = FALSE, max = 5,
                    inflection = 0.06, steepness = 0.5e+3) {
  y <- max / (1 + exp(-(x - inflection)*steepness))
  if (isTRUE(negative)) y <- -y
  return(y)
}

p <- ggplot() +
  geom_function(fun = penalty, n = 1000) +
  geom_vline(xintercept = 0.05, colour = "red") +
  theme_bw(base_size = 8, base_family = "sans") +
  #scale_x_continuous(expand = c(0, 0)) +
  xlim(c(0, 1)) + 
  coord_cartesian(xlim = c(0, 0.3)) +
  labs(x = expression(B[lim]~risk),
       y = "fitness penalty")
p
ggsave(filename = "output/plots/PA/Blim_penalty_curve.png",
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/Blim_penalty_curve.pdf",
      width = 8.5, height = 6, units = "cm", dpi = 600)



### ------------------------------------------------------------------------ ###
### multiplier all stocks ####
### ------------------------------------------------------------------------ ###

stocks_subset <- stocks$stock[1:29]

mult_all <- foreach(fhist = c("one-way", "random"), .combine = bind_rows) %:%
  foreach(stock = stocks_subset, .combine = bind_rows) %do% {
    #browser()
    runs <- readRDS(paste0("output/500_50/uncertainty_cap/", fhist, "/", stock, "/",
                   "multiplier--obj_ICES_MSYPA_runs_last10.rds"))
    
    runs <- lapply(runs, function(x) {
      s <- as.data.frame(rbind(x$stats[1:11], x$stats[12:22]))
      names(s) <- rownames(x$stats)[1:11]
      tmp <- cbind(t(x$pars), s)
      tmp <- as.data.frame(lapply(as.data.frame(tmp), unlist))
      tmp$stat_yrs <- c("all", "last10")
      return(tmp)
    })
    runs <- do.call(rbind, runs)
    row.names(runs) <- NULL
    runs$fhist <- fhist
    runs$stock <- stock
    
    return(runs)
}
mult_all <- mult_all %>% 
    left_join(stocks[, c("stock", "k")]) %>%
    mutate(stock_k = paste0(stock, "~(italic(k)==", k, ")")) %>%
    mutate(stock_k = factor(stock_k, levels = unique(stock_k)))



mult_all %>% 
  filter(stat_yrs == "last10" & fhist == "random") %>%
  ggplot(aes(x = multiplier, y = risk_Blim)) +
  geom_line() +
  geom_hline(yintercept = 0.05, colour = "red") +
  facet_wrap(~ stock_k, labeller = "label_parsed", scales = "free") +
  coord_cartesian(xlim = c(0, 1))



### ------------------------------------------------------------------------ ###
### MSYPA fitness function - all stocks ####
### ------------------------------------------------------------------------ ###

fhist <- "one-way"
stocks_subset <- stocks$stock[1:29]
res <- foreach(stock = stocks_subset, .combine = bind_rows) %:%
  foreach(obj = c("MSYPA"), .combine = bind_rows) %:%
  foreach(params = c("default", "multiplier", "full_cap"), 
          .combine = bind_rows) %:%
  foreach(stat_yrs = c("all"), .combine = bind_rows) %do% {#browser()
    #browser()
    ### load data
    par_file <- switch(params,
      "default" = "multiplier", 
      "multiplier" = "multiplier", 
      "cap" = "upper_constraint-lower_constraint", 
      "cap_multiplier" = "multiplier-upper_constraint-lower_constraint", 
      "full" = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                      "interval-multiplier"), 
      "full_cap" = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                      "interval-multiplier-upper_constraint-lower_constraint")
    )
    obj_file <- switch(obj,
      "PA" = "obj_ICES_PA2",
      "MSY" = "obj_SSB_C_risk_ICV",
      "MSYPA" = "obj_ICES_MSYPA"
    )
    path <- paste0("output/500_50/uncertainty_cap/", fhist, "/", stock, "/",
                   par_file, "--", obj_file)
    path_runs <- paste0(path, "_runs", 
                        ifelse(stat_yrs == "all", "", paste0("_", stat_yrs)),
                        ".rds")
    path_res <- paste0(path, "_res", 
                        ifelse(stat_yrs == "all", "", paste0("_", stat_yrs)),
                        ".rds")
    ### use GA paper results for "full" GA
    # if (isTRUE(obj == "MSY" & params == "full")) {
    #   path <- paste0("output/500_50/ms/trial/", fhist, "/pol/",
    #                par_file, "--", obj_file)
    # }
    if (!file.exists(path_runs)) return(NULL)
    print("found something")
    ga_res <- readRDS(path_res)
    ga_runs <- readRDS(path_runs)
    
    ### optimised parameters
    if (isFALSE(params == "default")) {
      pars <- ga_solution(ga_res)
    } else {
      pars <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0)
    }
    pars[which(is.nan(pars))] <- Inf
    tmp <- as.data.frame(t(pars))
    names(tmp) <- c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
                    "exp_r", "exp_f", "exp_b", "interval", "multiplier",
                    "upper_constraint", "lower_constraint")
    #if (is.nan(tmp$upper_constraint)) tmp$upper_constraint <- Inf
    tmp$obj <- obj
    tmp$fhist <- fhist
    tmp$stat_yrs_obj <- stat_yrs
    tmp$ga_obj <- params
    tmp$stock <- stock
    
    ### stats
    par_scn <- pars
    if (isTRUE(length(which(is.na(par_scn))) > 0)) {
      par_scn <- par_scn[-which(is.na(par_scn))]
    }
    par_scn <- paste0(par_scn, collapse = "_")
    stats_tmp <- ga_runs[[par_scn]]
    stats_tmp <- as.data.frame(lapply(as.data.frame(t(stats_tmp$stats)), unlist))
    
    ### combine pars and stats
    stats_tmp <- cbind(tmp, stats_tmp)
    
    ### if different stat_yrs period used, extract also default stats
    if (isFALSE(stat_yrs == "all")) {
      stats_tmp <- rbind(stats_tmp, stats_tmp)
      stats_tmp$stat_yrs <- c("all", stat_yrs)
      stats_names <- c("risk_Blim", "risk_Bmsy", "risk_halfBmsy", 
                       "risk_collapse", "SSB", "Fbar", "Catch", "SSB_rel", 
                       "Fbar_rel", "Catch_rel", "ICV")
      stats_tmp[2, stats_names] <- stats_tmp[2, paste0(stats_names, "_", stat_yrs)]
      stats_tmp[, paste0(stats_names, "_", stat_yrs)] <- NULL
    } else {
      stats_tmp$stat_yrs <- "all"
      ### remove redundant stats
      stats_tmp[, grep(x = names(stats_tmp), pattern = "_last10")] <- NULL
    }
    
    ### recreate fitness
    yr_suffix <- ""
    if (isTRUE(obj == "PA")) {
      stats_tmp$fitness <- sapply(seq(nrow(stats_tmp)), function(x) {
        sum(stats_tmp[x, "Catch_rel"]) -
          sum(penalty(x = stats_tmp[x, "risk_Blim"], negative = FALSE, max = 5,
                      inflection = 0.06, steepness = 0.5e+3))
      })
    } else if (isTRUE(obj == "MSY")) {
      stats_tmp$fitness <- sapply(seq(nrow(stats_tmp)), function(x) {
        -sum(abs(stats_tmp[x, "SSB_rel"] - 1),
                          abs(stats_tmp[x, "Catch_rel"] - 1),
                          stats_tmp[x, "ICV"], 
                          stats_tmp[x, "risk_Blim"])
      })
    } else if (isTRUE(obj == "MSYPA")) {
      stats_tmp$fitness <- sapply(seq(nrow(stats_tmp)), function(x) {
        -sum(abs(stats_tmp[x, "SSB_rel"] - 1),
                          abs(stats_tmp[x, "Catch_rel"] - 1),
                          stats_tmp[x, "ICV"], 
                          penalty(x = stats_tmp[x, "risk_Blim"], 
                                  negative = FALSE, max = 5, 
                                  inflection = 0.06, steepness = 0.5e+3))
      })
    }
    return(stats_tmp)
}
res <- res %>%
  mutate(ga_obj = factor(ga_obj, levels = c("default", "multiplier", "full_cap"),
                         labels = c("default", "GA multiplier", 
                                    "GA all"))) %>% 
  left_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_k = paste0(stock, "~(italic(k)==", k, ")")) %>%
  mutate(stock_k = factor(stock_k, levels = unique(stock_k))) %>%
  mutate(stock = factor(stock, levels = stocks$stock))
  


saveRDS(res, file = "output/500_50/uncertainty_cap/all_stocks_mult_full.rds")
res <- readRDS("output/500_50/uncertainty_cap/all_stocks_mult_full.rds")

### plot
res_plot <- res %>%
  select(obj, fhist, stat_yrs_obj, stat_yrs, ga_obj, risk_Blim, SSB_rel, 
         Fbar_rel, Catch_rel, ICV, fitness, stock) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV, fitness), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel",
                                       "risk_Blim", "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value")))
stats_targets <- data.frame(stat = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                     "B[lim]~risk", "ICV", "fitness~value"),
                            target = c(1, 1, 1, 0, 0, NA))

p_theme <- theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

p_stats_SSB <- res_plot %>% 
  filter(stat %in% c("SSB/B[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = ga_obj,
             colour = ga_obj)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = TRUE, colour = "black", size = 0.1) +
  scale_fill_discrete("") +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  theme(legend.key.height = unit(1, "lines"),
        legend.key.width = unit(0.5, "lines")) +
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))

p_stats_F <- res_plot %>% 
  filter(stat %in% c("F/F[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = ga_obj,
             colour = ga_obj)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))
p_stats_C <- res_plot %>% 
  filter(stat %in% c("Catch/MSY")) %>%
  ggplot(aes(x = stock, y = value, fill = ga_obj,
             colour = ga_obj)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))
p_stats_risk <- res_plot %>% 
  filter(stat %in% c("B[lim]~risk")) %>%
  ggplot(aes(x = stock, y = value, fill = ga_obj,
             colour = ga_obj)) +
  geom_hline(yintercept = ifelse(isTRUE(TRUE), 0.05, 0), 
             linetype = "solid", size = 0.5, 
             colour = ifelse(isTRUE(TRUE), "red", "grey")) +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
  scale_y_continuous(trans = trans_from(0), limits = c(0, 1))
p_stats_ICV <- res_plot %>% 
  filter(stat %in% c("ICV")) %>%
  ggplot(aes(x = stock, y = value, fill = ga_obj,
             colour = ga_obj)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
  scale_y_continuous(trans = trans_from(0), limits = c(0, 1))
p_stats_fitness <- res_plot %>% 
  filter(stat %in% c("fitness~value")) %>%
  ggplot(aes(x = stock, y = value, fill = ga_obj,
             colour = ga_obj)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.margin = unit(x = c(0, 3, 3, 3), units = "pt")) +
  scale_y_continuous(trans = trans_from(0), 
                     limits = c(NA, NA)#,
                     #breaks = c(0, -0.5, -1, -1.5), 
                     #minor_breaks = c(-0.25, -0.75, -1.25)
                     )
p_stats_comb <- plot_grid(
  plot_grid(p_stats_SSB + theme(legend.position = "none"), 
            p_stats_F, p_stats_C, p_stats_risk, p_stats_ICV,
            p_stats_fitness,
            ncol = 1, align = "v", rel_heights = c(1.25, 1, 1, 1, 1, 1.5)),
  get_legend(p_stats_SSB), rel_widths = c(1, 0.2), ncol = 2
)

ggsave(filename = "output/plots/PA/all_GA_params_MSYPA.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/all_GA_params_MSYPA.pdf",
      width = 17, height = 13, units = "cm", dpi = 600)





### ------------------------------------------------------------------------ ###
### rfb-rule & rb-rule with +20 -30% uncertainty cap - multipliers ####
### ------------------------------------------------------------------------ ###

res_cap <- foreach(stock = stocks$stock[1:29], .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows) %:%
  #foreach(obj = c("MSYPA"), .combine = bind_rows) %:%
  #foreach(params = c("default", "multiplier"), .combine = bind_rows) %:%
  foreach(rule = c("rfb", "rb"), .combine = bind_rows) %do% {#browser()
  #foreach(stat_yrs = c("all"), .combine = bind_rows) %do% {
    #browser()
    ### load data
    par_file <- "multiplier-upper_constraint1.2-lower_constraint0.7"
    obj_file <- "obj_ICES_MSYPA"
    scenario <- switch(rule, 
                       "rfb" = "cap_20_30",
                       "rb" = "cap_20_30_rb")
    path <- paste0("output/500_50/", scenario, "/", fhist, "/", stock, "/",
                   par_file, "--", obj_file)
    path_runs <- paste0(path, "_runs_last10.rds")
    path_res <- paste0(path, "_res_last10.rds")
    if (!file.exists(path_runs)) return(NULL)
    print("found something")
    #ga_res <- readRDS(path_res)
    ga_runs <- readRDS(path_runs)
    
    stats <- lapply(ga_runs, function(x) {
      tmp <- cbind(as.data.frame(t(x$pars)),
                   as.data.frame(lapply(as.data.frame(t(x$stats)), unlist)))
      tmp <- rbind(tmp, tmp)
      tmp$stat_yrs <- c("all", "last10")
      stats_names <- c("risk_Blim", "risk_Bmsy", "risk_halfBmsy", 
                       "risk_collapse", "SSB", "Fbar", "Catch", "SSB_rel", 
                       "Fbar_rel", "Catch_rel", "ICV")
      tmp[2, stats_names] <- tmp[2, paste0(stats_names, "_last10")]
      tmp[, paste0(stats_names, "_last10")] <- NULL
      return(tmp)
    })
    stats <- do.call(rbind, stats)
    stats$stock <- stock
    stats$fhist <- fhist
    stats$rule <- rule
    return(stats)
}

res_cap <- res_cap %>%
  left_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_k = paste0(stock, "~(italic(k)==", k, ")")) %>%
  mutate(stock_k = factor(stock_k, levels = unique(stock_k))) %>%
  mutate(stock = factor(stock, levels = stocks$stock))
  


saveRDS(res_cap, file = "output/500_50/rfb-rb-cap2030.rds")
res_cap <- readRDS("output/500_50/rfb-rb-cap2030.rds")


### plot
res_cap <- res_cap %>%
  #filter(stat_yrs == "all") %>%
  select(multiplier, risk_Blim, SSB_rel, Fbar_rel, Catch_rel, ICV,
         stat_yrs, stock, fhist, rule, stock_k) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel",
                                       "risk_Blim", "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value")))
stats_targets <- data.frame(stat = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                     "B[lim]~risk", "ICV"),
                            target = c(1, 1, 1, 0, 0))
### plot full period rfb-rule
p <- res_cap %>% 
  filter(stat_yrs == "all" & rule == "rfb") %>%
  ggplot(aes(x = multiplier, y = value,
             colour = as.factor(fhist), linetype = as.factor(fhist))) +
  geom_line(size = 0.3) +
  geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.05),
             aes(yintercept = y), colour = "red") +
  facet_grid(stat ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  scale_linetype_discrete("fishing\nhistory") +
  scale_colour_discrete("fishing\nhistory") +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA)) +
  scale_x_continuous(breaks = c(0, 1, 2)#,
                     #expand = expansion(mult = c(0.1, 0.1))
                     )
p
ggsave(filename = "output/plots/defcap_rfb_mult.pdf",
       width = 40, height = 10, units = "cm")
ggsave(filename = "output/plots/defcap_rfb_mult.png", type = "cairo",
       width = 40, height = 10, units = "cm", dpi = 600)


res_cap <- readRDS("output/500_50/rfb-rb-cap2030.rds")
res_sub <- res_cap %>% 
  filter(stat_yrs == "all" & rule == "rfb" & fhist == "one-way") %>%
  filter(k <= 0.19) %>%
  group_by(multiplier) %>%
  summarise(risk = median(risk_Blim))
res_sub %>% 
  ggplot(aes(x = multiplier, y = risk)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  ylim(0, NA)



res_sub <- res_cap %>% 
  filter(stat_yrs == "all" & rule == "rfb" & fhist == "one-way") %>%
  filter(k > 0.19 & k <= 0.32) %>%
  group_by(multiplier) %>%
  summarise(risk = median(risk_Blim))
res_sub %>% 
  ggplot(aes(x = multiplier, y = risk)) +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  ylim(0, NA)



res_sub <- res_cap %>% 
  filter(stat_yrs == "all" & rule == "rfb" & fhist %in% c("random","one-way")) %>%
  filter(k <= 0.32) %>%
  # mutate(group = ifelse(k > 0.19, "medium-italic(k)", "low-italic(k)")) %>%
  mutate(group = ifelse(TRUE, "medium-italic(k)", "low-italic(k)")) %>%
  group_by(multiplier, group) %>%
  summarise(risk = median(risk_Blim))
res_sub %>% 
  ggplot(aes(x = multiplier, y = risk)) +
  geom_line() +
  geom_hline(yintercept = 0.05, colour = "red") +
  ylim(0, NA) +
  facet_wrap(~ group, labeller = "label_parsed") +
  labs(x = "multiplier", y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/defcap_rfb_mult_groups.png", type = "cairo",
       width = 17, height = 8, units = "cm", dpi = 600)




NumDiffTable <- function(par, val, method = "central") {
  ### forward difference approximation
  if (identical(method, "forward")) {
    diff <- sapply(seq_along(par)[1:(length(par) - 1)], function(x) {
      (val[x + 1] - val[x])/(par[x + 1] - par[x])
    })
    diff <- c(diff, NA)
  ### backward difference approximation
  } else if (identical(method, "backward")) {
    diff <- sapply(seq_along(par)[2:(length(par))], function(x) {
      (val[x] - val[x - 1])/(par[x] - par[x - 1])
    })
    diff <- c(NA, diff)
  ### use central difference approximation
  } else if (identical(method, "central")) {
    diff <- sapply(seq_along(par)[2:(length(par) - 1)], function(x) {
      (val[x + 1] - val[x - 1])/(par[x + 1] - par[x - 1])
    })
    ### first/last values
    diff0 <- (val[2] - val[1])/(par[2] - par[1])
    diff1 <- (val[length(val)] - val[length(val) - 1]) /
      (par[length(val)] - par[length(val) - 1])
    diff <- c(diff0, diff, diff1)
  }
  return(diff)
}

res_sub <- res_cap %>% 
  filter(stat_yrs == "all" & rule == "rfb" & fhist %in% c("random","one-way")) %>%
  filter(k <= 0.32) %>%
  mutate(group = ifelse(k > 0.19, "medium-italic(k)", "low-italic(k)")) %>%
  #mutate(group = ifelse(TRUE, "medium-italic(k)", "low-italic(k)")) %>%
  group_by(multiplier, group) %>%
  summarise(risk = median(risk_Blim)) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(risk_diff = NumDiffTable(par = multiplier, val = risk,
                                       method = "central")) %>%
  pivot_longer(c(risk, risk_diff)) %>%
  mutate(name = factor(name, levels = c("risk", "risk_diff"),
                       labels = c("B[lim]~risk",
                                  "B[lim]~risk~bold('\\'')")))
p <- res_sub %>%
  ggplot(aes(x = multiplier, y = value)) +
  geom_smooth(method = "loess", span = 0.2, n = 10000, colour = "grey",
              level = 0, size = 0.5) +
  geom_point(size = 0.5, stroke = 0) +
  facet_grid(name ~ group, scales = "free", labeller = label_parsed) +
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
  # geom_vline(xintercept = 0.705, 
  #            colour = "black", linetype = "dashed", size = 0.4) +
  theme_bw(base_size = 8) +
  labs(x = "catch rule multiplier", y = "") #+
p


res_sub <- res_cap %>% 
  filter(stat_yrs == "all" & rule == "rfb" & fhist %in% c("random","one-way")) %>%
  filter(k <= 0.32) %>%
  mutate(group = ifelse(k > 0.19, "medium-italic(k)", "low-italic(k)")) %>%
  #mutate(group = ifelse(TRUE, "medium-italic(k)", "low-italic(k)")) %>%
  group_by(group, stock, fhist) %>%
  #summarise(risk = median(risk_Blim)) %>%
  # ungroup() %>%
  # group_by(group) %>%
  mutate(risk_diff = NumDiffTable(par = multiplier, val = risk_Blim,
                                       method = "central")) %>%
  pivot_longer(c(risk_Blim, risk_diff)) %>%
  mutate(name = factor(name, levels = c("risk", "risk_Blim"),
                       labels = c("B[lim]~risk",
                                  "B[lim]~risk~bold('\\'')")))
p <- res_sub %>%
  ggplot(aes(x = multiplier, y = value, group = interaction(stock, fhist))) +
  geom_smooth(method = "loess", span = 0.2, n = 10000, colour = "grey",
              level = 0, size = 0.5) +
  geom_point(size = 0.5, stroke = 0) +
  facet_grid(name ~ group, scales = "free", labeller = label_parsed) +
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
  # geom_vline(xintercept = 0.705, 
  #            colour = "black", linetype = "dashed", size = 0.4) +
  theme_bw(base_size = 8) +
  labs(x = "catch rule multiplier", y = "") #+
p
ggsave(filename = "output/plots/defcap_rfb_mult_groups_diff.png", type = "cairo",
       width = 17, height = 12, units = "cm", dpi = 600)

res_sub <- res_cap %>% 
  filter(stat_yrs == "all" & rule == "rfb" & fhist %in% c("random","one-way")) %>%
  filter(k <= 0.32) %>%
  mutate(group = ifelse(k > 0.19, "medium-italic(k)", "low-italic(k)")) %>%
  #mutate(group = ifelse(TRUE, "medium-italic(k)", "low-italic(k)")) %>%
  group_by(group, stock, fhist) %>%
  #summarise(risk = median(risk_Blim)) %>%
  # ungroup() %>%
  # group_by(group) %>%
  mutate(risk_diff = NumDiffTable(par = multiplier, val = risk_Blim,
                                       method = "central")) %>%
  ungroup() %>%
  group_by(group, multiplier) %>%
  summarise(risk_Blim = median(risk_Blim),
            risk_Blim_diff = median(risk_diff)) %>%
  pivot_longer(c(risk_Blim, risk_Blim_diff)) %>%
  mutate(name = factor(name, levels = c("risk_Blim", "risk_Blim_diff"),
                       labels = c("B[lim]~risk",
                                  "B[lim]~risk~bold('\\'')")))
p <- res_sub %>%
  ggplot(aes(x = multiplier, y = value, group = interaction(stock, fhist))) +
  geom_smooth(method = "loess", span = 0.2, n = 10000, colour = "grey",
              level = 0, size = 0.5) +
  geom_point(size = 0.5, stroke = 0) +
  facet_grid(name ~ group, scales = "free", labeller = label_parsed) +
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
  # geom_vline(xintercept = 0.705, 
  #            colour = "black", linetype = "dashed", size = 0.4) +
  theme_bw(base_size = 8) +
  labs(x = "catch rule multiplier", y = "") #+
p
ggsave(filename = "output/plots/defcap_rfb_mult_groups_diff_median.png", type = "cairo",
       width = 17, height = 12, units = "cm", dpi = 600)




### ------------------------------------------------------------------------ ###
### r(f)b-rule with +20 -30% cap when above Itrigger - multipliers ####
### ------------------------------------------------------------------------ ###
### basis for WKLIFE X decision on multipliers

### load results
res_cap <- foreach(stock = stocks$stock[1:29], .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows) %:%
  foreach(rule = c("rfb", "rb"), .combine = bind_rows) %do% {
    #browser()
    ### load data
    par_file <- "multiplier-upper_constraint1.2-lower_constraint0.7"
    obj_file <- "obj_ICES_MSYPA"
    scenario <- switch(rule, 
                       "rfb" = "cap_20_30_rfb_b",
                       "rb" = "cap_20_30_rb_b")
    path <- paste0("output/500_100/", scenario, "/", fhist, "/", stock, "/",
                   par_file, "--", obj_file)
    path_runs <- paste0(path, "_runs.rds")
    if (!file.exists(path_runs)) return(NULL)
    print("found something")
    ga_runs <- readRDS(path_runs)
    
    stats <- lapply(ga_runs, function(x) {
      tmp <- cbind(as.data.frame(t(x$pars)),
                   as.data.frame(lapply(as.data.frame(t(x$stats)), unlist)))
      return(tmp)
    })
    stats <- do.call(rbind, stats)
    stats$stock <- stock
    stats$fhist <- fhist
    stats$rule <- rule
    rownames(stats) <- NULL
    return(stats)
}
### split reporting periods
smry_stats <- c("risk_Blim", "risk_Bmsy", "risk_halfBmsy", "risk_collapse", 
                "SSB", "Fbar", "Catch", "SSB_rel", "Fbar_rel", "Catch_rel", "ICV")
periods <- c("first10", "41to50", "last10", "firsthalf", "lastfhalf", "11to50")
res_cap <- lapply(split(res_cap, seq(nrow(res_cap))), function(x) {
  stats_more <- t(sapply(periods, function(y) {
    as.numeric(x[grep(x = names(x), pattern = y)])
  }))
  stats <- rbind(all = as.numeric(x[smry_stats]), stats_more)
  stats <- as.data.frame(stats)
  colnames(stats) <- smry_stats
  stats$stat_yrs <- rownames(stats)
  rownames(stats) <- NULL
  ### add scenario definition
  stats <- cbind(x[c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
                     "exp_r", "exp_f", "exp_b", "interval", "multiplier",
                     "upper_constraint", "lower_constraint",
                     "stock", "fhist", "rule")],
                 stats)
  return(stats)
})
res_cap <- do.call(rbind, res_cap)
### add k
res_cap <- res_cap %>%
  left_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_k = paste0(stock, "~(italic(k)==", k, ")")) %>%
  mutate(stock_k = factor(stock_k, levels = unique(stock_k))) %>%
  mutate(stock = factor(stock, levels = stocks$stock))
  
saveRDS(res_cap, file = "output/500_100/rfb-rb-cap2030_b.rds")
res_cap <- readRDS("output/500_100/rfb-rb-cap2030_b.rds")


### rb-rule ####
### plot
res_plot <- res_cap %>%
  select(multiplier, risk_Blim, SSB_rel, Fbar_rel, Catch_rel, ICV,
         stat_yrs, stock, fhist, rule, stock_k) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel",
                                       "risk_Blim", "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value")))
### plot full period rb-rule
res_plot %>% 
  filter(stat_yrs == "all" & rule == "rb") %>%
  ggplot(aes(x = multiplier, y = value,
             colour = as.factor(fhist), linetype = as.factor(fhist))) +
  geom_line(size = 0.3) +
  geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.05),
             aes(yintercept = y), colour = "red") +
  facet_grid(stat ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  scale_linetype_discrete("fishing\nhistory") +
  scale_colour_discrete("fishing\nhistory") +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA)) +
  scale_x_continuous(breaks = c(0, 1, 2)#,
                     #expand = expansion(mult = c(0.1, 0.1))
                     )
ggsave(filename = "output/plots/cap_2030_b/all_stocks_stats_rb.pdf",
       width = 40, height = 10, units = "cm")
ggsave(filename = "output/plots/cap_2030_b/all_stocks_stats_rb.png", type = "cairo",
       width = 40, height = 10, units = "cm", dpi = 600)


### medians of stock groups
res_plot <- res_cap %>%
  filter(stat_yrs == "all" & rule == "rb") %>%
  mutate(ICV = ifelse(multiplier == 0 & ICV == 1, NA, ICV)) %>%
  # filter(stock %in% c("ang3", "rjc2", "smn", "wlf", "meg", "lin", "rjc", "syc",
  #                     "sdv", "ang", "ang2", "pol", "had", "nep", "mut", "sbb",
  #                     "ple", "syc2", "arg", "tur")) %>%
  #filter(stock %in% stocks$stock[1:29]) %>%
  select(multiplier, risk_Blim, SSB_rel, Fbar_rel, Catch_rel, ICV,
         fhist, stock, stock_k, k)
res_mult <- res_plot %>% 
  group_by(multiplier) %>%
  summarise(risk_Blim = median(risk_Blim),
            SSB_rel = median(SSB_rel),
            Fbar_rel = median(Fbar_rel),
            Catch_rel = median(Catch_rel),
            ICV = median(ICV)) %>%
  mutate(stock = "median")
res_mult %>%
  filter(risk_Blim >= 0.04 & risk_Blim <= 0.07) %>% print(n = Inf)
res_plot <- res_plot %>%
  full_join(res_mult)
res_plot <- res_plot %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel",
                                       "risk_Blim", "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value"))) %>%
  mutate(group = ifelse(stock == "median", "median", NA),
         group = ifelse(stock != "median" & fhist == "one-way", "one-way", group),
         group = ifelse(stock != "median" & fhist == "random", "random", group))

res_plot %>% 
  ggplot(aes(x = multiplier, y = value,
             colour = group, 
             alpha = group,
             group = interaction(stock, fhist))) +
  geom_line(size = 0.3) +
  geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.05),
             aes(yintercept = y), colour = "red") +
  facet_wrap(~ stat, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  #scale_colour_manual("", values = c("TRUE" = "black", "FALSE" = "blue")) +
  scale_colour_manual("", values = c("median" = "black", "one-way" = "blue",
                                     "random" = "green")) +
  scale_alpha_manual("", values = c("median" = 1, "one-way" = 0.2,
                                     "random" = 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA)) +
  scale_x_continuous(breaks = c(0, 0.5, 1))
ggsave(filename = "output/plots/cap_2030_b/multiplier_rb.pdf",
       width = 17, height = 10, units = "cm")
ggsave(filename = "output/plots/cap_2030_b/multiplier_rb.png", 
       type = "cairo", width = 17, height = 10, units = "cm", dpi = 600)


### rfb-rule ####
### plot
res_plot <- res_cap %>%
  select(multiplier, risk_Blim, SSB_rel, Fbar_rel, Catch_rel, ICV,
         stat_yrs, stock, fhist, rule, stock_k) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel",
                                       "risk_Blim", "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value")))
### plot full period rfb-rule
res_plot %>% 
  filter(stat_yrs == "all" & rule == "rfb") %>%
  ggplot(aes(x = multiplier, y = value,
             colour = as.factor(fhist), linetype = as.factor(fhist))) +
  geom_line(size = 0.3) +
  geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.05),
             aes(yintercept = y), colour = "red") +
  facet_grid(stat ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  scale_linetype_discrete("fishing\nhistory") +
  scale_colour_discrete("fishing\nhistory") +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA)) +
  scale_x_continuous(breaks = c(0, 1, 2)#,
                     #expand = expansion(mult = c(0.1, 0.1))
                     )
ggsave(filename = "output/plots/cap_2030_b/all_stocks_stats_rfb.pdf",
       width = 40, height = 10, units = "cm")
ggsave(filename = "output/plots/cap_2030_b/all_stocks_stats_rfb.png", type = "cairo",
       width = 40, height = 10, units = "cm", dpi = 600)



### medians
res_plot <- res_cap %>%
  filter(stat_yrs == "all" & rule == "rfb") %>%
  mutate(ICV = ifelse(multiplier == 0 & ICV == 1, NA, ICV)) %>%
  # filter(stock %in% c("ang3", "rjc2", "smn", "wlf", "meg", "lin", "rjc", "syc",
  #                     "sdv", "ang", "ang2", "pol", "had", "nep", "mut", "sbb",
  #                     "ple", "syc2", "arg", "tur")) %>%
  filter(stock %in% stocks$stock[1:20]) %>%
  select(multiplier, risk_Blim, SSB_rel, Fbar_rel, Catch_rel, ICV,
         fhist, stock, stock_k, k) %>%
  mutate(group_k = ifelse(k < 0.2, "low-k", "medium-k"))
res_mult <- res_plot %>% 
  group_by(group_k, multiplier) %>%
  summarise(risk_Blim = median(risk_Blim),
            SSB_rel = median(SSB_rel),
            Fbar_rel = median(Fbar_rel),
            Catch_rel = median(Catch_rel),
            ICV = median(ICV)) %>%
  mutate(stock = "median")
res_mult %>%
  filter(risk_Blim >= 0.04 & risk_Blim <= 0.07) %>% print(n = Inf)
res_plot <- res_plot %>%
  full_join(res_mult)
res_plot <- res_plot %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel",
                                       "risk_Blim", "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value"))) %>%
  mutate(group = ifelse(stock == "median", "median", NA),
         group = ifelse(stock != "median" & fhist == "one-way", "one-way", group),
         group = ifelse(stock != "median" & fhist == "random", "random", group))

res_plot %>% 
  ggplot(aes(x = multiplier, y = value,
             colour = group, 
             alpha = group,
             group = interaction(stock, fhist))) +
  geom_line(size = 0.3) +
  geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.05),
             aes(yintercept = y), colour = "red") +
  facet_grid(stat ~ group_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  #scale_colour_manual("", values = c("TRUE" = "black", "FALSE" = "blue")) +
  scale_colour_manual("", values = c("median" = "black", "one-way" = "blue",
                                     "random" = "green")) +
  scale_alpha_manual("", values = c("median" = 1, "one-way" = 0.2,
                                     "random" = 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA)) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2))
ggsave(filename = "output/plots/cap_2030_b/multiplier_rfb_fixed.pdf",
       width = 17, height = 10, units = "cm")
ggsave(filename = "output/plots/cap_2030_b/multiplier_rfb_fixed.png", 
       type = "cairo", width = 17, height = 10, units = "cm", dpi = 600)
### zoom
res_plot %>% 
  filter(multiplier >= 0.7 & multiplier <= 1.05) %>%
  filter(stat == "B[lim]~risk") %>%
  ggplot(aes(x = multiplier, y = value,
             colour = group, 
             alpha = group,
             group = interaction(stock, fhist))) +
  geom_line(size = 0.3) +
  geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.05),
             aes(yintercept = y), colour = "red") +
  facet_grid(stat ~ group_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  #scale_colour_manual("", values = c("TRUE" = "black", "FALSE" = "blue")) +
  scale_colour_manual("", values = c("median" = "black", "one-way" = "blue",
                                     "random" = "green")) +
  scale_alpha_manual("", values = c("median" = 1, "one-way" = 0.2,
                                     "random" = 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA))# +
  #scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2))
ggsave(filename = "output/plots/cap_2030_b/multiplier_rfb_zoom_fixed.pdf",
       width = 17, height = 10, units = "cm")
ggsave(filename = "output/plots/cap_2030_b/multiplier_rfb_zoom_fixed.png", 
       type = "cairo", width = 17, height = 10, units = "cm", dpi = 600)

### plot multiplier vs. k
m_vs_k <- res_cap %>%
  filter(stat_yrs == "all" & rule == "rfb") %>%
  group_by(fhist, stock) %>%
  filter(risk_Blim <= 0.05) %>%
  filter(multiplier == max(multiplier)) %>%
  filter(stock %in% stocks$stock[1:20])
p_m_vs_k <- m_vs_k %>%
  ggplot(aes(x = k, y = multiplier)) +
  geom_line(aes(group = stock), colour = "grey", size = 0.3) +
  geom_point(aes(colour = fhist), size = 0.4) +
  scale_color_brewer("fishing\nhistory", palette = "Set1") +
  labs(x = expression(italic(k)~"[year"^{-1}*"]"), 
       y = expression("multiplier ("*italic(m)*")")) +
  ylim(c(0, NA)) + xlim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        legend.key.width = unit(0.5, "lines"))
p_m_vs_k
ggsave(filename = "output/plots/cap_2030_b/multiplier_rfb_k.pdf",
       width = 10, height = 6, units = "cm")
ggsave(filename = "output/plots/cap_2030_b/multiplier_rfb_k.png", 
       type = "cairo", width = 10, height = 6, units = "cm", dpi = 600)
### add smoother
p_m_vs_k + geom_smooth(size = 0.3)
ggsave(filename = "output/plots/cap_2030_b/multiplier_rfb_k_smoother.png", 
       type = "cairo", width = 10, height = 6, units = "cm", dpi = 600)
### add linear regression
p_m_vs_k + geom_smooth(method = "lm", size = 0.3)
ggsave(filename = "output/plots/cap_2030_b/multiplier_rfb_k_lm.png", 
       type = "cairo", width = 10, height = 6, units = "cm", dpi = 600)
### add robust linear regression
library(MASS)
p_m_vs_k + geom_smooth(method = "rlm", size = 0.3)
ggsave(filename = "output/plots/cap_2030_b/multiplier_rfb_k_rlm.png", 
       type = "cairo", width = 10, height = 6, units = "cm", dpi = 600)
# lm(data = m_vs_k[, c("k", "multiplier")], formula = multiplier ~ k)
# m_rlm <- MASS::rlm(formula = m_vs_k$multiplier ~ m_vs_k$k)
# Call:
#   rlm(formula = m_vs_k$multiplier ~ m_vs_k$k)
# Converged in 7 iterations
# 
# Coefficients:
#   (Intercept)    m_vs_k$k 
# 1.0473833  -0.6469253 
# 
# Degrees of freedom: 40 total; 38 residual
# Scale estimate: 0.0755
#
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.205060 -0.047844 -0.002702  0.052150  0.630248 
# 
# Coefficients:
#   Value   Std. Error t value
# (Intercept)  1.0474  0.0482    21.7240
# m_vs_k$k    -0.6469  0.2631    -2.4592

### steps
p_m_vs_k +
  geom_line(data = data.frame(k = c(0.08, 0.195), multiplier = c(0.95, 0.95)),
            colour = "black") +
  geom_line(data = data.frame(k = c(0.2, 0.32), multiplier = c(0.9, 0.9)),
            colour = "black")
ggsave(filename = "output/plots/cap_2030_b/multiplier_rfb_k_step.png", 
       type = "cairo", width = 10, height = 6, units = "cm", dpi = 600)

### step + robust regression
p_m_vs_k +
  geom_smooth(method = "rlm", size = 0.3, se = FALSE) +
  geom_line(data = data.frame(k = c(0.08, 0.195), multiplier = c(0.95, 0.95)),
            colour = "black") +
  geom_line(data = data.frame(k = c(0.2, 0.32), multiplier = c(0.9, 0.9)),
            colour = "black")
ggsave(filename = "output/plots/cap_2030_b/multiplier_rfb_k_step_rlm.png", 
       type = "cairo", width = 10, height = 6, units = "cm", dpi = 600)
ggsave(filename = "output/plots/cap_2030_b/multiplier_rfb_k_step_rlm.pdf",
       width = 10, height = 6, units = "cm")



### check 10%
m_vs_k_10 <- res_cap %>%
  filter(stat_yrs == "all" & rule == "rfb") %>%
  group_by(fhist, stock) %>%
  filter(risk_Blim <= 0.10) %>%
  filter(multiplier == max(multiplier)) %>%
  filter(stock %in% stocks$stock[1:20])
m_vs_k_10 %>%
  ggplot(aes(x = k, y = multiplier)) +
  geom_line(aes(group = stock), colour = "grey", size = 0.3) +
  geom_point(aes(colour = fhist), size = 0.4) +
  scale_color_brewer("fishing\nhistory", palette = "Set1") +
  labs(x = expression(italic(k)~"[year"^{-1}*"]"), 
       y = expression("multiplier ("*italic(m)*")")) +
  ylim(c(0, NA)) + xlim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        legend.key.width = unit(0.5, "lines")) +
  geom_smooth(method = MASS::rlm, size = 0.3, se = FALSE)



### combine risk and multiplier plots
p_risk <- res_plot %>% 
  filter(multiplier >= 0.7 & multiplier <= 1.05) %>%
  filter(stat == "B[lim]~risk") %>%
  ggplot(aes(x = multiplier, y = value,
             colour = group, 
             alpha = group,
             group = interaction(stock, fhist))) +
  geom_line(size = 0.1) +
  geom_line(data = res_plot %>% 
              filter(multiplier >= 0.7 & multiplier <= 1.05) %>%
              filter(stat == "B[lim]~risk" & group == "median"),
            aes(x = multiplier, y = value), colour = "black", size = 0.4, 
            show.legend = FALSE) + 
  geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.05),
             aes(yintercept = y), colour = "black", linetype = "2121") +
  facet_grid(stat ~ group_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  scale_colour_manual("", values = c("median" = "black", "one-way" = "#E41A1C",
                                     "random" = "#377EB8")) +
  scale_alpha_manual("", values = c("median" = 1, "one-way" = 1,
                                    "random" = 1)) +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6),
        legend.position = c(0.7, 0.85),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA))
p_mult <- m_vs_k %>%
  ggplot(aes(x = k, y = multiplier)) +
  geom_line(aes(group = stock), colour = "grey", size = 0.3) +
  geom_point(aes(colour = fhist, shape = fhist), size = 0.7) +
  scale_color_brewer("fishing\nhistory", palette = "Set1") +
  scale_shape("fishing\nhistory") +
  labs(x = expression(italic(k)~"[year"^{-1}*"]"), 
       y = expression("multiplier")) +
  ylim(c(0, NA)) + xlim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        legend.key.width = unit(0.5, "lines"),
        legend.position = c(0.8, 0.85),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  geom_line(data = data.frame(k = c(0.08, 0.195), multiplier = c(0.95, 0.95)),
            colour = "black") +
  geom_line(data = data.frame(k = c(0.2, 0.32), multiplier = c(0.9, 0.9)),
            colour = "black")
plot_grid(p_risk, p_mult, labels = c("(a)", "(b)"), label_size = 9,
          nrow = 1, rel_widths = c(1.5, 1), align = "v", axis = "b")
ggsave(filename = "output/plots/cap_2030_b/multiplier_justification.pdf",
       width = 17, height = 6, units = "cm")
ggsave(filename = "output/plots/cap_2030_b/multiplier_justification.png", 
       type = "cairo", width = 17, height = 6, units = "cm", dpi = 600)
