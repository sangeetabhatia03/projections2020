library (readr)
library(dplyr)
library(ggplot2)
library(purrr)
library(drake)
library(incidence)
library(EpiEstim)
library(epitrix)
library(projections)
library(cowplot)

r_samples <- function(res, n = 1000) {

    in_last_window <- tail(res, 1)
    mean_r <- in_last_window$`Mean(R)`
    std_r <- in_last_window$`Std(R)`
    reparam <- epitrix::gamma_mucv2shapescale(
            mu = mean_r,
            cv = std_r / mean_r
    )
    out <- rgamma(
        n = n,
        shape = reparam$shape,
        scale = reparam$scale
    )
    out

}

windows_endpoints <- function(ndays, tw) {

    t_start <- 2:(ndays - tw)
    t_end <- t_start + tw

    list(
        t_start = t_start,
        t_end = t_end
    )

}


fix_manually <- function(ts) {

    ts$cum_cases[ts$`Country/Region` == "Japan" &
                 ts$date == "2020-01-23"] <- 2
    ## is 45 in the ts, checked against wikipedia
    ts$cum_cases[ts$`Country/Region` == "Japan" &
                 ts$date == "2020-02-06"] <- 25



    ts

}

source("params.R")
source("plan.R")
drake_config(plan)
