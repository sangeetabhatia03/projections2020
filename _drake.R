library (readr)
library(dplyr)
library(ggplot2)
library(purrr)
library(drake)
library(incidence)
library(EpiEstim)
library(epitrix)
library(projections)

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

cum_epicurve <- function(ts) {

    ggplot(ts,
           aes(date, cum_cases, col = `Country/Region`)) +
        geom_line() +
        theme(legend.position = "none") +
        theme_classic()

}

inurl <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv"

countries_considered <- c(
    "Italy", "Japan", "Iran (Islamic Republic of)",
    "France", "UK", "Spain"
)
mean_si <- 3.96
std_si <- 4.75
size <- 0.16
time_window <- c(4, 6, 8)
names(time_window) <- time_window
date_to_project_from <- "2020-03-03"
n_sim <- 1000
n_days <- 10

SItrunc <- 20

SI_Distr <- sapply(
  0:SItrunc,
  function(e) EpiEstim::discr_si(e, mean_si, std_si)
)

SI_Distr <- SI_Distr / sum(SI_Distr)

source("plan.R")
drake_config(plan)
