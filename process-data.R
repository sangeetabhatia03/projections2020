## Data from here:
## https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv
library(dplyr)
library(ggplot2)
library(purrr)

ts <- readr::read_csv("john_hopkins_timeseries.csv")
ts <- tidyr::gather(ts, key = date, value = cum_cases, `1/22/20`:`3/10/20`)
ts$date <- lubridate::mdy(ts$date)

ts$cum_cases[ts$`Country/Region` == "Japan" & ts$date == "2020-01-23"] <- 2
## is 45 in the ts, checked against wikipedia
ts$cum_cases[ts$`Country/Region` == "Japan" & ts$date == "2020-02-06"] <- 25


dplyr::filter(ts, `Country/Region` != "Mainland China") %>%
ggplot(aes(date, cum_cases, col = `Country/Region`)) +
    geom_line() +
    theme(legend.position="none")

## Daily Incidence

by_country <- split(ts, ts$`Country/Region`)
countries <- c("Italy", "Japan", "Iran (Islamic Republic of)", "France")
interesting <- by_country[countries]

incidence <- purrr::map(

    interesting,
    function(x) {
        x$cases <- c(0, diff(x$cum_cases))
        if (any(x$cases < 0)) message(x$`Country/Region`[1])
        x
    }

)


## SI
## mean serial interval = 4.5 days
##
### 1. Estimating R_t with EpiEstim
si_mean <- 3.96
si_std <- 4.75
CV_SI <- si_std / si_mean
SItrunc <- 60

SI_Distr <- sapply(
  0:SItrunc,
  function(e) EpiEstim::DiscrSI(e, si_mean, si_mean * CV_SI)
)

SI_Distr <- SI_Distr / sum(SI_Distr)
## # serial interval estimate used: mean = 3.96, sd =  4.75
## # from Du et al. from The University of Texas at Austin and University of Hong Kong
## # (https://www.medrxiv.org/content/10.1101/2020.01.28.20019299v4.full.pdf)
## # (https://www.medrxiv.org/content/10.1101/2020.02.19.20025452v2.full.pdf)
time_window <- 7

r0s <- purrr::map(
    incidence,
    function(df) {
        I <- df$cases
        start <- 2:(nrow(x) - time_window)
        end <- start + time_window
        res <- EpiEstim::EstimateR(
                             I,
                             T.Start = start,
                             T.End = end,
                             method = "NonParametricSI",
                             SI.Distr = SI_Distr,
                             plot = FALSE,
                             CV.Posterior = 1,
                             Mean.Prior = 1,
                             Std.Prior = 0.5
                         )
        res$R$date <- x$date[end]
        res$R
    }
)





## Projections
projections <- purrr::map2(
    incidence, r0s,
    function(df, r0) {
        I <- df$cases
        incid <- project(
            incid = matrix(I, ncol = 1),
            R = matrix(r0[["Median(R)"]][nrow(r0)], ncol = 1),
            si = SI_Distr,
            pij = matrix(1, nrow = 1, ncol = 1),
            n_days = 7
        )
        out <- data.frame(
            date = seq(
                from = tail(df$date, 1) + 1,
                length.out = 7,
                by = "1 day"
            ),
            incid = incid
        )
        out

    }
)


for (country in names(incidence)) {

    p <- ggplot(incidence[[country]], aes(date, cases)) +
        geom_col() +
        geom_line(data = projections[[country]], aes(date, incid)) +
        theme_classic()
    outfile <- snakecase::to_snake_case(country)
    ggsave(glue::glue("{outfile}.pdf"), p)

}




