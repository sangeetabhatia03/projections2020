## Data from here:
## https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv
library(dplyr)
library(ggplot2)
library(purrr)
source("project.R")
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
countries <- c(
    "Italy", "Japan", "Iran (Islamic Republic of)", "France", "UK", "Spain"
)
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
SItrunc <- 20

SI_Distr <- sapply(
  0:SItrunc,
  function(e) EpiEstim::DiscrSI(e, si_mean, si_mean * CV_SI)
)

SI_Distr <- SI_Distr / sum(SI_Distr)

size <- 0.16
## # serial interval estimate used: mean = 3.96, sd =  4.75
## # from Du et al. from The University of Texas at Austin and University of Hong Kong
## # (https://www.medrxiv.org/content/10.1101/2020.01.28.20019299v4.full.pdf)
## # (https://www.medrxiv.org/content/10.1101/2020.02.19.20025452v2.full.pdf)
time_window <- c(2, 4, 6, 8)
names(time_window) <- time_window

r0s <- purrr::map_dfr(
    incidence,
    function(df) {
        out <- map_dfr(
            time_window,
            function(tw) {
                I <- df$cases
                start <- 2:(nrow(df) - tw)
                end <- start + tw
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
                res$R$date <- df$date[end]
                res$R
            },
          .id = "time_window"
        )

    }, .id = "country"
)

####
#### Uncertain SI #########
##### Crude way of doing this for now ################
##### Look at various estimates of SI and get mean, and s.d. of mean
#####
## mean_estimates <- c(
##     7.5, 4.41, 4.0, 4.6, 3.96, 5.1, 3.3, 6.3, 4.56, 4.22, 6.37, 5.9, 5.2, 3.95
## )
## sd_mean <- sd(mean_estimates)
## sd_estimates <- c(
##     3.17, 4.75, 4.8
## )
## sd_sd <- sd(sd_estimates)

## r0s <- purrr::map_dfr(
##     incidence,
##     function(df) {
##         out <- map_dfr(
##             time_window,
##             function(tw) {
##                 I <- df$cases
##                 start <- 2:(nrow(df) - tw)
##                 end <- start + tw
##                 res <- EpiEstim::estimate_R(
##                                      I,
##                                      #T.Start = start,
##                                      #T.End = end,
##                                      method = "uncertain_si",
##                                      config = EpiEstim::make_config(
##                                          list(
##                                              t_start = start,
##                                              t_end = end,
##                                              mean_si = mean(mean_estimates),
##                                              std_mean_si = sd_mean,
##                                              min_mean_si = 1,
##                                              max_mean_si = 8,
##                                              std_si = mean(sd_estimates),
##                                              std_std_si = sd_sd,
##                                              min_std_si = 1,
##                                              max_std_si = 5,
##                                              n1 = 100,
##                                              n2 = 100
##                                          )
##                                      )
##                                  )
##                 res$R$date <- df$date[end]
##                 res$R
##             },
##           .id = "time_window"
##         )

##     }, .id = "country"
## )


ggplot() +
    geom_line(
        data = r0s, aes(x = date, y = `Median(R)`, col = time_window)
    ) +
    geom_ribbon(
        data = r0s,
        aes(
            x = date,
            ymin = `Quantile.0.025(R)`,
            ymax = `Quantile.0.975(R)`,
            fill = time_window
            ),
        alpha = 0.3
    ) +
    facet_wrap(~country, ncol = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.position = "bottom")


nsim <- 1000
n_days <- 21

r0s <- split(r0s, r0s$country)
date_to_project_from <- "2020-03-03"
## Projections
projections <- purrr::imap(
    incidence,
    function(df, country) {
        country_r0 <- r0s[[country]]
        country_r0_tw <- split(
            country_r0, country_r0$time_window
        )
        ## Ger R0 around a week back and project for a period
        ## that also has some observations.
        country_r0_tw <- map(
            country_r0_tw, ~ dplyr::filter(., date == date_to_project_from)
        )
        I <- df$cases[df$date <= as.Date(date_to_project_from)]
        incid <- map_dfr(
            country_r0_tw,
            function(r0) {
                R <- r0[["Median(R)"]][nrow(r0)]
                pij <- matrix(1, nrow = 1, ncol = 1)
                incid <- project_sims(
                    nsim, I, R, si, pij, size, n_days, date_to_project_from
                )
                incid
        },
        .id = "time_window"
   )
  }
 )

quantiles <- purrr::map(
    projections,
    function(x) {
        qntls <- apply(
            x[ , -c(1, 2)], 1, quantile, probs = c(0.025, 0.5, 0.975)
        )
        out <- data.frame(
            date = x$date,
            time_window = x$time_window,
            q50  = qntls[2, ],
            q025 = qntls[1, ],
            q975 = qntls[3, ]
        )
        out
    }


)

for (country in names(incidence)) {

    p <- ggplot() +
        geom_line(
            data = incidence[[country]],
            aes(date, cases)
        ) +
        geom_line(
            data = quantiles[[country]],
            aes(date, q50, col = time_window)
        ) +
        geom_ribbon(
            data = quantiles[[country]],
            aes(x = date, ymin = q025, ymax = q975, fill = time_window),
            alpha = 0.3
        ) +
        theme_classic() +
        theme(legend.position = "bottom")
    plog <- p +
        scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        annotation_logticks(sides = "l")

    outfile <- snakecase::to_snake_case(country)
    ggsave(glue::glue("{outfile}.pdf"), p)
    ggsave(glue::glue("{outfile}_logscale.pdf"), plog)

}




####### rmse

rmse <- map2(
    incidence,
    projections,
    function(obs, pred) {
        x <- tail(obs$cases, n_days)
        tw_y <- split(pred, pred$time_window)
        rel_mse <- map_dfr(
            tw_y,
            function(pred2) {
                pred2 <- as.matrix(pred2[ , -c(1, 2)])
                out <- assessr::rel_mse(x, pred2)
                out <- data.frame(
                    date = tail(obs$date, n_days),
                    relmse = out
                )
                out$date <- as.Date(out$date)
                out
            },
            .id = "time_window"
        )
        rel_mse

    }
)


for (country in names(incidence)) {

    p <- ggplot() +
        geom_line(
            data = rmse[[country]],
            aes(date, relmse, col = time_window)
        ) +
        theme_classic() +
        theme(legend.position = "bottom") +
        geom_vline(
            xintercept = as.numeric(as.Date(date_to_project_from)),
            linetype = 4
        )
    plog <- p +
        scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        annotation_logticks(sides = "l")

    outfile <- snakecase::to_snake_case(country)
    ggsave(glue::glue("{outfile}_rmse.pdf"), p)
    ggsave(glue::glue("{outfile}_rmse_logscale.pdf"), plog)

}
