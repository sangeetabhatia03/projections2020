inurl <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv"

countries_considered <- c(
    "Italy", "Japan", "Iran",
    "France", "United Kingdom", "Spain"
)
mean_si <- 3.96
std_si <- 4.75
size <- 0.16
time_window <- c(2, 4, 6, 8)
names(time_window) <- time_window
date_to_project_from <- "2020-03-08"
n_sim <- 1000
n_days <- 10

SItrunc <- 20

SI_Distr <- sapply(
  0:SItrunc,
  function(e) EpiEstim::discr_si(e, mean_si, std_si)
)

SI_Distr <- SI_Distr / sum(SI_Distr)
