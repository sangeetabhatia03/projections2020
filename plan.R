plan <- drake_plan(

    raw_data = readr::read_csv(inurl) %>%
        dplyr::filter(`Country/Region` != "Mainland China"),

    tall_data = tidyr::gather(
        raw_data, key = date, value = cum_cases, -(`Province/State`:Long)
    ) %>%
        mutate_at(vars("date"), lubridate::mdy) %>%
        fix_manually() %>%
        dplyr::count(
        `Country/Region`, date, wt = cum_cases, name = "cum_cases"
    ),

    incidence = tall_data[tall_data$`Country/Region`
                          %in%
                          countries_considered, ] %>%
        split(.$`Country/Region`) %>%
        purrr::map(function(x) {
            x$cases <- c(0, diff(x$cum_cases))
            if (any(x$cases < 0)) message(x$`Country/Region`[1])
            x <- tidyr::uncount(data = x, weight = cases)
            x <- incidence(x$date)
            x
        }
    ),

    rts = map(
        incidence,
        function(df) {
            df <- subset(df, to = as.Date(date_to_project_from))
            out <- map(
                time_window,
                function(tw) {
                    ndays <- length(incidence::get_dates(df))
                    eps <- windows_endpoints(ndays = ndays, tw = tw)
                    config <- make_config(
                        t_start = eps$t_start,
                        t_end = eps$t_end,
                        mean_si = mean_si,
                        std_si = std_si
                    )
                    res <- estimate_R(
                        df,
                        method = "parametric_si",
                        config = config
                    )
                    res$R
                }
            )

        }
    ),

    ## Samples from posterior distribution of R0

    rt_samples = map(
        rts,
        function(x) map(x, ~ r_samples(.))
    ),

    ## Projections using projections package with
    ## negative-binomial offspring distribution

    projections = map2(
        incidence,
        rt_samples,
        function(df, rt) {
            df <- subset(df, to = as.Date(date_to_project_from))
            out <- map(
                rt,
                function(rt_tw) {
                    project(
                        x = df,
                        R = rt_tw,
                        si = SI_Distr,
                        n_sim = n_sim,
                        n_days = n_days,
                        R_fix_within = TRUE,
                        model = "negbin",
                        size = size
                    )

                }
            )
        }
    )

    ## RMSE


)
