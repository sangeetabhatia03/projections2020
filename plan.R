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
    ),

    ## Central estimates for each country for each
    ## time window
    qntls = purrr::map(
        projections,
        function(country) {
           purrr::map_dfr(
               country,
               function(x) {
                   out <- apply(x, 1, quantile, probs = 0.5)
                   out <- as.data.frame(out)
                   out <- tibble::rownames_to_column(out, var = "date")
                   out
               },
               .id = "time_window"
            )

        }
     ),

    tables = purrr::map2(
        incidence,
        qntls,
        function(obs, pred) {

            pred$date <- as.Date(pred$date)
            obs <- as.data.frame(obs)
            x <- dplyr::full_join(
                x = obs, y = pred, by = c("dates" = "date")
                )
            x <- na.omit(x)
            x <- dplyr::arrange(x, time_window)

        }
    )

    ## RMSE from median
    rmse = purrr::map_dfr(
        tables,
        function(df) {
            by_tw <- split(df, df$time_window)
            out <- purrr::map_dfr(
                by_tw,
                function(both) {

                    obs <- matrix(both$counts, ncol = 1)
                    pred <- matrix(both$out, ncol = 1)
                    err <- assessr::rel_mse(obs, pred)

                    data.frame(
                        date = df$date,
                        rel_rmse = err
                    )

                }, .id = "time_window"
            )
            out
        }, .id = "country"
    )

)

