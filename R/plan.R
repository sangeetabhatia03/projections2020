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

    rts = purrr::map(
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
                        ##si_distr = SI_Distr
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
                   out <- apply(x, 1, quantile, probs = c(0.025, 0.5, 0.975))
                   out <- as.data.frame(t(out))
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
    ),

    ## RMAE from median
    rmae = purrr::map_dfr(
        tables,
        function(df) {
            by_tw <- split(df, df$time_window)
            out <- purrr::map_dfr(
                by_tw,
                function(both) {

                    obs <- matrix(both$counts, ncol = 1)
                    pred <- matrix(both$`50%`, ncol = 1)
                    err <- assessr::rel_mae(obs, pred)

                    data.frame(
                        date = both$date,
                        rmae = err
                    )

                }, .id = "time_window"
            )
            out
        }, .id = "country"
        ),

    rmae_mean = dplyr::group_by(rmae, country, time_window) %>%
        summarise(mean_rmae = mean(rmae)) %>%
        ungroup(),
    ## For each counrty, pick the time window that minimises
    ## the error.

    best_tw = split(rmae_mean, rmae_mean$country) %>%
        purrr::map_dfr(function(x) x[which.min(x$mean_rmae), ]),

    ## And for the best time window, extract projections

    projections_best_tw = purrr::imap(
        qntls,
        function(qntl, cntry) {
            tw <- best_tw$time_window[best_tw$country == cntry]
            message("Best time window for ", cntry, " is ", tw)
            qntl[qntl$time_window == tw, ]
        }
     ),

    ## For flexible plotting

    incidence_df = purrr::map(
        incidence, as.data.frame
    ),

    rt_best_tw =  purrr::imap(
        rts,
        function(rt, cntry) {
            tw <- best_tw$time_window[best_tw$country == cntry]
            message("Best time window for ", cntry, " is ", tw)
            rt[[tw]]
        }
     ),

    plots = purrr::walk(
        countries_considered,
        function(cntry) {

            projections_best_tw[[cntry]]$date <- as.Date(
                projections_best_tw[[cntry]]$date
            )

            limits <- c(
                min(incidence_df[[cntry]]$dates),
                max(projections_best_tw[[cntry]]$date)
            )


            message("Country ", cntry)

            p <- ggplot() +
                geom_col(
                    data = incidence_df[[cntry]],
                    aes(dates, counts)
                )

            p <- p +
                geom_line(
                    data = projections_best_tw[[cntry]],
                    aes(date, `50%`)
                ) +
                geom_ribbon(
                    data = projections_best_tw[[cntry]],
                    aes(x = date, ymin = `2.5%`, ymax = `97.5%`),
                    alpha = 0.3
                ) + theme_classic() +
                scale_x_date(limits = limits) +
                xlab("") +
                ylab("Daily Incidence") +
                ggtitle(cntry)

            rt <- rt_best_tw[[cntry]]
            rt$date <- incidence_df[[cntry]]$dates[rt$t_start]
            rt$date <- as.Date(rt$date)
            rt_used <- tail(rt, 1)
            rt_used <- rt_used[rep(seq_len(nrow(rt_used)), n_days), ]
            rt_used$date <- seq(
                from = tail(rt$date, 1) + 1,
                length.out = n_days,
                by = "1 day"
            )
            rt <- rbind(rt, rt_used)

            p2 <- ggplot() +
                geom_line(
                    data = rt,
                    aes(date, `Median(R)`)
                ) +
                geom_ribbon(
                    data = rt,
                    aes(x = date,
                        ymin = `Quantile.0.025(R)`,
                        ymax = `Quantile.0.975(R)`,
                        ),
                    alpha = 0.3
                ) +
                scale_x_date(limits = limits) +
                theme_classic() +
                xlab("") +
                geom_hline(yintercept = 1, linetype = "dashed") +
                ylim(0, 5)
                ##expand_limits(y = 0)

            p3 <- cowplot::plot_grid(
                p, p2, nrow = 2, align="hv", rel_heights = c(2,1)
            )

            cowplot::save_plot(
                filename = glue::glue("{cntry}.png"),
                plot = p3
            )

        }

    )



)

