plan <- drake_plan(

    raw_data = readr::read_rds(file_in(infile)),

    deaths_to_use = raw_data[["D_active_transmission"]],

    ## Convert to incidence object
    tall_deaths = tidyr::gather(
        deaths_to_use, key = country, value = deaths, -dates
        ) %>%
        split(.$country) %>%
        purrr::map(~ ts_to_incid(
                ts = ., date_col = "dates", case_col = "deaths"
     )) %>%
    ## Only keep those which have at least 7 days of non-zero deaths
    ## Discuss with Pierre
    keep(
        function(x) (dim(x)[1] - first_nonzero_incidence(x)) >= 7
    ),

    si_distrs = purrr::map2(
        raw_data[["si_mean"]],
        raw_data[["si_std"]],
        function(mu, sigma) {
            reparams <- gamma_mucv2shapescale(
                mu = mu,
                cv = sigma / mu
            )
            miss_at_most <-0.001
            cutoff <- ceiling(
                qgamma(
                    1 - miss_at_most,
                    shape = reparams$shape,
                    scale = reparams$scale
                )
            )

            EpiEstim::discr_si(k = 0:cutoff, mu = mu, sigma = sigma)

        }
    ),
    ######### Vanilla EpiEstim with 7 day window in the first ########
    ######### instance.                                       ########
    rts = purrr::map(
        tall_deaths,
        function(df) {
            df <- subset(df, to = as.Date(date_to_project_from))
            out <- map(
                si_distrs,
                function(si_distr) {
                    ndays <- length(incidence::get_dates(df))
                    eps <- windows_endpoints(df)
                    config <- make_config(
                        t_start = eps$t_start,
                        t_end = eps$t_end,
                        si_distr = si_distr
                    )
                    res <- estimate_R(
                        df,
                        method = "non_parametric_si",
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
        tall_deaths,
        rt_samples,
        function(df, rt) {
            df <- subset(df, to = as.Date(date_to_project_from))
            out <- map2(
                rt, si_distrs,
                function(rt_si, si) {
                    projections::project(
                        x = df,
                        R = rt_si,
                        si = si,
                        n_sim = n_sim,
                        n_days = n_days,
                        R_fix_within = TRUE,
                        model = "negbin",
                        size = size
                    )

                }
            )
            out
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
               .id = "si_distr"
            )

        }
     ),

    plots = purrr::map(
        names(tall_deaths),
        function(z) {
            x <- tall_deaths[[z]]
            y <- qntls[[z]]
            y$date <- as.Date(y$date)
            y$si_distr <- factor(y$si_distr)
            x <- as.data.frame(x)

            p <- ggplot() +
                geom_col(data = x, aes(dates, counts)) +
                geom_line(data = y, aes(date, `50%`, col = si_distr)) +
                geom_ribbon(data = y,
                            aes(x = date,
                                ymin = `2.5%`,
                                ymax = `97.5%`,
                                fill = si_distr
                                ),
                            alpha = 0.3) +
                ggpubr::theme_pubr() +
                xlab("") +
                ylab("") + theme(legend.position = "none")
            z <- snakecase::to_snake_case(z)
            message(getwd())
            ggsave(
                filename = glue::glue("figures/{z}.png"),
                p
            )

        }
    ),
    ## Reformat, 7 columns, 10000 rows, 1 df per country
    ## can transpose trivially as these are matrices
    projections_t = purrr::map(
        projections, function(x) purrr::map(x, ~ t(.))
    ),

    ## Save in format required
    out = readr::write_rds(
        x = list(
            I_active_transmission = raw_data[["I_active_transmission"]],
            D_active_transmission = raw_data[["D_active_transmission"]],
            Country = raw_data[["Country"]],
            Rt_last = rt_samples,
            Predictons = projections_t
        ), path = file_out(outfile)
     ),

    ## apeestim
    inftvty = purrr::map(
        tall_deaths,
        function(deaths) {
            out <- purrr::map(
                si_distrs,
                function(si_distr) {
                    incid <- as.numeric(incidence::get_counts(deaths))
                    inftvty <- EpiEstim::overall_infectivity(
                        incid, si_distr
                        )
                    inftvty
                }
                )
            out
        }
    ),

    r_apeestim = purrr::map(
        names(tall_deaths),
        function(country) {
            deaths <- tall_deaths[[country]]
            r_prior <- c(1, 5)
            a <- 0.025
            trunctime <- first_nonzero_incidence(deaths)
            message("Truncating for ", country, " at ", trunctime)
            out <- purrr::map(
                si_distrs,
                function(si_distr) {
                    incid <- as.numeric(incidence::get_counts(deaths))
                    inftvty <- EpiEstim::overall_infectivity(
                        incid, si_distr
                    )
                    apeEstim(
                        incid,
                        si_distr,
                        inftvty,
                        r_prior,
                        a,
                        trunctime,
                        country
                   )
                }
             )
            out
        }
    )


)

