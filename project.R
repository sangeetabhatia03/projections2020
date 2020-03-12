##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @param pij : row vector of length n where n is the
##' number of locations. Each entry is the probability of
##' flow from i to j.
##' @param r_t : t X n matrix. Entry ij is the reproduction
##' number at location j and time i.
##' @param i_t : t X n matrix. Entry ij is the incidence
##' count at location j and time i.
##' @param w_t : column vector of length t. Samples from the
##' serial distribution.
##' @title
##' @return the mean of the Poisson distribution determining
##' the distribution of incidences at location j at time t.
##' @author Sangeeta Bhatia
##' @export
lambda_j_t <- function(pij, r_t, incid, ws){
    out <- (ws  %*% incid) * r_t
    out <-  out %*% pij
    out
}

##' Project future incidence at a location
##' Simulates future incidence based on past incidence data at a
##' given set of locations, the Reproduction numbers and
##' the distribution of the serial interval.
##' Runs a single simulation.
##' @export
##' @param incid A data frame of incidence counts for n locations and t days.
##' Dates run down the rows and locations are across columns.
##' @param R either a 1 X n matrix of reproduction numbers
##' for each location or a matrix of n_days X n with R for
##' each time instant and each location for which projection
##' is to be made.
##' @param si A vector drawn from serial interval
##' distribution. If the length of the vector is less than
##' the nrow(incid) + n_days, it is padded with 0s.
##' @param pij a n x n matrix of probabilities where n is
##' the number of locations.
##' @param n_days number of days for which a simulation
##' should be run.
##' @return data frame containing projected incidence count
##' with n_days rows and n columns.
##' @author Sangeeta Bhatia
##' @export
project <- function(incid, R, si, pij, size, n_days = 7){

    if (ncol(R) != ncol(incid)){
        stop("R should be a either a 1 X N or T X N matrix.")
    }
    if (nrow(R) != nrow(incid)){
        if (nrow(R) != 1){
            stop("R should be a either a 1 X N or T X N
                  matrix.")
        }
    }

    n_loc <- ncol(incid)
    out   <- matrix(0, nrow = n_days, ncol = n_loc)
    out   <- rbind(incid, out)
    start <- nrow(incid) + 1
    end   <- nrow(out)

    if (nrow(R) == 1){
        R <- matrix(R,
                    nrow = end - start + 1,
                    ncol = n_loc, byrow = TRUE)
    }
    if (length(si) < end)
         ws <- rev(c(si, rep(0, end - length(si))))
    else ws <- rev(si)
    for (i in start:end){
        i_t <- out[1:i, ]
        w_t <- utils::tail(ws, i)
        r_t <- R[i - start + 1, ]
        mu  <- lambda_j_t(pij, r_t, i_t, w_t)
        ## not sure why, check project function in package
        ## projections
        size <- (w_t %*% i_t) * size
        if (size == 0) size <- 1
        out[i, ] <- stats::rnbinom(1, size = size, mu = mu)
    }

    out[start:nrow(out), ]
}

custom_join <- function(sim1, sim2) {
    dplyr::left_join(sim1, sim2, by = "date")
}


project_sims <- function(nsim, I, R, si, pij, size, n_days, date_to_project_from) {
    incid <- map(
        1:nsim,
        function(i) {
            incid <- project(
                incid = matrix(I, ncol = 1),
                R = matrix(R, ncol = 1),
                si = SI_Distr,
                pij = matrix(1, nrow = 1, ncol = 1),
                size = size,
                n_days = n_days
            )
            out <- data.frame(
                date = seq(
                    from = as.Date(date_to_project_from),
                    length.out = n_days,
                    by = "1 day"
                ),
                sim = incid
            )
            out
        }

    )
    incid <- Reduce(f = custom_join, x = incid)
    ## colnames(incid)[3:length(incid)] <- paste0(
    ##     "sim", length(3:length(incid))
    ## )
    incid

}
