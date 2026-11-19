# Resample loci and re-estimate the essentiality cutoff used in Figure 2A.
# Fitting follows Bio-Tradis bin/tradis_essentiality.R (sanger-pathogens/Bio-Tradis).
# The published diagnostic evaluated its tail-odds score on 5,000 index values.
# Usage: Rscript --vanilla figure_2a_bootstrap.R sites.tsv [replicates] [seed]
#
# The input is the original 6,206-locus TraDIS gene-insertion export with an
# ins_index column. The percentile interval measures sensitivity to sampling
# loci from this pooled library, conditional on treating loci as exchangeable.
# It does not measure biological replicate variability or per-gene significance.

cutoff_fit <- function(index) {
    # Preserve the original indexing used to choose the loess fitting window.
    initial_hist <- hist(index, breaks = 200, plot = FALSE)
    density <- initial_hist$density
    peak <- which.max(density[10:length(density)])
    last_bin <- floor(initial_hist$mids[peak + 3L] * 2000)
    detail <- hist(index[index < last_bin / 2000],
                   breaks = (0:last_bin) / 2000, plot = FALSE)
    position <- seq_along(detail$density)
    smooth <- loess(detail$density ~ position)
    division <- detail$mids[which.min(predict(smooth))]

    # Fit the exponential to low indices, including zeros, and the gamma to
    # indices in the high mode up to the Bio-Tradis tail-truncation rule.
    full_hist <- hist(index, breaks = "FD", plot = FALSE)
    last_dense <- max(which(full_hist$counts > 5))
    upper <- full_hist$mids[last_dense]
    low <- index[index < division]
    high <- index[index >= division & index < upper]
    if (!length(low) || !length(high) || any(high <= 0)) {
        stop("The component fitting samples are empty or invalid.")
    }
    exp_rate <- 1 / mean(low)
    high_mean <- mean(high)
    # Gamma maximum likelihood, with rate = shape / mean. Profiling the rate
    # reduces this to the same two-parameter fit as MASS::fitdistr.
    dispersion <- log(high_mean) - mean(log(high))
    shape <- uniroot(function(k) log(k) - digamma(k) - dispersion,
                     c(0.001, 100000), tol = 1e-9)$root
    rate <- shape / high_mean
    score <- function(x) {
        log_odds_gamma <- pgamma(x, shape, rate, log.p = TRUE) -
            pgamma(x, shape, rate, lower.tail = FALSE, log.p = TRUE)
        log_odds_exp <- pexp(x, exp_rate, log.p = TRUE) -
            pexp(x, exp_rate, lower.tail = FALSE, log.p = TRUE)
        (log_odds_gamma + log_odds_exp) / log(2)
    }
    max_index <- max(index)
    step <- max_index / 4999
    continuous <- uniroot(function(x) score(x) + 2,
                          c(step / 1000, max_index), tol = 1e-11)$root
    ambiguous <- uniroot(function(x) score(x) - 2,
                         c(step / 1000, max_index), tol = 1e-11)$root
    # The original plot adds zero-index loci a second time in the low-mode
    # scaling factor. These density multipliers should not be used as priors.
    c(essential = floor(continuous / step) * step,
      ambiguous = ceiling(ambiguous / step) * step,
      exp_rate = exp_rate, gamma_shape = shape, gamma_rate = rate,
      exp_weight = (length(low) + sum(index == 0)) / length(index),
      gamma_weight = length(high) / length(index))
}

main <- function(args) {
    if (!length(args) || length(args) > 3L) {
        stop("Usage: Rscript figure_2a_bootstrap.R sites.tsv [replicates] [seed]")
    }
    replicates <- if (length(args) >= 2) as.integer(args[2]) else 20000L
    seed <- if (length(args) >= 3) as.integer(args[3]) else 9446L
    if (!is.finite(replicates) || replicates < 2L || !is.finite(seed)) {
        stop("Supply an integer replicate count of at least two and an integer seed.")
    }
    header <- readLines(args[1], n = 1, warn = FALSE)
    sep <- if (grepl("\t", header, fixed = TRUE)) "\t" else ","
    sites <- read.table(args[1], header = TRUE, sep = sep, quote = "\"", comment.char = "")
    index <- sites$ins_index
    if (length(index) != 6206L || any(!is.finite(index)) || any(index < 0)) {
        stop("Expected the original 6,206 nonnegative, finite insertion indices.")
    }
    fitted <- cutoff_fit(index)
    if (round(fitted["essential"], 5) != 0.03209 ||
        round(fitted["ambiguous"], 5) != 0.03975) {
        stop("Refit did not reproduce both published model thresholds.")
    }
    set.seed(seed)
    boot <- replicate(replicates, {
        draw <- sample.int(length(index), size = length(index), replace = TRUE)
        cutoff_fit(index[draw])[["essential"]]
    })
    limits <- quantile(boot, c(0.025, 0.975, 0.005, 0.995), names = FALSE, type = 7)
    # In the original plot the 90 bins are centered at steps from 0 to max.
    # Edges = ((0:90) - 0.5) * histogram_bin_width.
    histogram_bins <- 90L
    result <- data.frame(replicates, successful = length(boot), seed, loci = length(index),
                         fitted_cutoff = unname(fitted["essential"]),
                         fitted_ambiguous = unname(fitted["ambiguous"]),
                         lower_95 = limits[1], upper_95 = limits[2],
                         exp_rate = unname(fitted["exp_rate"]),
                         gamma_shape = unname(fitted["gamma_shape"]),
                         gamma_rate = unname(fitted["gamma_rate"]),
                         exp_weight = unname(fitted["exp_weight"]),
                         gamma_weight = unname(fitted["gamma_weight"]),
                         histogram_bins,
                         histogram_bin_width = max(index) / (histogram_bins - 1),
                         model_grid_points = 5000L,
                         lower_99 = limits[3], upper_99 = limits[4])
    write.csv(result, stdout(), row.names = FALSE)
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
