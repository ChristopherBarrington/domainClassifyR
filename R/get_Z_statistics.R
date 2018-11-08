#' Calculates Z statistics and p-values
#' 
#' Calculates the probability that a domain sector contains more high-score fend pairs than expected by chance. Uses the
#' observed data and theoretical fend pairs as expected proportions.
#' 
#' @param contacts List of domain entries
#' 
#' @details Converts the Z-statistic to a p-value of enrichment using pnorm(-Z). For comparison, a p-value is
#' calculated using pnorm(OBS, mean=MEAN, sd=SDEV, lower.tail=FALSE). No correction for multiple testing is made.
#' 
#' @return The contacts list, with additional information. Contains `P`, derived from the Z-statistic and `p` derived from the observed high-score count.
#'
#' @export
get_Z_statistics <- function(contacts) {
	message('[get_Z_statistics] Calculating Z-statistics and p-values')
	progress_bar <- txtProgressBar(min=0, max=length(contacts), initial=0, style=3) 

	for(i in 1:length(contacts)) {
		contacts[[i]]$Z_STATISTICS <- lapply(contacts[[i]]$CONTACTS$RANDOM_SAMPLE_TABLES, function(x) {
			z <- list(MEAN=plyr::aaply(x, 2, mean),
	    	          SDEV=plyr::aaply(x, 2, sd),
	    	          OBS=contacts[[i]]$FEND.SCORE.DISTRIBUTION$HIGH_SCORE)
			z$Z <- (z$OBS-z$MEAN)/z$SDEV
			z$Z[is.na(z$Z)] <- 0
			z$P <- pnorm(-z$Z)

			z$p <- pnorm(z$OBS, mean=z$MEAN, sd=z$SDEV, lower.tail=FALSE)
			z$p[is.na(z$p)] <- 1

			z})

		setTxtProgressBar(progress_bar, i)
	}
	message()

	contacts
}

