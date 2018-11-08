#' Correct for multiple testing in a data.frame
#' 
#' Uses p.adjust to correct for multiple testing on all p-values in the data.frame.
#' 
#' @param contacts The list of domains generated from get.contacts()
#' @param types Correct the p-values calculated from which null distribution(s)? `POTENTIAL` uses all theoretically possible fend pairs or `OBSERVED` uses the fend pairs that were identified.
#' @param method The correction method to apply, see p.adjust()
#' 
#' @return The contacts list, with additional information in the Z_STATISTICS.
#' 
#' @export
get_Q_values <- function(contacts, types=c('POTENTIAL','OBSERVED'), method='BH') {
	message('[get_Q_values] Correcting for multiple testing')
	for(type in types) {
		x <- ldply(contacts, function(x) x$Z_STATISTICS[[type]]$P)[-1]

		if(nrow(x)==length(contacts)) {
			adj <- p.adjust(unlist(x), method=method)
			dim(adj) <- dim(x)
			colnames(adj) <- colnames(x)

			for(i in seq(length(contacts))) {
				contacts[[i]]$Z_STATISTICS[[type]]$Q <- adj[i,]
			}
		}
	}
	contacts
}
