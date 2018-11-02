#' Select domains with a minimum number of contacts
#' 
#' Filters the domain list to retain intervals that have sufficiently high contacts recorded in N_CONTACTS
#' 
#' @param contacts The list of domains generated from count.contacts()
#' @param min Minimum number of fend pairs required in a domain
#' @param measure Element of the N_CONTACTS list to filter on
#' 
#' @return The contacts list, with additional information
#' 
#' @export
get_min_fend_pair_domains <- function(contacts, min=100, measure='HIGH_SCORE') {
	message(sprintf('[get_min_fend_pair_domains] Selecting domains with at least %d fend pairs in %s', min, measure))
	n.contacts <- as.numeric(lapply(contacts, function(x) x$N_CONTACTS[[measure]]))
	contacts[n.contacts>=min]
}
