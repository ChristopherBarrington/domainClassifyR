#' Counts the number of extracted fend pairs
#' 
#' For every element of the CONTACTS list, count the number of fend pairs
#' 
#' @param contacts The list of domains generated from get.contacts()
#' 
#' @return The contacts list, with additional information
#' 
#' @export
count.contacts <- function(contacts) {
	message('[count.contacts] Counting contacts')
	progress_bar <- txtProgressBar(min=0, max=length(contacts), initial=0, style=3) 
	for(i in 1:length(contacts)) {
		if(!is.null(contacts[[i]]$CONTACTS))
			contacts[[i]]$N_CONTACTS <- lapply(contacts[[i]]$CONTACTS, nrow)
		
		if(!is.null(contacts[[i]]$FENDS))
			contacts[[i]]$N_FENDS <- nrow(contacts[[i]]$FENDS)
		
		if(!is.null(contacts[[i]]$FENDS.2D))
			contacts[[i]]$N_FENDS.2D <- nrow(contacts[[i]]$FENDS.2D)
		
		setTxtProgressBar(progress_bar, i)
	}
	message()

	contacts
}

count.contacts.bak <- function(contacts) {
	message('[count.contacts] Counting contacts')
	llply(contacts, function(domain) {
		if(!is.null(domain$CONTACTS))
			domain$N_CONTACTS <- lapply(domain$CONTACTS, nrow)
		if(!is.null(domain$FENDS))
			domain$N_FENDS <- lapply(domain$FENDS, nrow)
		if(!is.null(domain$FENDS.2D))
			domain$N_FENDS.2D <- lapply(domain$FENDS.2D, nrow)
		domain}, .progress='text')
}
