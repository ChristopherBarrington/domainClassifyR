#' Subset high-scoring fend pairs from all observed contacts
#' 
#' Select the fend pairs that are above a score threshold and save their information into the HIGH_SCORE slot of CONTACTS for each domain
#' 
#' @param contacts The list of domains generated from get.contacts() and add.domain.features()
#' @param min_score The minimum score that a fend pair must have to be retained
#' 
#' @return The contacts list, with additional information
#' 
#' @export
get.high.score.fends <- function(contacts, min_score=50) {
	message(sprintf('[get.high.score.fends] Identifying high-scoring fends [>=%d]', min_score))
	# llply(contacts, function(x) {x$CONTACTS$HIGH_SCORE <- x$CONTACTS$ALL[x$CONTACTS$ALL$SCORE>=min_score,]; x}, .progress='text')

	progress_bar <- txtProgressBar(min=0, max=length(contacts), initial=0, style=3) 
	for(i in 1:length(contacts)) {
		contacts[[i]]$CONTACTS$HIGH_SCORE <- contacts[[i]]$CONTACTS$ALL[contacts[[i]]$CONTACTS$ALL$SCORE>=min_score,]
		setTxtProgressBar(progress_bar, i)
	}
	message()

	contacts
}
