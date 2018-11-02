#' Get the values for sectors of each domain
#' 
#' Generates a set of values that could be used to describe the distribution of fend pairs among the 4 domain sectors
#' 
#' @param contacts The list of domains generated from get.contacts()
#' 
#' @return The contacts list, with additional information
#' 
#' @export
get.domain.position.classifications <- function(contacts) {
	message('Getting distribution of fend positions')
	if(!is.factor(contacts[[1]]$CONTACTS[[1]]$POSITION))
		stop('The POSITION of CONTACTS is not a factor!')

	progress_bar <- txtProgressBar(min=0, max=length(contacts), initial=0, style=3) 
	for(i in 1:length(contacts)) {
		domain <- contacts[[i]]

		domain$FEND.SCORE.DISTRIBUTION <- lapply(domain$CONTACTS, function(x) table(x$POSITION))

		domain$FEND.SCORE.DISTRIBUTION$SCORES_dALL_mHS <- domain$FEND.SCORE.DISTRIBUTION$HIGH_SCORE/domain$FEND.SCORE.DISTRIBUTION$ALL*domain$N_CONTACTS$HIGH_SCORE
		domain$FEND.SCORE.DISTRIBUTION$SCORES_dHS_mALL <- domain$FEND.SCORE.DISTRIBUTION$HIGH_SCORE/domain$N_CONTACTS$HIGH_SCORE*domain$FEND.SCORE.DISTRIBUTION$ALL
		domain$FEND.SCORE.DISTRIBUTION$SCORES_dALL <- domain$FEND.SCORE.DISTRIBUTION$HIGH_SCORE/domain$FEND.SCORE.DISTRIBUTION$ALL
		domain$FEND.SCORE.DISTRIBUTION$SCORES_dHS <- domain$FEND.SCORE.DISTRIBUTION$HIGH_SCORE/domain$N_CONTACTS$HIGH_SCORE

		if(!is.null(domain$FENDS.2D)) {
			# domain$FENDS.MAX <- lapply(domain$FENDS.2D, function(x) table(factor(x$POSITION, levels=c(CORNER=0,FORWARD=0,REVERSE=0,OTHER=0))))
			domain$FENDS.MAX <- table(domain$FENDS.2D$POSITION)
			domain$FEND.SCORE.DISTRIBUTION$SCORES_dMAX <- domain$FEND.SCORE.DISTRIBUTION$HIGH_SCORE/domain$FENDS.MAX
		}

		domain$FEND.SCORE.DISTRIBUTION$SCORES <- domain$FEND.SCORE.DISTRIBUTION$SCORES_dHS

		domain$FEND.SCORE.DISTRIBUTION <- lapply(domain$FEND.SCORE.DISTRIBUTION, function(x) {x[is.na(x)] <- 0; x})

		contacts[[i]] <- domain
		setTxtProgressBar(progress_bar, i)
	}
	message()
	domain <- NULL
	gc()

	contacts
}

get.domain.position.classifications.bak <- function(contacts) {
	message('Getting distribution of fend positions')
	if(!is.factor(contacts[[1]]$CONTACTS[[1]]$POSITION))
		stop('The POSITION of CONTACTS is not a factor!')
	llply(contacts, function(domain) {

		domain$FEND.SCORE.DISTRIBUTION <- lapply(domain$CONTACTS, function(x) table(x$POSITION))

		# domain$FEND.SCORE.DISTRIBUTION <- lapply(domain$CONTACTS, function(x) table(factor(x$POSITION, levels=c('CORNER','FORWARD','REVERSE','OTHER'))))

		# domain$FEND.SCORE.DISTRIBUTION <- lapply(domain$CONTACTS, function(x) {
		# 	full.table <- setNames(rep(0,4), c('CORNER','FORWARD','REVERSE','OTHER'))
		# 	pos.table <- table(x$POSITION)
		# 	full.table[names(pos.table)] <- pos.table
		# 	full.table})

		domain$FEND.SCORE.DISTRIBUTION$SCORES_dALL_mHS <- domain$FEND.SCORE.DISTRIBUTION$HIGH_SCORE/domain$FEND.SCORE.DISTRIBUTION$ALL*domain$N_CONTACTS$HIGH_SCORE
		domain$FEND.SCORE.DISTRIBUTION$SCORES_dHS_mALL <- domain$FEND.SCORE.DISTRIBUTION$HIGH_SCORE/domain$N_CONTACTS$HIGH_SCORE*domain$FEND.SCORE.DISTRIBUTION$ALL
		domain$FEND.SCORE.DISTRIBUTION$SCORES_dALL <- domain$FEND.SCORE.DISTRIBUTION$HIGH_SCORE/domain$FEND.SCORE.DISTRIBUTION$ALL
		domain$FEND.SCORE.DISTRIBUTION$SCORES_dHS <- domain$FEND.SCORE.DISTRIBUTION$HIGH_SCORE/domain$N_CONTACTS$HIGH_SCORE

		if(!is.null(domain$FENDS.2D)) {
			domain$FENDS.MAX <- lapply(domain$FENDS.2D, function(x) table(factor(x$POSITION, levels=c(CORNER=0,FORWARD=0,REVERSE=0,OTHER=0))))
			# domain$FENDS.MAX <- lapply(domain$FENDS.2D, function(x) {
			# 	full.table <- c(CORNER=0,FORWARD=0,REVERSE=0,OTHER=0)
			# 	pos.table <- table(x$POSITION)
			# 	full.table[names(pos.table)] <- pos.table
			# 	full.table})
			domain$FEND.SCORE.DISTRIBUTION$SCORES_dMAXall <- domain$FEND.SCORE.DISTRIBUTION$HIGH_SCORE/domain$FENDS.MAX$ALL
			domain$FEND.SCORE.DISTRIBUTION$SCORES_dMAXvalid <- domain$FEND.SCORE.DISTRIBUTION$HIGH_SCORE/domain$FENDS.MAX$VALID
		}

		domain$FEND.SCORE.DISTRIBUTION$SCORES <- domain$FEND.SCORE.DISTRIBUTION$SCORES_divMAXvalid

		domain}, .progress='text')
}

