#' Determine the group(s) to which a domain belongs
#' 
#' A threshold of significance is applied to the calculated p-values, and a domain classification attributed.
#' 
#' @param contacts The list of domains generated from get_Q_values()
#' @param types The names of the elements to iterate over; these are produced by get.fend.samples() and are usually POTENTIAL (if get.fends was TRUE) and OBSERVED
#' @param measures The names of the p-value elements in Z_STATISTICS to iterate over; these are produced by get_Z_statistics()
#' @param thresholds A list of significance thresholds to apply
#' @param p.threshold The threshold to set for all p-value measures
#' @param z.quantile The quantile of Z scores above which are considered significant
#' @param z.threshold A threshold to apply, overriding z.quantile and ignoring `types`
#' 
#' @details
#' The thresholds list is a 2-level list, the first level indicates the null method used (eg 'OBSERVED') and the second level specifies the measure (eg 'Q' or 'Z').
#'
#' @return
#' The contacts list, with additional information in the Z_STATISTICS
#' 
#' @export
classify_domains.double <- function(contacts, types, measures=c('P','Q','p','Z'), thresholds, p.thresholds=c(1e-4,1e-8), z.quantiles=c(0.5,0.75), z.thresholds=NULL) {
	message('[classify_domains] Defining parameters to classify domains:')
	if(missing(types))
		types <- unique(sapply(contacts, function(x) names(x$Z_STATISTICS)))

	if(missing(thresholds))
		thresholds <- list()

	if(class(thresholds)!='list')
		thresholds <- list()

	if(missing(z.thresholds))
		z.thresholds <- NULL

	p.measures <- measures[! measures %in% 'Z']
	thresholds <- lapply(SelfName(types), function(type) {
		if(is.null(z.thresholds))
			z.threshold <- get_Z_quantiles(contacts, quantiles=z.quantiles, measures=type)[[type]]
		append(MergeLists(thresholds[[type]], setNames(rep(list(p.thresholds),length(p.measures)), p.measures)),
		       list(Z=z.thresholds))})
	print(thresholds)
	print(plyr::ldply(lapply(thresholds, function(x) data.frame(Threshold=x, row.names=NULL)), .id='Null method'))

	message('[classify_domains] Classifying domains')
	progress_bar <- txtProgressBar(min=0, max=length(contacts), initial=0, style=3) 

	for(i in seq(contacts)) {
		contacts[[i]]$CLASSIFICATIONS$Z_STATISTICS <- lapply(SelfName(types), function(type) {
			lapply(SelfName(measures), function(measure) {
				classes <- names(contacts[[i]]$Z_STATISTICS[[type]][[measure]])
				significant <- rep('AMBIGUOUS', length(classes))

				if(measure=='Z') {
					sig <- contacts[[i]]$Z_STATISTICS[[type]][[measure]]>=thresholds[[type]][[measure]][2]
					nonsig <- contacts[[i]]$Z_STATISTICS[[type]][[measure]]<thresholds[[type]][[measure]][1]
				} else {
					sig <- contacts[[i]]$Z_STATISTICS[[type]][[measure]]<=thresholds[[type]][[measure]][2]
					nonsig <- contacts[[i]]$Z_STATISTICS[[type]][[measure]]>thresholds[[type]][[measure]][1]
				}

				return_list <- list(SECTORS=list(SIGN=classes[sig],
				                                 NSIG=classes[nonsig],
				                                 AMBI=classes[!sig & !nonsig]))

				classes <- return_list$SECTORS$SIGN
				if(!all(sig | nonsig))
					classes <- 'AMBIGUOUS'
		 		if(all(nonsig))
		 			classes <- 'NONE'

		 		return_list$CLASSIFICATION <- paste(classes, collapse='_')
		 		return_list})})
		setTxtProgressBar(progress_bar, i)
	}
	message()

	contacts
}

classify_domains.single <- function(contacts, types, measures=c('P','Q','p','Z'), thresholds, p.threshold=1e-5, z.quantile=0.75, z.threshold=NULL) {
	message('[classify_domains] Defining parameters to classify domains:')
	if(missing(types))
		types <- unique(sapply(contacts, function(x) names(x$Z_STATISTICS)))

	if(missing(thresholds))
		thresholds <- list()

	if(class(thresholds)!='list')
		thresholds <- list()

	if(missing(z.threshold))
		z.threshold <- NULL

	p.measures <- measures[! measures %in% 'Z']
	thresholds <- lapply(SelfName(types), function(type) {
		if(is.null(z.threshold))
			z.threshold <- get_Z_quantiles(contacts, quantiles=z.quantile[1], measures=type)[[type]]
		append(MergeLists(thresholds[[type]], as.list(setNames(rep(p.threshold,length(p.measures)), p.measures))),
		       list(Z=z.threshold))})
	
	print(plyr::ldply(lapply(thresholds, function(x) data.frame(Threshold=x, row.names=NULL)), .id='Null method'))

	message('[classify_domains] Classifying domains')
	progress_bar <- txtProgressBar(min=0, max=length(contacts), initial=0, style=3) 

	for(i in seq(contacts)) {
		contacts[[i]]$CLASSIFICATIONS$Z_STATISTICS <- lapply(SelfName(types), function(type) {
			lapply(SelfName(measures), function(measure) {
				if(measure=='Z') {
					classes <- names(contacts[[i]]$Z_STATISTICS[[type]][[measure]])[contacts[[i]]$Z_STATISTICS[[type]][[measure]]>=thresholds[[type]][[measure]]]
				} else {
					classes <- names(contacts[[i]]$Z_STATISTICS[[type]][[measure]])[contacts[[i]]$Z_STATISTICS[[type]][[measure]]<=thresholds[[type]][[measure]]]
				}

				return_list <- list(SECTORS=classes)

				if(length(classes)==0)
					classes <- 'NONE'
		 		
		 		# paste(classes, collapse='_')
		 		return_list$CLASSIFICATION <- paste(classes, collapse='_')
		 		return_list})})
		setTxtProgressBar(progress_bar, i)
	}
	message()

	contacts
}

classify_domains.bak <- function(contacts, measures=c('P','Q','p','Z'), p.threshold=1e-5, z.threshold=20) {
	message('[classify_domains] Classifying domains')

	progress_bar <- txtProgressBar(min=0, max=length(contacts), initial=0, style=3) 

	for(i in seq(contacts)) {
		contacts[[i]]$CLASSIFICATIONS$Z_STATISTICS <- lapply(contacts[[i]]$Z_STATISTICS, function(type_list) {
			lapply(SelfName(measures), function(measure) {
				classes <- names(type_list[[measure]])[type_list[[measure]]<=p.threshold]
				if(measure=='Z')
					classes <- names(type_list[[measure]])[type_list[[measure]]>=z.threshold]
				if(length(classes)==0)
					classes <- 'NONE'
		 		paste(classes, collapse='_')})})
		setTxtProgressBar(progress_bar, i)
	}
	message()

	contacts
}

get_Z_quantiles <- function(contacts, quantiles=0.75, measures) {
	if(missing(measures))
		measures <- unique(sapply(contacts, function(x) names(x$Z_STATISTICS)))
	lapply(SelfName(measures), function(M) quantile(unlist(plyr::ldply(contacts, function(x) x$Z_STATISTICS[[M]]$Z)[-1]), probs=quantiles))
}


#' Determine the group(s) to which a domain belongs
#' 
#' A data.frame of domains and the calls of significant sectors is translated into the set of sectors that are enriched.
#' 
#' @param domain_calls A data.frame with domains as rows and sectors as columns.
#' 
#' @return
#' Character vector of length nrow(domain_calls) with `_` separated enriched sectors.
classify_domains.bak <- function(domain_calls) {
	classes <- plyr::alply(domain_calls, .margins=1, function(x) names(x)[x])
	names(classes) <- rownames(domain_calls)
	classes <- llply(classes, function(x) {if(length(x)==0) {x<-'NONE'}; x})
	groups <- plyr::ldply(classes, paste, collapse='_')
	groups
}
