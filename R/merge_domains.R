#' Combine a list of domain contacts entries
#' 
#' Concatenates the fend pairs for ALL and HIGH_SCORE from a list of domains and
#' defines the interval as the region covered by the merged domains
#' 
#' @param contacts List of domain contacts to be merged
#' @param intervalID Character string that will be used in the DOMAIN entry
#' 
#' @return A list for a new meta-domain
#'
#' @examples
#' contacts <- append(contacts, merge_domains(contacts=contacts[1:3], intervalID='foo+bar+baz'))
#' print.fixed(plot_domain_contacts(contacts, ncol=3))
#'
#' @export
merge_domains <- function(contacts, intervalID='merged') {
	domains <- ldply(contacts, function(x) x$DOMAIN)[-1]
	output_list <- list(list(CONTACTS=list(ALL=ldply(contacts, function(x) x$CONTACTS$ALL)[-1],
	                                       HIGH_SCORE=ldply(contacts, function(x) x$CONTACTS$HIGH_SCORE)[-1]),
	                         DOMAIN=data.frame(chrom1=domains$chrom1[1],
	                                           start1=min(domains$start1),
	                                           end1=max(domains$end2),
	                                           chrom2=domains$chrom2[1],
	                                           start2=min(domains$start1),
	                                           end2=max(domains$end2),
	                                           intervalID=intervalID)))
	names(output_list) <- intervalID
	output_list
}
