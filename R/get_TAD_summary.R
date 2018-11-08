#' @export
get_TAD_summary <- function(contacts, min_n_fends=100, z.thresholds=c(5,15)) {
	if(is.element('character', class(contacts)))
		load(contacts)

	contacts <- domainClassifyR::get_min_fend_pair_domains(contacts, min=min_n_fends)
	contacts <- domainClassifyR::get_Q_values(contacts)

	classified_domains <- domainClassifyR::classify_domains.double(contacts, z.thresholds=z.thresholds)
	z_groups <- sapply(classified_domains, function(x) x$CLASSIFICATIONS$Z_STATISTICS$OBSERVED$Z$CLASSIFICATION)

	domain_intervals <- classified_domains %>%
	                    ldply(function(x) mishaHelpR::convert_intervals_2D_to_1D(x$DOMAIN, set=1, extra_cols=c('SIZE1','intervalID')), .id=NULL) %>%
	                    # [c('chrom','start','end','intervalID','SIZE1')] %>%
	                    mutate(Z_GROUP=z_groups, COMPARTMENT=0)

	domain_intervals <- domain_intervals %>%
	                    dplyr::bind_rows(., mutate(., Z_GROUP='ALL'))

}
