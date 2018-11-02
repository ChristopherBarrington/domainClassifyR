#' Sample fragment ends randomly to get a null distribution
#' 
#' Randomly samples the total observed fend pairs and potential fend pairs to randomly select the same number of fend pairs that are
#' high-scoring. The random selection is repeated and a null distribution for domain sector scores determined.
#' 
#' @param contacts The list of domains generated from get.contacts()
#' @param n Number of times random selection is repeated
#' @param cluster_size Number of parallel threads to run
#' 
#' @details Parallelisation is achieved using the registerDoMC package and the cluster_size option.
#' 
#' @return The contacts list, with additional information
#' 
#' @export
get.fend.samples <- function(contacts, random_seed=FALSE, n=1000, cluster_size=20, distribute_method='SGE', parse_result=TRUE, sleeptime=300, sge_resources=list(), save.command=TRUE) {
	message(sprintf('[get.fend.samples] Randomly selecting contacts from ALL [x %s]', format(n, big.mark=',', scientific=FALSE)))

	# commands <- sprintf('domainClassifyR:::get.fend.samples.function(domain=contacts[[%d]], random_seed=%s, n=%d)', seq(length(contacts)), random_seed, n)
	# contacts <- distributR::distributr(commands=commands, distribute_method=distribute_method,
	#                                    parse_result=parse_result, sleeptime=sleeptime, sge_resources=sge_resources, # SGE options
	#                                    cluster_size=cluster_size, parallel_options=list(.options.multicore=list(preschedule=FALSE)), # PARALLEL options
	#                                    save.command=save.command)
	
	doMC::registerDoMC(cluster_size)
	contacts <- llply(contacts, function(D) domainClassifyR:::get.fend.samples.function(domain=D, n=n, random_seed=random_seed), .parallel=ifelse(cluster_size>1,TRUE,FALSE), .progress=ifelse(cluster_size==1,'text','none'))

	names(contacts) <- sapply(contacts, function(x) as.character(x$DOMAIN$intervalID))
	contacts
}

get.fend.samples.function <- function(domain, n, random_seed=FALSE) {
	if(!random_seed)
		set.seed(gsub(':','.',domain$DOMAIN$intervalID))

	sample_size <- max(domain$N_CONTACTS$HIGH_SCORE, 1)

	# define 2 lists to fill with random-sample fend distributions
	domain$CONTACTS$RANDOM_SAMPLE_TABLES <- list()

	# sample the observed fend pairs, if there are any
	if(domain$N_CONTACTS$ALL>1) {
		domain$CONTACTS$RANDOM_SAMPLE_TABLES$OBSERVED <- vector('list', n)
			for(i in seq(n))
				domain$CONTACTS$RANDOM_SAMPLE_TABLES$OBSERVED[[i]] <- table(sample(domain$CONTACTS$ALL$POSITION, size=sample_size))
	} else {
		domain$CONTACTS$RANDOM_SAMPLE_TABLES$OBSERVED <- list(table(domain$CONTACTS$ALL$POSITION))
	}

	# sample the potential fend pairs, if there are any
	if(!is.null(domain$N_FENDS.2D)) {
		domain$CONTACTS$RANDOM_SAMPLE_TABLES$POTENTIAL <- vector('list', n)
		for(i in seq(n))
			domain$CONTACTS$RANDOM_SAMPLE_TABLES$POTENTIAL[[i]] <- table(sample(domain$FENDS.2D$POSITION, size=sample_size))
	}
	
	# convert the random distribution list to a data.frame
	domain$CONTACTS$RANDOM_SAMPLE_TABLES <- lapply(domain$CONTACTS$RANDOM_SAMPLE_TABLES, do.call, what=rbind)

	domain
}


get.fend.samples.old <- function(contacts, n=1000, cluster.size=20, ...) {
	message(sprintf('[get.fend.samples] Randomly selecting contacts from ALL [x%d]', n))

	message(sprintf('[get.fend.samples] Making a local cluster of %d nodes...', cluster.size), appendLF=FALSE)
	snow.cluster <- snow::makeCluster(cluster.size, 'SOCK')
	on.exit({message('[get.fend.samples] Shutting the cluster down...');snow::stopCluster(snow.cluster)}, add=TRUE)
	message('DONE')

	parLapply(snow.cluster, contacts, function(domain, n) {
		domain$CONTACTS$RANDOM_SAMPLE_TABLES <- vector('list', n)
		for(i in seq(n)) {
			R <- table(domain$CONTACTS$ALL$POSITION[sample(x=seq(domain$N_CONTACTS$ALL), size=max(domain$N_CONTACTS$HIGH_SCORE, 1))])
			domain$CONTACTS$RANDOM_SAMPLE_TABLES[[i]] <- c(CORNER=0, FORWARD=0, REVERSE=0, OTHER=0)
			domain$CONTACTS$RANDOM_SAMPLE_TABLES[[i]][names(R)] <- R
		}
		domain$CONTACTS$RANDOM_SAMPLE_TABLES <- do.call(rbind, domain$CONTACTS$RANDOM_SAMPLE_TABLES)

		domain$Z_STATISTICS <- list(MEAN=plyr::aaply(domain$CONTACTS$RANDOM_SAMPLE_TABLES, 2, mean),
		                            SDEV=plyr::aaply(domain$CONTACTS$RANDOM_SAMPLE_TABLES, 2, sd),
		                            OBS=domain$FEND.SCORE.DISTRIBUTION$HIGH_SCORE)
		domain$Z_STATISTICS$Z <- (domain$Z_STATISTICS$OBS-domain$Z_STATISTICS$MEAN)/domain$Z_STATISTICS$SDEV
		domain}, n=n)

}



