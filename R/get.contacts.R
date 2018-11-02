#' Create a list of information for a set of intervals
#'
#' Extracts interaction frequency information for a set of 2D intervals. The domain information is cached
#' for quick loading, or loaded from the cache if the parameters have been run (and cached) before.
#'
#' @param domains A data.frame of 2D intervals
#' @param track_nm Name of the score track from which fend pairs are extracted
#' @param band Genomic distance between fend pairs to consider
#' @param cache.path Absolute path to the directory that cached domains can be found or saved
#' @param distribute_method Specify where the jobs to define each interval/domain should be distributed (using distributR)
#' @param sge.resources,max_jobs,cluster_size Parameters passed to distributR
#' @param get.fends Extracting the valid fragment ends and identifying potential 2D interactors (uses masses of memory!)
#' 
#' @return The contacts list, with additional information
#'
#' @seealso \code{distributR} for information on distributed jobs.
#'
#' @export
get_contacts <- function(domains,
                         track_nm='hic.Genome_Structure_Development.AST_NS_WT.tracks.observed_new_score',
                         band=c(10e3,15e6),
                         cache.path='/home/cbarring/CLASSIFY_DOMAINS.CACHE',
                         distribute_method='PARALLEL', sge.resources=list(h_rt='00:30:00'), max_jobs=500, cluster_size=10,
                         get.fends=FALSE) {
	if(class(domains)!='data.frame')
		stop('[get.contacts] The domains argument needs to be a data.frame of gintervals (1D or 2D)!')

	message(sprintf('[get.contacts] Getting contact matrices for %d domains', nrow(domains)))

	if(!dir.exists(cache.path))
		dir.create(cache.path, recursive=TRUE, showWarnings=FALSE)

	if(is.null(domains$intervalID))
		domains$intervalID <- seq(nrow(domains))

	domains$track_nm <- track_nm
	domains$min.dist <- band[1]
	domains$max.dist <- band[2]
	domains$get.fends <- get.fends

	cache.files.key <- apply(domains, 1, function(x) paste(gsub(' ','',x), collapse='_'))
	cache.files <- sprintf('%s/%s.RData', cache.path, sapply(cache.files.key, digest::digest))
	cache.files.missing <- !file.exists(cache.files)

	domains$cache.file.key <- cache.files.key
	domains$cache.file <- cache.files

	ret.val <- list()
	if(any(cache.files.missing)) {
		message(sprintf('[get.contacts] Saving new domains to cache [n=%d]', sum(cache.files.missing)))
		# ret.val <- distributR::distributr(parameters_table=domains[cache.files.missing,],
		#                                   function_name='domainClassifyR:::get.contact.function',
		#                                   distribute_method=distribute_method,
		#                                   sge_resources=sge.resources,
		#                                   cluster_size=cluster_size,
		#                                   jobs_name='CLASSIFY_DOMAINS',
		#                                   max_jobs=max_jobs,
		#                                   dry_run=FALSE)

		check_parameters_table_format <- function(parameters_table) {
			factor.cols <- sapply(parameters_table, is.factor)
			parameters_table[factor.cols] <- lapply(parameters_table[factor.cols], function(x) levels(x)[x])
			
			numeric.args <- sapply(parameters_table, is.numeric)
			parameters_table[!numeric.args] <- lapply(parameters_table[!numeric.args], sprintf, fmt="'%s'")

			parameters_table
		}

		collapse_parameters_table <- function(parameters_table)
			apply(parameters_table, 1, function(p) paste(names(p), unlist(p), sep='=', collapse=','))




		parameters_table <- check_parameters_table_format(domains[cache.files.missing,])
		parameters_set <- collapse_parameters_table(parameters_table)
		commands <- sprintf('%s(%s)', 'domainClassifyR:::get.contact.function', parameters_set)
		ret.val <- plyr::llply(commands, function(cmd) eval(expr=parse(text=cmd)), .parallel=FALSE, .progress='text')
	} else {
		message('[get.contacts] All domains saved in cache')
	}

	if(any(!file.exists(cache.files))) {
		print(cache.files[!file.exists(cache.files)])
		stop(sprintf('[get.contacts] Some (%d) cache files were not found!!', sum(!file.exists(cache.files))))
	}

	message(sprintf('[get.contacts] Preallocating a %d element list', length(cache.files)))
	ret.val <- vector('list', length(cache.files))

	message(sprintf('[get.contacts] Loading matrices from cache (%s)', cache.path))
	pb <- txtProgressBar(min=0, max=length(cache.files), initial=0, style=3) 
	for(i in 1:length(cache.files)) {
		load(cache.files[i])
		ret.val[[i]] <- contacts
		ret.val[[i]]$CACHE.FILE <- cache.files[i]
		setTxtProgressBar(pb,i)
	}
	contacts <- NULL
	names(ret.val) <- sapply(ret.val, function(x) as.character(x$DOMAIN$intervalID))
	message()

	if(any(sapply(ret.val, length)==1))
		stop('[get.contacts] One or more of the domains did not have all the information!')

	ret.val
}

get.contact.function <- function(domain,
                                 chrom1=domain$chrom1, start1=domain$start1, end1=domain$end1,
                                 chrom2=domain$chrom2, start2=domain$start2, end2=domain$end2,
                                 intervalID=domain$intervalID,
                                 track_nm='hic.Genome_Structure_Development.AST_NS_WT.tracks.observed_new_score',
                                 # band=c(10e3,15e6),
                                 min.dist=10e3, max.dist=15e6,
                                 cache.file=NULL, cache.file.key=NULL,
                                 force=FALSE,
                                 get.fends=TRUE,
                                 ...) {
	if(missing(domain))
		domain <- data.frame(chrom1=chrom1, start1=start1, end1=end1,
		                     chrom2=chrom2, start2=start2, end2=end2,
		                     intervalID=intervalID,
		                     list(...))
	band <- c(min.dist, max.dist)
	
	cache.file.log <- file(sprintf('%s.log', cache.file), open='w')
	sink(cache.file.log, type='message')		
	sink(cache.file.log, type='output', append=TRUE)		

	on.exit(sink(NULL, type='message'), add=TRUE)
	on.exit(sink(NULL, type='output'), add=TRUE)
	on.exit(close(cache.file.log), add=TRUE)

	message('[get.contact.function] Querying for: ', cache.file)
	message(sprintf('[get.contact.function] Host: %s', Sys.getenv('HOSTNAME')))
	message(sprintf('[get.contact.function] IntervalID: %s', as.character(domain$intervalID)))
	message(sprintf('[get.contact.function] Interval 1: %s %d %d', as.character(domain$chrom1), domain$start1, domain$end1))
	message(sprintf('[get.contact.function] Interval 2: %s %d %d', as.character(domain$chrom2), domain$start2, domain$end2))
	message(sprintf('[get.contact.function] Hi-C track: %s', track_nm))
	message(sprintf('[get.contact.function] Band: %d-%d', min.dist, max.dist))
	message(sprintf('[get.contact.function] Getting fends: ', get.fends))

	# get scores of observed contacts
	contacts <- gextract(track_nm,
	                     intervals=domain,
	                     colnames='SCORE',
	                     band=rev(-band))[c('start1','start2','SCORE')]

	if(is.null(contacts))
		contacts <- data.frame(start1=numeric(0), start2=numeric(0), SCORE=numeric(0))
	colnames(contacts) <- c('FEND.X','FEND.Y','SCORE')
	message(sprintf('[get.contact.function] Contacts extracted: %d', nrow(contacts)))

	if(get.fends) {
		# get theoretical number of interacting fends
		fend.trackname <- 'redb.HindIII'
		fends <- gscreen(sprintf('%s==1', fend.trackname), intervals=gintervals(domain$chrom1, c(domain$start1,domain$start2), c(domain$end1,domain$end2)))

		if(is.null(fends)) { # if there are no fends at all
			message('[get.contact.function] No fragment ends identified in the domain...')
			fends <- as.matrix(data.frame(start=numeric(),end=numeric()))
			fends.2D <- as.matrix(data.frame(FEND.X=numeric(), FEND.Y=numeric()))
		} else {
			message(sprintf('[get.contact.function] Total fragment ends identified: %d', nrow(fends)))		
			message('[get.contact.function] Identifying potential 2D interacting fragment ends within domain')
			get.fends.2D <- function(x) {
				i <- gintervals.interactors(x, min_dist=band[1], max_dist=band[2])
				if(nrow(i)>0) { # if there are fend interactors
					i <- i[c('start1','start2')]
				} else {
					i <- data.frame(start1=numeric(), start2=numeric())
				}
				colnames(i) <- c('FEND.X','FEND.Y')
				as.matrix(i)}
			fends.2D <- get.fends.2D(fends)
			fends <- as.matrix(fends[c('start','end')])
		}

		message(sprintf('[get.contact.function] Total possible fragment end interactions: %d', nrow(fends.2D)))
	}

	message('[get.contact.function] Saving contacts...', appendLF=FALSE)
	contacts <- list(DOMAIN=domain, CONTACTS=list(ALL=contacts))

	if(get.fends)
		contacts <- append(contacts, list(FENDS=fends, FENDS.2D=fends.2D))

	save(contacts, file=cache.file)
	message('DONE')
	message(sprintf('[get.contact.function] Size of contacts list: %s', format(object.size(contacts), units='auto')))
	message('[get.contact.function] Finished!')

	contacts
}

remove_2D_potential_fends <- function(contacts) {
	message('[count.contacts] Removing 2D fend pairs to save space')
	progress_bar <- txtProgressBar(min=0, max=length(contacts), initial=0, style=3) 

	for(i in 1:length(contacts)) {
		contacts[[i]]$FENDS.2D <- NULL

		setTxtProgressBar(progress_bar, i)
	}
	message()

	contacts

}
