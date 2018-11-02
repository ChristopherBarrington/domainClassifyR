#' Load information on domain intervals
#' 
#' @param project,dataset Name of the project and dataset that contain compaction tracks
#' @param band_from,band_to Size of the band used to define the compaction track
#' 
#' @return A list with intervals for the CLUSTERS, DOMAINS and BORDERS
#' 
#' @export
load_domains <- function(dataset, project, band_from=100e3, band_to=400e3) {
	gdb.root <- sprintf('hic.%s.%s', project, dataset)
	clusters.track <- sprintf('%s.clusters.%d_%d_n_2', gdb.root, band_from, band_to)
	domains.intervals <- sprintf('%s.compaction.%d_%d_domains_top5', gdb.root, band_from, band_to)
	borders.intervals <- sprintf('%s.compaction.%d_%d_domain_borders_top5', gdb.root, band_from, band_to)

	if(!gtrack.exists(clusters.track))
		clusters.track <- sprintf('%s.clusters.n_2', gdb.root)

	message(sprintf('[load_domains] Loading borders: %s', borders.intervals))
	message(sprintf('[load_domains] Loading domains: %s', domains.intervals))
	message(sprintf('[load_domains] Loading clusters: %s', clusters.track))

	data <- list(CLUSTERS=gextract(clusters.track, intervals=ALLGENOME, colnames='CLASS'),
	             DOMAINS=gintervals.load(domains.intervals),
	             BORDERS=gintervals.load(borders.intervals))
	data
}
