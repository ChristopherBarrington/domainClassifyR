# domainClassifyR

## MWE:

library(shaman)
library(misha)

intervals <- data.frame(chrom1='chr10', start1=50e6, end1=50.4e6, chrom2='chr10', start2=50e6, end2=50.4e6, intervalID='1', intervalType='TAD')

contacts <- get_contacts(domains=intervals,
                         track_nm='a.shaman.score.map',
                         cache.path=cache_path,
                         band=c(10e3,15e6),
                         get.fends=FALSE)
contacts <- add.domain.features(contacts, resolution=60e3)
contacts <- get.high.score.fends(contacts, min_score=40)
contacts <- count.contacts(contacts)
contacts <- get.domain.position.classifications(contacts)
contacts <- get.fend.samples(contacts, n=1000, cluster_size=30, distribute_method='none')
contacts <- get_Z_statistics(contacts)
contacts <- get_min_fend_pair_domains(contacts, min=100)
contacts <- get_Q_values(contacts)
