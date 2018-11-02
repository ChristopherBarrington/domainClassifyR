# domainClassifyR

## Installation

Please have 'misha' R package installed prior to running run_domainClassifyR.R, with a GENOME_DB configured (minimum example provided in example_data).

Install `domainClassifyR` using `devtools::github('ChristopherBarrington/domainClassifyR')`. Tested under Linux R 3.5.1 with:

```
Packages ----------------------------------------------------------------------
 package         * version    date
 devtools        * 1.13.6     2018-06-27
 domainClassifyR * 0.0.0.9000 2018-11-02
 doMC              1.3.5      2017-12-12
 dplyr           * 0.7.7      2018-10-16
 ggplot2         * 3.1.0      2018-10-25
 magrittr          1.5        2014-11-22
 misha           * 4.0.4      2018-09-24
 plyr            * 1.8.4      2016-06-08
 purrr             0.2.5      2018-05-29
 scales            1.0.0      2018-08-09
 snow            * 0.4-2      2016-10-14
 tibble            1.4.2      2018-01-22
 tidyselect        0.2.5      2018-10-11
```

## MWE:

Using example data, script should complete in minutes.

```
library(misha)
library(domainClassifyR)

# set the misha GENOME_DB root (eg to the mm10 directory included in the example_data)
gsetroot('/absolute/path/to/GENOME_DB/mm10/')

# load data / set parameters
intervals_2d <- gintervals.load("intervs.chr19_domains")
track_nm <- 'hic.example_project.example_dataset.tracks.score'
cache_path <- file.path(getwd(), 'domainclassifyr_cache')

# run the parts of the package
contacts <- get_contacts(domains=intervals_2d,
                         track_nm=track_nm,
                         cache.path=cache_path,
                         band=c(10e3,15e6),
                         get.fends=FALSE)
contacts <- add.domain.features(contacts, resolution=60e3)
contacts <- get.high.score.fends(contacts, min_score=40)
contacts <- count.contacts(contacts)
contacts <- get.domain.position.classifications(contacts)
contacts <- get.fend.samples(contacts, n=1000)
contacts <- get_Z_statistics(contacts)
contacts <- get_min_fend_pair_domains(contacts, min=100)
contacts <- get_Q_values(contacts)

# get a data.frame of qvalues for each sector of the domain
> ldply(contacts, function(d) d$Z_STATISTICS$OBSERVED$Q) %>% tibble::as.tibble()
# # A tibble: 77 x 5
#    .id         CORNER  FORWARD   REVERSE    OTHER
#    <chr>        <dbl>    <dbl>     <dbl>    <dbl>
#  1 2012:2012 1.00e+ 0 1.48e- 4 0.        1.00e+ 0
#  2 2013:2013 1.00e+ 0 1.00e+ 0 1.00e+  0 3.30e-21
#  3 2015:2015 1.00e+ 0 4.75e-65 1.00e+  0 1.00e+ 0
#  4 2017:2017 1.32e-52 1.00e+ 0 8.29e-  2 8.92e- 1
#  5 2019:2019 1.27e- 8 1.00e+ 0 9.79e-209 1.00e+ 0
#  6 2021:2021 4.54e-14 3.54e- 3 1.00e+  0 1.00e+ 0
#  7 2022:2022 1.00e+ 0 1.00e+ 0 1.00e+  0 6.35e-60
#  8 2023:2023 1.00e+ 0 1.00e+ 0 1.00e+  0 1.33e-11
#  9 2024:2024 4.34e-62 1.00e+ 0 1.15e-186 1.00e+ 0
# 10 2026:2026 5.88e-55 2.54e-25 1.00e+  0 1.00e+ 0
# # ... with 67 more rows

# plot the contacts in one domain
ggplot2::ggplot(contacts[['2106:2106']]$CONTACTS$HIGH_SCORE)+
  aes(x=FEND.X,y=FEND.Y,colour=SCORE)+
  geom_point(shape=20)+coord_fixed()+theme(aspect.ratio=1)
```
