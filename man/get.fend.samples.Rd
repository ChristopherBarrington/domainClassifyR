% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.fend.samples.R
\name{get.fend.samples}
\alias{get.fend.samples}
\title{Sample fragment ends randomly to get a null distribution}
\usage{
get.fend.samples(contacts, random_seed = FALSE, n = 1000,
  cluster_size = 20, distribute_method = "SGE", parse_result = TRUE,
  sleeptime = 300, sge_resources = list(), save.command = TRUE)
}
\arguments{
\item{contacts}{The list of domains generated from get.contacts()}

\item{n}{Number of times random selection is repeated}

\item{cluster.size}{Number of parallel threads to run}
}
\value{
The contacts list, with additional information
}
\description{
Randomly samples the total observed fend pairs and potential fend pairs to randomly select the same number of fend pairs that are
high-scoring. The random selection is repeated and a null distribution for domain sector scores determined.
}
\details{
Parallelisation is achieved using the distributR package and the PARALLEL distribution option.
}
