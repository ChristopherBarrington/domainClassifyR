% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_min_fend_pair_domains.R
\name{get_min_fend_pair_domains}
\alias{get_min_fend_pair_domains}
\title{Select domains with a minimum number of contacts}
\usage{
get_min_fend_pair_domains(contacts, min = 100, measure = "HIGH_SCORE")
}
\arguments{
\item{contacts}{The list of domains generated from count.contacts()}

\item{min}{Minimum number of fend pairs required in a domain}

\item{measure}{Element of the N_CONTACTS list to filter on}
}
\value{
The contacts list, with additional information
}
\description{
Filters the domain list to retain intervals that have sufficiently high contacts recorded in N_CONTACTS
}
