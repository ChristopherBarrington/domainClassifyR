#' Add information to the domain contacts entry
#' 
#' @param contacts The list of domains generated from get.contacts()
#' @param resolution The size of the region at the start/end of the domain to define a corner
#' 
#' @return The contacts list, with additional information
#' 
#' @export
add.domain.features <- function(contacts, resolution=50e3) {
	message('[add.domain.features] Adding domain features')

	check_order_2D_intervals(ldply(contacts, function(x) x$DOMAIN)[-1])

	progress_bar <- txtProgressBar(min=0, max=length(contacts), initial=0, style=3) 
	for(i in 1:length(contacts)) {
		domain <- contacts[[i]]
		d <- domain$DOMAIN

		if(is.null(domain$DOMAIN$SIZE1))
			domain$DOMAIN$SIZE1 <- domain$DOMAIN$end1-domain$DOMAIN$start1+1
		if(is.null(domain$DOMAIN$SIZE2))
			domain$DOMAIN$SIZE2 <- domain$DOMAIN$end2-domain$DOMAIN$start2+1

		domain$DOMAIN$AREA <- domain$DOMAIN$SIZE1*domain$DOMAIN$SIZE2

		# d$MIN_BORDER_INTERACTOR_FEND_DISTANCE <- d$end-d$start-resolution-resolution
		
		# define psotions of observed contacts (scores)
		domain$CONTACTS$ALL$GROUP.X <- floor((domain$CONTACTS$ALL$FEND.X-d$start1)/resolution) # resolution bin starting from DOMAIN.START
		domain$CONTACTS$ALL$GROUP.Y <- floor((d$end2-domain$CONTACTS$ALL$FEND.Y)/resolution) # resolution bin starting from DOMAIN.END
		domain$CONTACTS$ALL$FEND.DISTANCE <- domain$CONTACTS$ALL$FEND.Y-domain$CONTACTS$ALL$FEND.X

		domain$CONTACTS$ALL$POSITION <- factor(rep('OTHER', nrow(domain$CONTACTS$ALL)), levels=c('CORNER','FORWARD','REVERSE','OTHER'))
		domain$CONTACTS$ALL$POSITION[domain$CONTACTS$ALL$GROUP.X==0 & domain$CONTACTS$ALL$GROUP.Y==0] <- 'CORNER'
		domain$CONTACTS$ALL$POSITION[domain$CONTACTS$ALL$GROUP.X==0 & domain$CONTACTS$ALL$GROUP.Y>0] <- 'FORWARD'
		domain$CONTACTS$ALL$POSITION[domain$CONTACTS$ALL$GROUP.X>0 & domain$CONTACTS$ALL$GROUP.Y==0] <- 'REVERSE'

		# define positions of potential contacts (fends)
		if(!is.null(domain$FENDS.2D)) {
			x <- as.data.frame(domain$FENDS.2D)
			x <- data.frame(x,
			                GROUP.X=floor((x$FEND.X-d$start1)/resolution), # resolution bin starting from DOMAIN.START
			                GROUP.Y=floor((d$end2-x$FEND.Y)/resolution), # resolution bin starting from DOMAIN.END
			                FEND.DISTANCE=x$FEND.Y-x$FEND.X,
			                POSITION=factor(rep('OTHER', nrow(x)), levels=c('CORNER','FORWARD','REVERSE','OTHER')))
			
			x$POSITION[x$GROUP.X==0 & x$GROUP.Y==0] <- 'CORNER'
			x$POSITION[x$GROUP.X==0 & x$GROUP.Y>0] <- 'FORWARD'
			x$POSITION[x$GROUP.X>0 & x$GROUP.Y==0] <- 'REVERSE'

			domain$FENDS.2D <- x
		}
		contacts[[i]] <- domain
		setTxtProgressBar(progress_bar, i)
	}
	message()
	domain <- NULL
	gc()

	contacts
}

#' Check that the 2D intervals have interval1 is upstream of interval2
#'
check_order_2D_intervals <- function(intervals) {
	mid.intervals1 <- rowMeans(intervals[c('start1','end1')])
	mid.intervals2 <- rowMeans(intervals[c('start2','end2')])
	if(any(mid.intervals1>mid.intervals2))
		stop('Check that the domains are ordered such that interval1 is upstream of interval2')
}





add.domain.features.bak <- function(contacts, resolution=50e3) {
	message('[add.domain.features] Adding domain features')

	check_order_2D_intervals(ldply(contacts, function(x) x$DOMAIN)[-1])

	llply(contacts, function(domain) {
		d <- domain$DOMAIN
		o <- domain$CONTACTS$ALL

		if(is.null(d$SIZE1))
			d$SIZE1 <- d$end1-d$start1+1
		if(is.null(d$SIZE2))
			d$SIZE1 <- d$end2-d$start2+1

		d$AREA <- d$SIZE1*d$SIZE2

		# d$MIN_BORDER_INTERACTOR_FEND_DISTANCE <- d$end-d$start-resolution-resolution
		
		# define psotions of observed contacts (scores)
		domain$CONTACTS$ALL$GROUP.X <- floor((domain$CONTACTS$ALL$FEND.X-d$start1)/resolution) # resolution bin starting from DOMAIN.START
		domain$CONTACTS$ALL$GROUP.Y <- floor((d$end2-domain$CONTACTS$ALL$FEND.Y)/resolution) # resolution bin starting from DOMAIN.END
		domain$CONTACTS$ALL$FEND.DISTANCE <- domain$CONTACTS$ALL$FEND.Y-domain$CONTACTS$ALL$FEND.X

		domain$CONTACTS$ALL$POSITION <- factor(rep('OTHER', nrow(domain$CONTACTS$ALL)), levels=c('CORNER','FORWARD','REVERSE','OTHER'))
		domain$CONTACTS$ALL$POSITION[domain$CONTACTS$ALL$GROUP.X==0 & domain$CONTACTS$ALL$GROUP.Y==0] <- 'CORNER'
		domain$CONTACTS$ALL$POSITION[domain$CONTACTS$ALL$GROUP.X==0 & domain$CONTACTS$ALL$GROUP.Y>0] <- 'FORWARD'
		domain$CONTACTS$ALL$POSITION[domain$CONTACTS$ALL$GROUP.X>0 & domain$CONTACTS$ALL$GROUP.Y==0] <- 'REVERSE'

		# define positions of potential contacts (fends)

		if(!is.null(domain$FENDS.2D)) {
			x <- as.data.frame(domain$FENDS.2D)

			x$GROUP.X <- floor((x$FEND.X-d$start1)/resolution) # resolution bin starting from DOMAIN.START
			x$GROUP.Y <- floor((d$end2-x$FEND.Y)/resolution) # resolution bin starting from DOMAIN.END
			x$FEND.DISTANCE <- x$FEND.Y-x$FEND.X
			x$POSITION <- factor(rep('OTHER', nrow(x)), levels=c('CORNER','FORWARD','REVERSE','OTHER'))
			x$POSITION[x$GROUP.X==0 & x$GROUP.Y==0] <- 'CORNER'
			x$POSITION[x$GROUP.X==0 & x$GROUP.Y>0] <- 'FORWARD'
			x$POSITION[x$GROUP.X>0 & x$GROUP.Y==0] <- 'REVERSE'

			domain$FENDS.2D <- x
		}



#		if(!is.null(domain$FENDS.2D))
#			domain$FENDS.2D <- lapply(domain$FENDS.2D, function(x) {
#				x$GROUP.X <- floor((x$FEND.X-d$start1)/resolution) # resolution bin starting from DOMAIN.START
#				x$GROUP.Y <- floor((d$end2-x$FEND.Y)/resolution) # resolution bin starting from DOMAIN.END
#				x$FEND.DISTANCE <- x$FEND.Y-x$FEND.X
#				x$POSITION <- factor(rep('OTHER', nrow(x)), levels=c('CORNER','FORWARD','REVERSE','OTHER'))
#				x$POSITION[x$GROUP.X==0 & x$GROUP.Y==0] <- 'CORNER'
#				x$POSITION[x$GROUP.X==0 & x$GROUP.Y>0] <- 'FORWARD'
#				x$POSITION[x$GROUP.X>0 & x$GROUP.Y==0] <- 'REVERSE'
#				x})

		# define a list of distance measures that can be density'd etc
###		domain$DISTANCES <- lapply(domain$CONTACTS, function(contact.class) {
###			contact.class.distances <- list(FEND.DISTANCE=contact.class$FEND.Y-contact.class$FEND.X,     # distance between fends (specific interactor=high values)
###				                            RECIPROCAL.FEND.DISTANCE=d$SIZE-(contact.class$FEND.Y-contact.class$FEND.X), # amount of domain not spanned by fend pair (specific interactor=low values)
###			                                RECIPROCAL.FORWARD.DISTANCE=contact.class$FEND.X-d$start,    # distance from upstream domain border (forward extrusion=low values)
###			                                RECIPROCAL.REVERSE.DISTANCE=d$end-contact.class$FEND.Y)      # distance from downstream border (reverse extrusion=low values)
###			contact.class.distances})
		domain$DOMAIN <- d

		domain}, .progress='text')
}


