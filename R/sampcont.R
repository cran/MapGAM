#***********************************************************************************
#
# Unmatched Control Sampling
# Copyright (C) 2016, The University of California, Irvine
#
# This library is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this library??? if not, <http://www.gnu.org/licenses/>.
#
#*******************************************************************************
sampcont <- function(rdata, type="stratified", casecol=1, Xcol=2, Ycol=3, regions=NULL, 
                     addstrat=NULL, times=NULL, n=1, nrow=100, ncol=100) { 
  if (type=="stratified") {
	if (!is.null(times)) {
	  warning("times argument has been renamed; please use addstrat instead in the future")
	  if (!is.null(addstrat)) stop("use either addstrat or times, not both")
	  addstrat = times
	}
  	if (!is.null(addstrat) && length(addstrat)!=dim(rdata)[1]) 
	  stop("addstrat argument must be NULL or a vector/factor of length = rows in rdata")
    if (is.null(regions)) {
      XYnames = names(rdata)[c(Xcol,Ycol)]
      Xrange = range(rdata[,Xcol])
      Yrange = range(rdata[,Ycol])
      polyGrid = PBSmapping::makeGrid(x=seq(Xrange[1],Xrange[2],length.out=nrow),
      		y=seq(Yrange[1],Yrange[2],length.out=ncol))
      names(rdata)[c(Xcol,Ycol)] = c("X","Y")
      rdata$EID = 1:length(rdata$X)
      idpolys = PBSmapping::findCells(PBSmapping::as.EventData(rdata),polyGrid)
      rdata = merge(idpolys,rdata)
      regions = paste((1:nrow)[rdata$PID],(1:ncol)[rdata$SID],sep=",")
      rdata = rdata[,!(names(rdata) %in% c("EID","PID","SID","Bdry"))]
      names(rdata)[c(Xcol,Ycol)] = XYnames
      gridsize = c(nrow,ncol)
    } else {
 	  gridsize = NA 
	  if (length(regions)!=dim(rdata)[1]) 
	 		 stop("regions argument must be NULL or a vector/factor of length = rows in rdata")
	  }
    strata = paste(regions,addstrat,sep="")
    eligible = rdata[,casecol]==0		  	# eligible controls
    counts = table(strata[eligible])		# number of eligible controls in each stratum
    if (median(counts)<=n) warning(paste("Most strata include",n,"or fewer controls;", 
    		"consider alternative sampling designs if this is unexpected"))
    estrat = names(counts)			# strata with at least 1 eligible control
    ns = pmin(n,counts)	    			# number to sample from each stratum
    dsamp = rdata[rdata[,casecol]==1,]	# take all cases  	
    ncases = dim(dsamp)[1]				# number of cases
    dsamp[ncases + 1:sum(ns),] = NA		# placeholder for controls 
    w = c(rep(1,ncases),rep(counts/ns,times=ns))	# weights are 1 for cases, IPW for controls
    # sample controls from each non-empty stratum and add to dsamp data frame
    rownum = ncases
    for (i in 1:length(counts)) {
      ind = sample(1:counts[i],ns[i])		# sample without replacement
      controls = rdata[strata==estrat[i] & eligible,][ind,]
      dsamp[rownum + 1:ns[i],] = controls
      rownum = rownum + ns[i]
    }
	n = length(w)-ncases		# redefine n as total numer of controls
    cat(paste(n,"controls selected from",sum(eligible),
	"eligibles in",length(estrat),"strata."),fill=T)
  }
  if (type=="simple") {
    eligible = rdata[,casecol]==0		  		# eligible controls
    if (n > sum(eligible)) stop(paste("rdata contains only ",n," eligible controls"))
    dsamp = rdata[rdata[,casecol]==1,]			# take all cases  	
    ind = sample(1:dim(rdata)[1],n,prob=eligible)	# simple random sample
	dsamp = rbind(dsamp,rdata[ind,])
	w = c(rep(1,sum(!eligible)),rep(n/sum(eligible),n))
	gridsize = c(1,1)
  }
  return(list(rdata=dsamp,w=w,ncont=n,type=type,gridsize=gridsize))
}
