############################################################################################
# RSEQREP: RNA-Seq Reports, an open-source cloud-enabled framework for reproducible
# RNA-Seq data processing, analysis, and result reporting
# 
# https://github.com/emmesgit/RSEQREP
#
# Copyright (C) 2017 The Emmes Corporation 
# 
# This program is free software that contains third party software subject to various licenses, 
# namely, the GNU General Public License version 3 (or later), the GNU Affero General Public License 
# version 3 (or later), and the LaTeX Project Public License v.1.3(c). A list of the software contained 
# in this program, including the applicable licenses, can be accessed here: 
# 
# https://github.com/emmesgit/RSEQREP/SOFTWARE.xlsx  
# 
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# Program:  pvclust-max-dist.R
# Version:  RSEQREP 1.0.0
# Author:   Johannes Goll, Travis Jensen
# Purpose:  modifications to the pvclust algorithm.  Max distance cutoff added.
# Input:    N/A
# Output:  	N/A
############################################################################################

hc2split = function (x) 
{
	A <- x$merge
	n <- nrow(A) + 1
	B <- list()
	for (i in 1:(n - 1)) {
		ai <- A[i, 1]
		if (ai < 0) 
			B[[i]] <- -ai
		else B[[i]] <- B[[ai]]
		ai <- A[i, 2]
		if (ai < 0) 
			B[[i]] <- sort(c(B[[i]], -ai))
		else B[[i]] <- sort(c(B[[i]], B[[ai]]))
	}
	CC <- matrix(rep(0, n * (n - 1)), nrow = (n - 1), ncol = n)
	for (i in 1:(n - 1)) {
		bi <- B[[i]]
		m <- length(bi)
		for (j in 1:m) CC[i, bi[j]] <- 1
	}
	split <- list(pattern = apply(CC, 1, paste, collapse = ""), 
			member = B)
	return(split)
}


pvrect <- function(x, alpha=0.95, pv="au", type="geq", max.only=TRUE, border=2,max.dist=0.5, ...)
  {
    len <- nrow(x$edges)
    member <- hc2split(x$hclust)$member
    order  <- x$hclust$order
    usr <- par("usr")
    xwd <- usr[2] - usr[1]
    ywd <- usr[4] - usr[3]
    cin <- par()$cin

    ht <- c()
    j <- 1

    if(is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt"))))
       stop("Invalid type argument: see help(pvrect)")
    
    for(i in (len - 1):1)
      {
		 
        if     (pm==1) wh <- (x$edges[i,pv] >= alpha & x$hclust$height[i]<=max.dist) # Greater than or EQuals
        else if(pm==2) wh <- (x$edges[i,pv] <= alpha & x$hclust$height[i]<=max.dist) # Lower than or EQuals
        else if(pm==3) wh <- (x$edges[i,pv] >  alpha & x$hclust$height[i]<=max.dist) # Greater Than
        else if(pm==4) wh <- (x$edges[i,pv] >  alpha & x$hclust$height[i]<=max.dist) # Lower Than

        if(wh)
          {
            mi <- member[[i]]
            ma <- match(mi, order)
            
            if(max.only == FALSE || (max.only && sum(match(ma, ht, nomatch=0)) == 0))
              {
                xl <- min(ma)
                xr <- max(ma)
                yt <- x$hclust$height[i]
                yb <- usr[3]
                
                mx <- xwd / length(member) / 3
                my <- ywd / 200
                
                rect(xl - mx, yb + my, xr + mx, yt + my, border=border, shade=NULL, ...)
                
                j <- j + 1
              }
            ht <- c(ht, ma)
          }
      }
  }


pvpick <- function(x, alpha=0.95, pv="au", type="geq", max.only=TRUE,max.dist=0.5)
  {
    len <- nrow(x$edges)
    member <- hc2split(x$hclust)$member
    order  <- x$hclust$order
    
    ht <- c()
    a  <- list(clusters=list(), edges=c()); j <- 1

    if(is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt"))))
       stop("Invalid type argument: see help(pickup)")
    
    for(i in (len - 1):1)
      {
	  if     (pm==1) wh <- (x$edges[i,pv] >= alpha & x$hclust$height[i]<=max.dist) # Greater than or EQuals
	  else if(pm==2) wh <- (x$edges[i,pv] <= alpha & x$hclust$height[i]<=max.dist) # Lower than or EQuals
	  else if(pm==3) wh <- (x$edges[i,pv] >  alpha & x$hclust$height[i]<=max.dist) # Greater Than
	  else if(pm==4) wh <- (x$edges[i,pv] >  alpha & x$hclust$height[i]<=max.dist) # Lower Than


        if(wh)
          {
            mi <- member[[i]]
            ma <- match(mi, order)

            if(max.only == FALSE || (max.only && sum(match(ma, ht, nomatch=0)) == 0))
              {
                a$clusters[[j]] <- x$hclust$labels[mi]
                a$edges <- c(a$edges,i)
                
                j <- j + 1
              }
            ht <- c(ht, ma)
          }
      }
    
    a$edges <- a$edges[length(a$edges):1]
    a$clusters <- a$clusters[length(a$edges):1]

    return(a)
  }