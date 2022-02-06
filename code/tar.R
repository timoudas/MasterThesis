##  tar.R
##  
##  Bruce E. Hansen
##  Department of Economics
##  Social Science Building
##  University of Wisconsin
##  Madison, WI 53706-1393
##  behansen@wisc.edu
##  http://www.ssc.wisc.edu/~bhansen/
##  
##  
##  This program estimates a p'th order threshold autoregression and tests the
##  hypothesis of a linear autoregression, using the statistics described in
##  "Inference when a nuisance parameter is not identified under the null
##  hypothesis."
##  
##  The procedure is a function of eight inputs, and takes the form
##  
##    output <- tar(y,p,omit,r1,r2,rep,out,graph)
##    output$tests
##    output$pvalues 
##  
##  The inputs are:
##  
##  y     =  the vector of observations
##  p     =  the order of the autoregression
##  omit  =  lags (below p) to omit from autoregression [0 implies an AR(p)]
##  r1    =  the lower quantile
##  r2    =  the upper quantile
##  rep   =  number of simulation replications
##  out   =  indicator variable for screen output
##  graph =  indicator variable for graphical output
##  
##  
##  Additional detail on inputs:
##  
##  p       should be a positive integer
##  
##  omit    if 0, then an AR(p) model is used
##          if a positive integer, then that AR lag is omitted from model
##          if a vector of positive integers, then all these lags are omitted
##  
##  out     should be either 0 or 1
##          if 0, then nothing is printed to the screen nor to an output file
##          if 1, then OLS estimates and test statistics are printed to the
##              screen and to an output file (if previously declared)
##  
##  graph   should be either 0 or 1
##          if 0, then no graphs are generated
##          if 1, then graphs of the test statistic sequences are displayed
##  
##  There are two output vectors, the first contains the test statistics, and the
##  second the estimated asymptotic p-values.  The statistics (in order) are:
##  SupLM, ExpLM, AveLM, SupLMs, ExpLMs, AveLMs.

tar <- function(dat,p,omit,r1,r2,rep,ot,gr){ 

    dat <- as.matrix(dat) 
    omit <- as.matrix(omit)

    # Construct Data #
    t <- nrow(dat)-p
    if (omit[1,1]==0) {k=p 
    }else{ro <- nrow(omit) 
          k <- p-ro}
    k1 <- k+1
    tk1 <- t-k1
    tk2 <- t-k1*2
    y <- as.matrix(dat[(1+p):(t+p)]) 
    x <- matrix(1,nrow=t,ncol=1)
    for (j in 1:p){
        if (sum(omit==j)==0){
            x <- cbind(x,dat[(1+p-j):(t+p-j)])
        }    
    }

    # Null Model #
    s <- t(x)%*%y
    m <- t(x)%*%x
    mi <- solve(m)
    bnull <- mi%*%s
    u <- y-x%*%bnull
    sig <- (t(u)%*%u)/tk1  

    if (ot==1) {
        kx <- as.matrix(colSums(t(x)*(mi%*%t(x)))) 
        xu <- x*((u/(1-kx))%*%matrix(1,nrow=1,ncol=k1))
        xu <- xu-colMeans(xu) 
        v <- t(xu)%*%xu       
        se <- sqrt(as.matrix(diag(mi%*%v%*%mi)))       
 
        ar <- seq(1,p,by=1)
        if (sum(omit>0)>0){ 
            sel <- matrix(0,nrow=p,ncol=1)
            sel[omit] <- matrix(1,nrow=ro,ncol=1)
            ar <- ar[sel==0]
        }   
        indv <- matrix(0,nrow=k,ncol=1)
        for (j in 1:k){
            indv[j] <- paste(c("y(t-"),ar[j],c(")"),sep="")
        }        
        indv <- rbind("C",indv)
        tindv <- format(indv, digits=6)
        tbnull <- format(bnull, digits=6)
        tse <- format(se, digits=6) 

        for (j in 1:k1){
            if (j==1) {
                cat(" OLS Estimation of Null Linear Model", "\n",            
                    "-----------------------------------", "\n",
                    " Variable","   ","Estimate","   ","S.E.", "\n")
            }
            cat(" ",tindv[j],"    ",tbnull[j],"   ",tse[j], "\n")
            if (j==k1) {
            cat("   ","\n", 
                "Residual Variance =", sig, "\n",
                "   ","\n") 
            }
        }                    
    }
   
    # Test Statistic Storage #

    t1 <- round(t*r1)
    t2 <- round(t*r2)
    tn <- t2-t1

    lm <- matrix(0,nrow=tn,ncol=k) 
    lms <- lm  
    sigr <- lm
    lmmax <- matrix(0,nrow=rep,ncol=1) 
    lmsum <- lmmax 
    lmsume <- lmmax
    smax <- lmmax 
    ssum <- lmmax 
    ssume <- lmmax

    sim <- matrix(rnorm(t*(rep)),t,rep)   
    xu <- x*(u%*%matrix(1,nrow=1,ncol=k1))
    xmi <- x%*%mi
    misim <- t(xmi)%*%sim
    missim <- mi%*%(t(xu)%*%sim)

    for (d in 1:k){   # Pick Threshold Variable #
        
        if (ot==1) cat("Searching over Threshold Variable:", d,"\n") 
        
        z <- rank((x[,d+1]))
        zr <- order(z,seq(1,t,by=1))
        
        xr <- x*((z<=t1)%*%matrix(1,nrow=1,ncol=k1))
        xz1 <- zr[1:t1]
        x1 <- x[xz1,]
        mr <- t(x1)%*%x1
        sim1 <- sim[xz1,]

        xrsim <- t(x1)%*%sim1
        xusim <- t(xu[xz1,])%*%sim1

        for (r in 1:tn){   # Pick Threshold #
            
            t1r <- t1+r 
            xz <- zr[1:t1r]
            xz2 <- zr[(t1r+1):t]
            zrt <- zr[t1r]
            xsr <- x[zrt,]
            xsr <- t(as.matrix(xsr))
            xr[zrt,1:k1] <- xsr
            mr <- mr+t(xsr)%*%xsr
            xt <- xr-xmi%*%mr
            mrir <- mr%*%mi%*%mr 
            xti <- solve(mr-mrir)
            xtu <- t(xr)%*%u
            br <- xti%*%xtu
            ur <- u-xt%*%br
            sigr[r,d] <- t(ur)%*%ur

            xrulm <- xt*(u%*%matrix(1,nrow=1,ncol=k1))
            vrlmt <- solve(t(xrulm)%*%xrulm)
            lm[r,d] <- t(xtu)%*%vrlmt%*%xtu   # LM Statistic #
            lms[r,d] <- (t(xtu)%*%xti%*%xtu)/sig   # Standard LM Statistic #  

            # Simulation Manipulations #

            xrsim <- xrsim+t(xsr)%*%sim[zrt,] 
            sims <- xrsim-mr%*%misim
            wrep <- as.matrix(colSums(sims*(xti%*%sims)))
            smax <- apply(rbind(t(smax),t(wrep)),2,max)
            ssum <- ssum+wrep
            ssume <- ssume+exp(wrep/2)
               
            xusim <- xusim+as.matrix(xu[zrt,])%*%t(as.matrix(sim[zrt,]))
            sims <- xusim-mr%*%missim
            wrep <- as.matrix(colSums(sims*(vrlmt%*%sims)))
            lmmax <- apply(rbind(t(lmmax),t(wrep)),2,max) 
            lmsum <- lmsum+wrep
            lmsume <- lmsume+exp(wrep/2)
        }
    }

    rm(sim,misim,missim)

    # Test Statistics #
    tests <- rbind(max(lm),log(mean(exp(lm/2))),mean(lm),
                   max(lms),log(mean(exp(lms/2))),mean(lms))

    # P-values #
    tnk <- tn*k
    ssum <- ssum/tnk
    lmsum <- lmsum/tnk
    ssume <- log(ssume/tnk)
    lmsume <- log(lmsume/tnk)

    pvalues <- as.matrix(colMeans(cbind(as.matrix(lmmax),lmsume,lmsum,
               as.matrix(smax),ssume,ssum)>
               matrix(1,nrow=rep,ncol=1)%*%t(tests)))
    rm(lmmax,lmsume,lmsum,smax,ssume,ssum)

    if (ot==1){
        dm <- order(apply(sigr,2,min))
        dm <- dm[1]
        rm <- order(sigr[,dm])
        rm <- rm[1]
        yx <- cbind(y,x)
        yx <- yx[order(x[,dm+1]),]
        ys <- yx[,1]
        xs <- yx[,(2:(k+2))]
        xr <- cbind(rbind(xs[1:(t1+rm),],matrix(0,nrow=(t-t1-rm),ncol=k1)),
                    rbind(matrix(0,nrow=(t1+rm),ncol=k1),xs[(t1+rm+1):t,]))  
        dri <- solve(t(xr)%*%xr)
        br <- dri%*%(t(xr)%*%ys)
        ur <- ys-xr%*%br
        kx <- colSums(t(xr)*(dri%*%t(xr)))
        xu <- xr*(ur/(1-kx))%*%matrix(1,nrow=1,ncol=ncol(xr)) 
        xu <- xu-colMeans(xu)
        vrw <-t(xu)%*%xu
        se <- sqrt(as.matrix(diag(dri%*%vrw%*%dri))) 
        sig <- (t(ur)%*%ur)/(t-2*k1)
        sig1 <- (t(ur[1:(t1+rm)])%*%ur[1:(t1+rm)])/(t1+rm-k1)
        sig2 <- (t(ur[(t1+rm+1):t])%*%ur[(t1+rm+1):t])/(t-t1-rm-k1)
        
        th <- yx[t1+rm,dm+2]
        cat("   ","\n",
            " Global Estimates", "\n",
            "-----------------------------------", "\n",
            " Threshold Variable Lag  ", dm, "\n",
            " Threshold Estimate      ", th, "\n",
            " Error Variance          ", sig, "\n", 
            "   ","\n")
 
        thn <- format(th, digits=6)
        reg1 <-  paste(c("y(t-"),dm,c(") < "),thn,sep="")
        reg2 <-  paste(c("y(t-"),dm,c(") > "),thn,sep="")
        tbr <- format(br, digits=6)
        tse <- format(se, digits=6) 

        for (j in 1:k1){
            if (j==1) {
                cat(" Regime 1: ", reg1, "\n",            
                    "-----------------------------------", "\n",
                    " Variable","   ","Estimate","   ","S.E.", "\n")
            }
            cat(" ",tindv[j],"    ",tbr[j],"   ",tse[j], "\n")
            if (j==k1) {
            cat("   ","\n", 
                "Regime 1 Error Variance =", sig1, "\n",
                "   ","\n")  
            }
        }
        
        tbr <- format(br[(k+2):(2*k1),], digits=6)
        tse <- format(se[(k+2):(2*k1),], digits=6) 
        for (j in 1:k1){
            if (j==1) {
                cat(" Regime 2: ", reg2, "\n",            
                    "-----------------------------------", "\n",
                    " Variable","   ","Estimate","   ","S.E.", "\n")
            }
            cat(" ",tindv[j],"    ",tbr[j],"   ",tse[j], "\n")
            if (j==k1) {
            cat("   ","\n", 
                "Regime 2 Error Variance =", sig2, "\n",
                "   ","\n")  
            }
        }
            
        cat(" Test Statistics and Estimated Asymptotic P-Values", "\n",
            "-------------------------------------------------", "\n",
            " Robust LM Statistics", "\n",
            "     SupLM", "  ", tests[1], "  ", pvalues[1], "\n", 
            "     ExpLM", "  ", tests[2], "  ", pvalues[2], "\n",  
            "     AveLM", "  ", tests[3], "  ", pvalues[3], "\n", 
            "      ", "\n",  
            " Standard LM Statistics", "\n",
            "     SupLMs", "  ", tests[4], "  ", pvalues[4], "\n",  
            "     ExpLMs", "  ", tests[5], "  ", pvalues[5], "\n",  
            "     AveLMs", "  ", tests[6], "  ", pvalues[6], "\n") 
    }
          
    if (gr==1) {
        for (d in 1:k){
            z <- sort(yx[,d+2])
            zt <- z[(t1+1):t2]
            xxlim <- range(zt)
            yylim <- range(cbind(lm[,d],lms[,d]))   
            mtitle <- paste(c("Figure "),d,sep="")
            x11()
            plot(zt,lms[,d],lty=1,xlim=xxlim,ylim=yylim,type="l",ann=0)
            lines(zt,lm[,d],lty=2)            
            title(main=mtitle,xlab=tindv[d+1],ylab="Pointwise LM Statistics") 
            legend(xxlim[1],yylim[2],c("Standard","Hetero-Robust"),lty=c(1,2))          
        }
    }
    list(tests=tests,pvalues=pvalues)
}


