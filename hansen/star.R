##  star.R
##
##  This program is to calculate the approximate LM test for a STAR model 
##  
##  Bruce E. Hansen
##  Department of Economics
##  Social Science Building
##  University of Wisconsin
##  Madison, WI 53706-1393
##  behansen@wisc.edu
##  http://www.ssc.wisc.edu/~bhansen/

star <- function(dat,p,omit){ 

    dat <- as.matrix(dat) 
    omit <- as.matrix(omit)

    # Construct Data #
    t <- nrow(dat)-p
    if (omit[1,1]==0) {k=p 
    }else{ro <- nrow(omit) 
          k <- p-ro}
    k1 <- k+1
    y <- as.matrix(dat[(1+p):(t+p)]) 
    x <- matrix(1,nrow=t,ncol=1)
    for (j in 1:p){
        if (sum(omit==j)==0){
            x <- cbind(x,dat[(1+p-j):(t+p-j)])
        }    
    }

    # Null Model #
    u <- y-x%*%solve(t(x)%*%x)%*%t(x)%*%y
    uu <- t(u)%*%u

    f <- cbind(x,(x[,2:k1]^2),(x[,2:k1]^3))
    if (k>1){ 
        j <- 2
        while (j <= k1){
            j1 <- j+1 
            while (j1 <= k1){
                f=cbind(f,(x[,j]*x[,j1]))
            j1 <- j1+1
            }
        j <- j+1
        }
    } 
    e <- u-f%*%solve(t(f)%*%f)%*%t(f)%*%u
    ee <- t(e)%*%e
    s3 <- t*(uu-ee)/uu
    ps3 <- 1-pchisq(s3,(k1*k/2)+k)
    list(s3=s3,ps3=ps3)  
}


