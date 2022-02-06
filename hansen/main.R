##  main.R
##  
##  Bruce E. Hansen
##  Department of Economics
##  Social Science Building
##  University of Wisconsin
##  Madison, WI 53706-1393
##  behansen@wisc.edu
##  http://www.ssc.wisc.edu/~bhansen/
##   
##  This program replicates the empirical work reported in
##  "Inference when a nuisance parameter is not identified
##  under the null hypothesis."
## 
##  It calls the procedures tar.R and star.R and loads
##  the data gnp.dat.
##  
##  Refer to the procedure files for explanations.
 
gnp <- read.table("z:\\gnp.dat")
yg <- log(gnp[,1])
y <- (yg[2:175]-yg[1:174])*400
omit <- matrix(c(3,4),2,1)

source("Z:\\tar.R")
out1 <- tar(y,5,omit,0.15,0.85,1000,1,1)
source("Z:\\star.R")
out2 <- star(y,5,omit)
cat("STAR test: ", out2$s3, out2$ps3, "\n")
save.image(file = "tarout.RData")



