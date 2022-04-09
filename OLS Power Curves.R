library(MASS)
library(nlme)
library("lmtest")
library("sandwich")
library(xtable)
library(ggplot2)

rm(list = ls())

# define given parameters: 
{
  n= c(100,250,500,1000,5000, 10000, 100000)
  R= 1000
  
  mu <- c(1, -0.7, 3.4, 0)
  covMat <- matrix(c(1.1, -0.2, -0.1, 0, 
                   -0.2, 0.9, -0.2, 0,
                   -0.1, -0.2, 1, 0, 
                   0,0,0,1), nrow=4, ncol=4, byrow=TRUE)
  beta <- c(-0.92, 0, 2.11, 1.3)
  
  criticalT <- rep(NA, 7)
  
  c <- seq(-1, 1, 0.1)
}


# define output table: 
#       c=-1 | c=-0.9 | ... | c=1|
#      ---------------------------
#n=100 |     |        |     |    |
# ...  |     |        |     |    |
#n=10^5|     |        |     |    |
outputTable <- as.data.frame(matrix(nrow=length(n), ncol=length(c), dimnames=list(paste("N=", n), paste("C=",c)) ))



#loop over sample sizes
for(i in 1:length(n)){
  
  #critical T wont change once we set a sample size. 
  criticalT[i] <- qt(p=0.05/2, df=n[i]- 4, lower.tail=FALSE)
  
  
  # loop over values in null hypothesis
  for(k in 1:length(c)){
  
    
    
    #this vector will contain 0/1 values depending on if the null H0: B1 = c is rejected or not (1=reject)
    # the proportion of times the null is rejected is the mean of this vector
    hypTest_temp <-matrix(nrow=R, ncol=1)
    
    
      #loop over R samples 
      for(j in 1:R){
        
        #generate X, Y
        X <- mvrnorm(n[i],mu,covMat)
        X <- cbind(rep(1, n[i] ), X)
        y <- X[,1:4]%*%beta + X[,5] #coefficients + error
        
        #get beta and se: (Assignment instructions do not prevent us from using lm(.) by the way)
        beta1 <- summary(lm(y~X[, 2:4]))$coefficients[2,1]
        se1<- summary(lm(y~X[, 2:4]))$coefficients[2,2]
        
        #generate calculated t-statistic for the k'th null: 
        calcT <- (beta1-c[k])/se1
        hypTest_temp[j] <- abs(calcT)>criticalT[i]
      }
    
    outputTable[i,k] <- mean(hypTest_temp)
    print(cat("done with null=", c[k]))
  }  
  print(cat("done with sample size=",n[i]))
}

plot1 <- ggplot() + geom_line(aes(x=as.vector(c), y=as.numeric(outputTable[1,]), colour="N=100" ), size=1) +
        geom_line(aes(x=as.vector(c), y=as.numeric(outputTable[2,]), colour="N=250" ), size=1) +
        geom_line(aes(x=as.vector(c), y=as.numeric(outputTable[3,]), colour="N=500" ), size=1) +
        geom_line(aes(x=as.vector(c), y=as.numeric(outputTable[4,]), colour="N=1000" ), size=1) +
        geom_line(aes(x=as.vector(c), y=as.numeric(outputTable[5,]), colour="N=5000" ), size=1) +
        geom_line(aes(x=as.vector(c), y=as.numeric(outputTable[6,]), colour="N=10000" ), size=1) +
        geom_line(aes(x=as.vector(c), y=as.numeric(outputTable[6,]), colour="N=10^5" ), size=1) +
        scale_color_manual(name = "Sample Sizes", values = c("N=100" = "red", "N=250" = "blue","N=500" = "green",
                                                         "N=1000" = "orange","N=5000" = "purple", "N=10000" = "cyan",
                                                         "N=10^5"="black")) + 
        theme_light() +
        ylab("Proportion of times the null is rejected") +
        xlab("C: null hypothesis value of B1")

xtable(t(outputTable))

