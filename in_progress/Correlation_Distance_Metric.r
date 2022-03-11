# efg, Stowers Institute for Medical Research
# 29 Dec 2005
#http://research.stowers.org/mcm/efg/R/Visualization/cor-cluster/index.htm
# Explore creating hierarchical clusters using data with "known"
# correlation values.

library(mvtnorm)

# Based on R-News posting by Robin Hankin, 16 Dec 2005
# http://tolstoy.newcastle.edu.au/R/help/05/12/17693.html
CorrelatedXY <- function(N, mean1,mean2, variance1,variance2, correlation)
{
  corvar <- correlation * sqrt(variance1*variance2)
  rmvnorm(n=N,mean=c(mean1,mean2),
          sigma=matrix(c(variance1, corvar,
                         corvar,    variance2), 2,2))
}


N <- 1000
mean1 <- 2
mean2 <- 5
variance1 <- 1
variance2 <- 2

set.seed(17)

Raw <- CorrelatedXY(N, mean1, mean2, variance1, variance2, 0.0)
for (i in 1:5)
{
  Raw <- cbind(Raw, CorrelatedXY(N, mean1, mean2, variance1, variance2,
                                 (-1)^i * 0.2*i) )
  
}
colnames(Raw) <- paste(rep(c("P", "P", "M", "M"), 3),
                       sprintf("%2.2d", c(0,0,2,2,4,4,6,6,8,8,10,10)),
                       rep( c("A", "B"), 6), sep="")


# Look at correlation matrix
round( cor(Raw), 4)

corRaw <- cor(Raw)

library(spatstat)  # "im" function
plot(im(corRaw[nrow(corRaw):1,]), main="Correlation Matrix Map")


## Use correlations between variables as "distance"

dissimilarity <- 1 - cor(Raw)
distance <- as.dist(dissimilarity)
round(distance, 4)
plot(hclust(distance), main="Dissimilarity = 1 - Correlation", xlab="")

# only changes scaling
distance <- as.dist( (1 - cor(Raw))/2 )
round(distance, 4)
plot( hclust(distance), main="Dissimilarity = (1 - Correlation)/2", xlab="")

dissimilarity <- 1 - abs(cor(Raw))
distance <- as.dist(dissimilarity)
round(distance, 4)
plot(hclust(distance), main="Dissimilarity = 1 - Abs(Correlation)", xlab="")

dissimilarity <- sqrt(1 - cor(Raw)^2)
distance <- as.dist(dissimilarity)
round(distance, 4)
plot(hclust(distance), main="Dissimilarity = Sqrt(1 - Correlation^2)", xlab="")


# Go back to "best" distance metric
dissimilarity <- 1 - abs(cor(Raw))
distance <- as.dist(dissimilarity)


library(cluster)
plot(agnes(distance))

plot( pam(distance, 2))
plot( pam(distance, 6))


library(Hmisc)

plot( varclus(Raw, similarity="spearman") )
plot( varclus(Raw, similarity="pearson") )


library(pvclust)

cluster.bootstrap <- pvclust(Raw, nboot=1000, method.dist="abscor")
plot(cluster.bootstrap)
pvrect(cluster.bootstrap)

# Default "correlation" methods does not work well here
cluster.bootstrap <- pvclust(Raw, nboot=1000, method.dist="correlation")
plot(cluster.bootstrap)
pvrect(cluster.bootstrap)