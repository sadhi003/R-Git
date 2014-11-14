### R code from vignette source 'Presentation.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: init
###################################################
options(width=60)


###################################################
### code chunk number 2: Introduction
###################################################
3 + 3
# This is a comment
x <- 3 ; y = 3
x <- y <- 6
x^2 # square
x = c(1, 2, 3, 4) ; x = 1:4
x
x+10 # R is vectorized
x[4] # Accessing the fourth element


###################################################
### code chunk number 3: Introduction 2
###################################################
x>2
x[x>2]
x[1:3]
is.na(x)
cumsum(x)
summary(x)


###################################################
### code chunk number 4: introPass
###################################################
x <- 1:10
x
x * 2
x
y <- x
x * 2
y


###################################################
### code chunk number 5: math
###################################################
mean(x)
median(x)
sd(x)
sum(x)
sqrt(x)


###################################################
### code chunk number 6: math2
###################################################
summary(x)


###################################################
### code chunk number 7: help (eval = FALSE)
###################################################
help.start()


###################################################
### code chunk number 8: help (eval = FALSE)
###################################################
?plot


###################################################
### code chunk number 9: packages (eval = FALSE)
###################################################
install.packages('ggplot2')


###################################################
### code chunk number 10: bioc (eval = FALSE)
###################################################
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("ArrayExpress")


###################################################
### code chunk number 11: bioc2 (eval = FALSE)
###################################################
 # R
 update.packages(checkBuilt=TRUE, ask=FALSE)
 # bioconductor
 biocLite(character(), ask=FALSE)


###################################################
### code chunk number 12: vector
###################################################
 # Vector of strings
 x <- c("Lincoln", "Roosevelt", "Jackson")
 ## Replace the second element of the vector
 x[2] <- NA
 x
 x[!is.na(x)]


###################################################
### code chunk number 13: vector
###################################################
 ## character vectors
 x <- c("Lincoln", "Roosevelt", "Jackson"); class(x)
 x <- c(1.2, 1.5); class(x)
 x <- 1:4; class(x)
 x <- x>2; class(x)


###################################################
### code chunk number 14: matrix
###################################################
 matrix(1:10, ncol=2)


###################################################
### code chunk number 15: matrix2
###################################################
 matrix(1:10, nrow=2, byrow=TRUE)


###################################################
### code chunk number 16: data.frame
###################################################
 df <- data.frame(
   patientID = c("001", "002", "003"), 
   treatment = c("drug", "placebo", "drug"),
   age = c(20, 30, 24)
 )
 df
 subset(df, treatment=="drug")


###################################################
### code chunk number 17: list
###################################################
 my.list <- list(
   fruits = c("oranges", "bananas", "apples"), 
   mat = matrix(1:10, nrow=2)
 )
 my.list


###################################################
### code chunk number 18: read.table
###################################################
URL <- "http://www.stanford.edu/~druau/pivot_table.csv"
pivot <- read.table(URL, sep=',', header=TRUE)
head(pivot)
class(pivot)


###################################################
### code chunk number 19: foreign (eval = FALSE)
###################################################
help(package="foreign")


###################################################
### code chunk number 20: xlsx (eval = FALSE)
###################################################
?read.xlsx


###################################################
### code chunk number 21: html (eval = FALSE)
###################################################
 eq <- readHTMLTable("http://www.iris.edu/seismon/last30.html")


###################################################
### code chunk number 22: import1
###################################################
URL <- "http://www.stanford.edu/~druau/pivot_table.csv"
pivot <- read.table(URL, sep=',', header=TRUE)
pivot$value <- round(pivot$value, digits=3)
head(pivot, 3)
tail(pivot, 3)


###################################################
### code chunk number 23: format
###################################################
library(reshape)
head((pivot <- cast(pivot, gene ~ condition)), 2)


###################################################
### code chunk number 24: format
###################################################
head(subset(pivot, cheetos>=0.2), 2)


###################################################
### code chunk number 25: format5
###################################################
library(sqldf)
head(sqldf('SELECT * FROM pivot WHERE cheetos >= 0.2'), 2)


###################################################
### code chunk number 26: format6
###################################################
library(doBy)
summaryBy(Sepal.Width + Petal.Width ~ Species, data=iris, FUN=c(mean))


###################################################
### code chunk number 27: database (eval = FALSE)
###################################################
 require("RMySQL")
 con <- dbConnect(MySQL(), user="david", password="will_not_tell_you", 
 dbname="db_name", host="mysql_server.stanford.edu")
 results <- dbGetQuery(con, "SELECT * FROM pat_table LIMIT 5")


###################################################
### code chunk number 28: export (eval = FALSE)
###################################################
 save(list=c(pivot, mat), file="data_exoprt.Rda")


###################################################
### code chunk number 29: export1 (eval = FALSE)
###################################################
 load("data_export.Rda")


###################################################
### code chunk number 30: export1 (eval = FALSE)
###################################################
 write.table(mat, file='matrix.csv', sep=',')


###################################################
### code chunk number 31: export1 (eval = FALSE)
###################################################
 library(foreign)
 ?write.arff
 ?write.dta
 library(SASxport)
 ?write.xport


###################################################
### code chunk number 32: microarray (eval = FALSE)
###################################################
 library(affy); library(GEOquery); library(mouse4302cdf)
 getGEOSuppFiles("GSE12499")


###################################################
### code chunk number 33: microarray2 (eval = FALSE)
###################################################
 system('tar -xf GSE12499/GSE12499_RAW.tar -C GSE12499/')
 system('rm GSE12499/*.CHP*; rm GSE12499/*.tar')


###################################################
### code chunk number 34: microarray3 (eval = FALSE)
###################################################
 da <- ReadAffy(celfile.path="./GSE12499/", compress=TRUE)


###################################################
### code chunk number 35: microarray3
###################################################
da


###################################################
### code chunk number 36: ArrayExpress1 (eval = FALSE)
###################################################
 pneumoHS = queryAE(keywords = "pneumonia", species = "homo+sapiens")


###################################################
### code chunk number 37: ArrayExpress2
###################################################
pneumoHS[1:3, 1:3]


###################################################
### code chunk number 38: ArrayExpress3 (eval = FALSE)
###################################################
 EGEOD1724 <- getAE("E-GEOD-1724", type='processed')
 cnames = getcolproc(EGEOD1724) # annotation
 EGEOD1724.da <- procset(EGEOD1724, cnames[2]) # build expression set


###################################################
### code chunk number 39: ArrayExpress4
###################################################
EGEOD1724.da


###################################################
### code chunk number 40: IF example
###################################################
i <- 1
if(i == 1){
  	print("i is equal 1")
} else{
	print("i is NOT equal to 1")
}


###################################################
### code chunk number 41: ifelse
###################################################
ifelse(i == 1, "i is equal 1", "i is NOT equal to 1")


###################################################
### code chunk number 42: FOR example
###################################################
for(i in 1:5){
	# do something
	print(i)
}


###################################################
### code chunk number 43: apply
###################################################
mat <- matrix(1:10, nrow=2, byrow=T)
mat
# SUMMING THE COLUMNS
apply(mat, 2, sum)
# SUMMING THE ROWS
apply(mat, 1, sum)
rowSums(mat) # the same


###################################################
### code chunk number 44: foreach example
###################################################
library(foreach)
library(doMC)
library(multicore)
ncore = multicore:::detectCores()
registerDoMC(cores = ncore)
results <- foreach(i = 1:5, .combine=c) %dopar% {
	i+i
}
results


###################################################
### code chunk number 45: apply2
###################################################
library(multicore)
library(rbenchmark)
n <- rep(100, 100)
benchmark(
    x <- mclapply(n, rnorm, mc.cores=ncore),
    x <- lapply(n, rnorm),
    columns = c("test", "replications", "elapsed", "relative"),
    order = "relative",
    replications = 20
)


###################################################
### code chunk number 46: simple function
###################################################
Mr.Euclide <- function(x, y){
	dist <- sqrt(sum((x - y)^2))
	return(dist)
}


###################################################
### code chunk number 47: Euclide distance
###################################################
x <- c(1, 1)
y <- c(2, 2)
Mr.Euclide(x, y)


###################################################
### code chunk number 48: OOeuclide
###################################################
setClass("myObject", representation(vec = "numeric"), 
    prototype = prototype(vec = 0))
(vectorA <- new("amia", vec=c(1, 1, 1)))


###################################################
### code chunk number 49: OOeuclide1
###################################################
Mr.Euclide <- function(x, y){
    if(class(x)!="myObject" & class(y)!="myObject") stop("error") 
    dist <- sqrt(sum((x@vec - y@vec)^2))
    return(dist)
}


###################################################
### code chunk number 50: OOeuclide2
###################################################
vectorB <- new("myObject", vec=c(2, 2, 2))
Mr.Euclide(vectorA, vectorB)
library(bioDist)
euc(matrix(c(vectorA@vec, vectorB@vec), nrow=2, byrow=T))

