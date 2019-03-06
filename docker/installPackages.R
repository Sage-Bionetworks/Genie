install.packages("synapser", repos=c("https://sage-bionetworks.github.io/ran", "http://cran.fhcrc.org"))
install.packages("dplyr", repo="https://ftp.osuosl.org/pub/cran/")
install.packages("argparse", repo="https://ftp.osuosl.org/pub/cran/")
install.packages("rmarkdown", repo="https://ftp.osuosl.org/pub/cran/")
install.packages("UpSetR", repo="https://ftp.osuosl.org/pub/cran/")

source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")

library(synapser)
library(dplyr)
library(argparse)
library(UpSetR)
library(VariantAnnotation)