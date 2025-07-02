# Check if required packages are already installed, and install if missing
#install.packages('/Users/tom/Downloads/vegan_1.6-0.tar.gz', repos = NULL, type="source")
#update.packages(ask=FALSE, checkBuilt=TRUE,repos = "http://cran.us.r-project.org")
#withr::with_makevars(c(PKG_LIBS = "-liconv"), install.packages("vegan",repos = "http://cran.us.r-project.org"), assignment = "+=")


library("vegan")
#library("dyld")


args = commandArgs(trailingOnly=TRUE)

  
print (args)

getwd()
setwd(getwd())
#print(paste0("Current working dir: ", getwd()))
  


maldi_data = args[1]

path_wanted <- dirname(normalizePath(maldi_data))

raw.data <- read.table (file = maldi_data , check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")

data <- t(raw.data)


d_Binary_matrix<- as.matrix(vegdist(data, method='jaccard'))

write.table(d_Binary_matrix, file ="Sample-Consortia-distances.tsv", sep = "\t")
