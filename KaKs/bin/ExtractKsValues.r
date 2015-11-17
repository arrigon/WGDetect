### Extracting Ks values from tree
# Nils Arrigo, Unil 2015.
########################################################

#Import params from bash command line
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }
  
  
# list files to analyse
listfiles = list.files(infolder, pattern = "2NG.dS", full.names = T, recursive = T)


# function to convert phylip file into distance matrix
splitfun = function(x) strsplit(x, split = "\\s+")[[1]]


# loop over files and extract Ks values
cat("", file = outfile, append = F) #start fresh
Ks = NULL
for(infile in listfiles){
  # load data
  data = readLines(infile)

  if(length(data) == 0) next()
  
  # convert to dist
  ntax = as.numeric(splitfun(data[1])[2]) + 1
  MAT = matrix(NA, ntax, ntax)
  for(i in 2:length(data)){
    x = rep(NA, ntax)
    tmp = splitfun(data[i])
    x[1:length(tmp)] = tmp
    MAT[i, ] = x
    }
    
  # compute quick WGMA tree and extract Ks values of nodes
  rownames(MAT) = colnames(MAT) = MAT[, 1]
  MAT = MAT[-1, -1]
  MAT = as.data.frame(MAT, stringsAsFactors = F)
  MAT = as.dist(MAT)
  tre = hclust(MAT, method = "average")

  # extract node heights and save to output
  heights = tre$height
  heights = heights[heights >= 0]
  if(length(heights) > 0){
    write.table(cbind(heights), file = outfile, append = T, sep = "\t", quote = F, row.names = F, col.names = F)
    Ks = c(Ks, heights)
    }
  }

