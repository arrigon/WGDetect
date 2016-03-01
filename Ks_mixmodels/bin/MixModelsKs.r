### Running mixture models from Ks files
# Nils Arrigo, Unil 2015.
########################################################

## UNCOMMENT if willing to pass args from bash command line:
#Import params from bash command line
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }
# 
# # Typical bash command line to invoke R and giving arguments:
# # R CMD BATCH '--args infile="final_ks_values_cat_67.scafSeq" outfolder="."' MixModelsKs.r &
# # You can insert that command in a perl script, e.g. to go parallel in running the thing.
# # Check out RunMixModels.pl to have an idea.


## Script params COMMENT the arguments that are passed from bash
# # I/O (leave as commented if using RunMixModels.pl to run the pipeline)
# infile = "../data.in/final_ks_values_cat_67.scafSeq"
# outfolder = getwd() #by default, output goes in getwd(); set full path to redirect output where you want
# 
# # KS window
# ksmin = 1e-4 #min ks
# ksmax = 2 #max ks
# 
# # Mixture model params
# kmax = 2 #max number of peaks being expected, WARNING: analysis time increases with k
# boots = 5 #bootstrapping effort during search for optimal number of peaks, WARNING: this is time consuming. Advised value is 1000
# epsilon = 1e-3 #convergence criterion; heuristics are stopped when loglik is improved by less than epsilon
# 
# # Graph params
# breaks = 50 #number of breaks on histogram
# 


###################################
########### Analysis; do not modify
## Prepare file names
bsn = gsub("\\..*$", "", basename(infile))
pdfname = paste(outfolder, "/", bsn, ".pdf", sep = "")
outname = paste(outfolder, "/", bsn, ".mixmodels.txt", sep = "")
binname = paste(outfolder, "/", bsn, ".Rdata", sep = "")

# check param values
#I/O
infile
outfolder
pdfname
outname

# Ks window
ksmin
ksmax

# Mixture model params
kmax
boots
epsilon

# Graph params
breaks

## require packages
require(mixtools)
source("bin/plot.mixEM.r")
# if not yet installed: install.packages("mixtools")
# check out documentation here: 
# http://exploringdatablog.blogspot.ch/2011/08/fitting-mixture-distributions-with-r.html
# http://cran.r-project.org/web/packages/mixtools/vignettes/mixtools.pdf


## Import and format data
data = read.delim(infile, header = F, as.is = T)
dat = data.frame(No = 1:nrow(data), Node.Ks = as.numeric(data[, 1]))

#discard nodes with missing Ks values
dat = dat[!is.na(dat$Node.Ks), ]

# focus on Ks values ranging between 0.0001 and 2 (this range can be changed)
dat = dat[dat$Node.Ks >= ksmin & dat$Node.Ks <= ksmax, ]   
  
  
## Compute mixture model
# Determine number of peaks, using parametric bootstrapping
a <- boot.comp(y = dat$Node.Ks, 
	      max.comp = kmax, 
	      B = boots, 
	      mix.type = "normalmix", 
	      epsilon = epsilon,
	      hist = F)
k.best = length(a$p.values) # as in boot.comp


# Estimate parameters (mixture of normal distributions)
mix.object = normalmixEM(dat$Node.Ks, 
			  lambda = NULL, #no prior on peak contributions
			  mu = NULL, #no prior on peak centers
			  sigma = NULL, #no prior on peak stdevs
			  k = k.best, #prior on number of peaks
			  epsilon = epsilon) #convergence param (changes in loglik among iterations)

centers = mix.object$mu
stdevs = mix.object$sigma
contribs = mix.object$lambda
loglik = mix.object$loglik
params = data.frame(k.best, centers, stdevs, contribs, loglik)


## Produce graph
pdf(file = pdfname)
  graph.params = hist(mix.object$x, plot = F)
  plot.mixEM(mix.object, which = 2, breaks = breaks, main2 = infile, xlab2 = "Ks values", ylab2 = "Number of duplicates")
  text(x = params$centers, 
       y = max(graph.params$density), 
       labels = paste(round(params$centers, 2), "+/-", round(params$stdevs, 2), "SD"), 
       col = 2:(nrow(params)+1),
       srt = 90)
  dev.off()

  
## Save outputs
write.table(params, file = outname, sep = "\t", quote = F, col.names = T, row.names = F)
#save(objects(), file = binname)


