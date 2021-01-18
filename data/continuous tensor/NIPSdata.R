
####################################################
# data preprocessing ###############################
papers = read.table("~/mode-1-papers.map")
papers = papers$V1
authors = read.table("~/mode-2-authors.map")
authors = authors$V1
words = read.table("~/mode-3-words-2.map")
words = words$V1
years = read.table("~/mode-4-years.map")
years = years$V1
length(papers) #2482
length(authors) #2862
length(words) #14036
length(years) #17


tns  = read.table("~/nips.tns")
str(tns)

# Authors and papers are corresponding each other (one variable is enough)
# Merge authors because multiple authors make words counted multiple times.

ntns = aggregate(tns[2],tns[-2],unique) 
# V1: papers, V2: words, V3: years, V4: counts, 
NIPS = ntns[1:4]
names(NIPS) = c("papers","words","years","counts")
range(NIPS$papers)
range(NIPS$words)
range(NIPS$years)

# reshape dataframe to tensor
library(reshape2)
tnsNIPS = acast(NIPS,papers~words~years,value.var = "counts")

save(tnsNIPS,NIPS,papers,words,years,file = "NIPS.RData")

 
#################################################################
### getting smaller size tensor (2482,14036,17) =>(500,510,17)

# want to reduce the size

############################
# reducing # of words to 500

par(mfrow = c(2,1))
num_w = aggregate(NIPS[4],NIPS[2],sum)
hist(num_w$counts,breaks = 100,xlab = "# of words",main = "words (full)")
range(num_w$counts) #6~23921
words[order(num_w$counts,decreasing = T)[1:500]]
rwindex = sort(order(num_w$counts,decreasing = T)[1:500])
rwords = words[rwindex]

hist(num_w[rwindex,]$counts,breaks = 100,xlab = "# of words",main = "words (reduced)")
range(num_w[rwindex,]$counts) #1373~23921

#############################
# reducing # of papers to 510
num_w_in_p = aggregate(NIPS[4],NIPS[1],sum)
hist(num_w_in_p$counts,breaks = 100,xlab = "# of words in a paper",main = "paper (full)",
     xlim = range(num_w_in_p$counts))

p_y = NULL
for(i in 1:17){
  p_y = c(p_y,length(unique(NIPS[NIPS$year==i,1])))
}

# randomly sample 30 papers for each year
a = cumsum(p_y)
rpindex = NULL
s = 0
for(i in a){
  s = s+1
  rpindex = c(rpindex,sort(sample((c(0,a)[s]+1):c(0,a)[s+1],30)))
}

length(rpindex)
rpapers = papers[rpindex]
hist(num_w_in_p[rpindex,]$counts,breaks = 100,xlab = "# of words in a paper",
     main = "paper (reduced)",xlim = range(num_w_in_p$counts))

## Reducing tensor
rtnsNIPS = tnsNIPS[rpindex,rwindex,]
dimnames(rtnsNIPS) = NULL
rwords
rpapers
ryears = years

save(rtnsNIPS,rwords,rpapers,ryears,file = "rNIPS.RData")
