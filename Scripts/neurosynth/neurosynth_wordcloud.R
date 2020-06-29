library(wordcloud)
library(scales)
# load in words for cloud
# column 1: words, column 2: weights
base.path = '~/Documents/Education/Cambridge/fMRI_NSPN/'
words <- read.csv(paste0(base.path, 'Scripts/prep_neurosynth/analyse_terms/MI_terms_neurosynth_no_atomical.csv',header=FALSE))
nWords = 50 # choose how many words to display

words.pos = words[words$V2>0,] # choosing only positive weights
freq.pos = words.pos$V2/sum(words.pos$V2)*100 # transform into frequency
pdf(paste0(base.path, '/Results/maturational_index/plots/wordcloud_pos.pdf')) # save image as pdf file
wordcloud(words=words.pos$V1, 
          freq=freq.pos, 
          max.words = nWords, 
          min.freq = min(freq.pos), 
          random.order = FALSE,
          scale=c(3,.2), 
          colors = '#C90C0C')
dev.off()

words.neg = words[words$V2<0,]
words.neg = words.neg[((dim(words.neg)[1]-nWords):dim(words.neg)[1]),]
freq.neg = abs(words.neg$V2)/sum(abs(words.neg$V2))*100
pdf(paste0(base.path, '/Results/maturational_index/plots/wordcloud_neg.pdf')) # save image as pdf file
wordcloud(words=words.neg$V1, 
          freq=freq.neg, 
          max.words = nWords, 
          random.order = FALSE, 
          min.freq = min(freq.neg), 
          scale=c(2.5,0.05), 
          colors = '#0C13C9',
          family='Helvetica') 
dev.off()


