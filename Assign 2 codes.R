#script for taxonomic classification using randomForest


#load libraries
library(rentrez)
library(seqinr)
library(Biostrings)
library(tidyverse)
library(randomForest)
library(ggplot2)
library(cvms)

--------####Part 1------------

#search database ("nuccore") for 16S rRNA of Bifidobacterium and Collinsella
#Both are gram-positive bacteria with a high G + C content

ncbi_search_result <- entrez_search(
  db = "nuccore",
  term = "Bifidobacterium OR Collinsella AND 16S rRNA AND 500:1000[SLEN])",
  use_history=TRUE)

#search return count
ncbi_search_result$count

mxhit <- 2000

#second search using maximum search return count
ncbi_search_result2 <- entrez_search(db = "nuccore",
  term = "Bifidobacterium OR Collinsella  AND 16S rRNA AND 500:1000[SLEN])",
  retmax = mxhit, use_history=TRUE)

#search return count
ncbi_search_result2$retmax

#split search_result into groups of 200 hits
searchsplit <- split(seq(1,length(ncbi_search_result2$ids)),
                     ceiling(seq_along(seq(1,length(ncbi_search_result2$ids)))/200))


ncbi_summ <- list() ##empty list to hold entrez_summary contents

# loop through the list of ids in "searchsplit" to obtain summary, and pause for 5 seconds in between splits
for (i in searchsplit){
  ncbi_summ[unlist(i)] <- entrez_summary("nuccore", id = ncbi_search_result2$ids[unlist(i)])
  Sys.sleep(5)
}

#check names of organisms using extract elements from summary
organism <- sapply(ncbi_summ, function(x) extract_from_esummary(x, elements = "organism"))

unique(organism)


----------####Part 2---------------

#fetch fasta files

#search result ids
ncbi_ids <- ncbi_search_result2$ids

#obtain fasta files using entrez_fetch
ncbi_fetch1 <- entrez_fetch(db="nuccore", id=ncbi_ids[1:300], rettype="fasta")
ncbi_fetch2 <- entrez_fetch(db="nuccore", id=ncbi_ids[301:600], rettype="fasta")
ncbi_fetch3 <- entrez_fetch(db="nuccore", id=ncbi_ids[601:900], rettype="fasta")
ncbi_fetch4 <- entrez_fetch(db="nuccore", id=ncbi_ids[901:1200], rettype="fasta")
ncbi_fetch5 <- entrez_fetch(db="nuccore", id=ncbi_ids[1201:1500], rettype="fasta")
ncbi_fetch6 <- entrez_fetch(db="nuccore", id=ncbi_ids[1501:1800], rettype="fasta")
ncbi_fetch7 <- entrez_fetch(db="nuccore", id=ncbi_ids[1801:2000], rettype="fasta")

#check fasta file
class(ncbi_fetch1)
head(ncbi_fetch1)

# Write to file.
write(ncbi_fetch1, "ncbi1.fasta", sep = "\n")
write(ncbi_fetch2, "ncbi2.fasta", sep = "\n")
write(ncbi_fetch3, "ncbi3.fasta", sep = "\n")
write(ncbi_fetch4, "ncbi4.fasta", sep = "\n")
write(ncbi_fetch5, "ncbi5.fasta", sep = "\n")
write(ncbi_fetch6, "ncbi6.fasta", sep = "\n")
write(ncbi_fetch7, "ncbi7.fasta", sep = "\n")

#Read in fasta files
fastaFiles <- list.files(pattern = "ncbi")
l_fastaFiles <- lapply(fastaFiles, readDNAStringSet)

#Convert to dataframe
l_dfFastaFiles <- lapply(
  l_fastaFiles, function(x) data.frame(Title = names(x), Sequence = paste(x) ))

# Combine the list of dataframes into one dataframe.
dfSeq <- do.call("rbind", l_dfFastaFiles)


#rm(l_fastaFiles, l_dfFastaFiles, fastaFiles)
#rm(bif_fetch1, bif_fetch2, bif_ids)

#length of unique names in the dataframe
length(unique(dfSeq$Title))
unique(dfSeq$Title)


----####Edit1####---------------------------------------------------------------
# I would recommend using the helper functions "Entrez_Functions.R" to fetch 
# data from online database. These functions are provided on Courselink.
# Edit 1 replaces codes from line 12 - 107

ncbi_search_result <- entrez_search(db = "nuccore",
  term = "Bifidobacterium OR Collinsella AND 16S rRNA AND 500:1000[SLEN])",
  use_history=TRUE)

#search return count
ncbi_search_result$count

mxhit <- 2000

# If you would like to specify your own maximum number of hits, I recommend 
# changing it in "Entrez_Functions.R". Then upload the modified version to your 
# own repository.

source("Entrez_Functions.R")

# fetch fasta files from NCBI
FetchFastaFiles(
  searchTerm = "Bifidobacterium OR Collinsella  AND 16S rRNA AND 500:1000[SLEN]", 
  seqsPerFile = 100, fastaFileName = "entrez_fetch.fasta")

# merge all fasta files and store it as a data frame called dfSeq
dfSeq <- MergeFastaFiles(filePattern = "entrez_fetch.fasta*")
--------------------------------------------------------------------------------
  

---------####Part 3-----------

#cleaning the dataframe

dfSeq$genus_name <- word(dfSeq$Title, 2L)

# Rearrange the columns.
dfSeq <- dfSeq[, c("Title", "genus_name", "Sequence")]

head(dfSeq)
unique(dfSeq$genus_name)
unique(dfSeq$Title)
count(dfSeq, vars = genus_name)

##filter
dfSeq2 <- dfSeq %>%
  filter(genus_name %in% c("Bifidobacterium", "Collinsella"))

#check
unique(dfSeq2$genus_name)
head(dfSeq2)


---------####Part 3-------

#checking sequence length of genus

#summary of sequence length
summary(nchar(dfSeq2$Sequence[dfSeq2$genus_name == "Bifidobacterium"]))
summary(nchar(dfSeq2$Sequence[dfSeq2$genus_name == "Collinsella"]))


#filter sequences
df16s <- dfSeq2 %>%
  filter(!is.na(Sequence)) %>%
  mutate(Sequence2 = Sequence) %>%
  mutate(Sequence2 = str_remove(Sequence2, "^[-N]+")) %>%
  mutate(Sequence2 = str_remove(Sequence2, "[-N]+$")) %>%
  mutate(Sequence2 = str_remove_all(Sequence2, "-+")) %>%
  filter(str_count(Sequence2, "N") <= (0.02 * str_count(Sequence2)))

#check
sum(is.na(df16s$Sequence2))
unique(df16s$genus_name)
count(df16s, vars = genus_name)


#summary of sequence length after removing "N"s
summary(nchar(df16s$Sequence2))


#calculating first quartile length value and third quartile value.

q1 <- quantile(nchar(df16s$Sequence2), probs = 0.25, na.rm = TRUE)
q1

q3 <- quantile(nchar(df16s$Sequence2), probs = 0.75, na.rm = TRUE)
q3

#do not want to remove Collinsella sequence length less than q1 to avoid having low number
#so, firstly filtering Bifidobacterium to keep sequence length greater than or equal to q1 and those less than or equal to q3, and all Collinsella

df16s <- df16s %>% filter((
  str_count(Sequence2) >= q1 & str_count(Sequence2) <= q3 & 
    genus_name == "Bifidobacterium") | genus_name == "Collinsella")


#secondly, keeping Collinsella sequence length less than or equal to q3, and all Bifidobacterium

df16s <- df16s %>% filter((str_count(Sequence2) <= q3 & 
                             genus_name == "Collinsella") | 
                             genus_name == "Bifidobacterium")

#check count of each genus
count(df16s, vars = genus_name)

#histogram and summary of sequence lengths

pdf(file = "./My Assign2Plot.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5) # The height of the plot in inches

hist(nchar(df16s$Sequence2[df16s$genus_name == "Bifidobacterium"]),
     xlab = "Sequence Length",
     ylab = "Frequency",
     main = "Frequency Histogram of Bifidobacterium Sequence Lengths")


hist(nchar(df16s$Sequence2[df16s$genus_name == "Collinsella"]),
     xlab = "Sequence Length",
     ylab = "Frequency",
     main = "Frequency Histogram of Collinsella Sequence Lengths")

dev.off()

summary(nchar(df16s$Sequence2[df16s$genus_name == "Collinsella"]))
summary(nchar(df16s$Sequence2[df16s$genus_name == "Bifidobacterium"]))

----####Edit2####--------------------------------------------------------------
# I would recommend creating histograms using ggplot2 instead of base R. 
# Consider using R Markdown to add figure captions.
# Edit 2 replaces codes from line 216 - 225

# Figure 1.Histogram Illustrating the Distribution of Bifidobacterium 16s Sequence Length
df16s %>%
  filter(genus_name == "Bifidobacterium") %>%
  ggplot(aes(x = nchar(Sequence2))) +
  geom_histogram(bins = 15, fill = "mediumpurple", color = "black") +
  ggtitle("Frequency Distribution of Bifidobacterium 16s rRNA Sequence Lengths") +
  xlab('Sequence Length') +
  ylab('Frequency') +
  # solid vertical lines represent the median
  geom_vline(aes(xintercept = median(nchar(Sequence2))), 
             color = "orange", linewidth = 1) +
  # dashed lines represent the 1st & 3rd quartiles
  geom_vline(aes(xintercept = quantile(nchar(Sequence2), probs = c(0.25))), 
             color = "orange", linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(nchar(Sequence2), probs = c(0.75))), 
             color = "orange", linewidth = 1, linetype = "dashed") +
  theme(plot.title = element_text(size = 15))


# Figure 2.Histogram Illustrating the Distribution of Collinsella 16s Sequence Length
df16s %>%
  filter(genus_name == "Collinsella") %>%
  ggplot(aes(x = nchar(Sequence2))) +
  geom_histogram(bins = 15, fill = "mediumpurple", color = "black") +
  ggtitle("Frequency Distribution of Collinsella 16s rRNA Sequence Lengths") +
  xlab('Sequence Length') +
  ylab('Frequency') +
  # solid vertical lines represent the median
  geom_vline(aes(xintercept = median(nchar(Sequence2))), 
             color = "orange", linewidth = 1) +
  # dashed lines represent the 1st & 3rd quartiles
  geom_vline(aes(xintercept = quantile(nchar(Sequence2), probs = c(0.25))), 
             color = "orange", linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(nchar(Sequence2), probs = c(0.75))), 
             color = "orange", linewidth = 1, linetype = "dashed") +
  theme(plot.title = element_text(size = 15))
--------------------------------------------------------------------------------


---------####Part 4---------

#Calculating sequence features and building a classifier 

#convert class
df16s <- as.data.frame(df16s)
df16s$Sequence2 <- DNAStringSet(df16s$Sequence2)

#calculate nucleotide frequencies. 
df16s <- cbind(df16s, as.data.frame(letterFrequency(df16s$Sequence2, letters = c("A", "C","G", "T"))))

#Adding A, T, and G proportions as new columns.
df16s$Aprop <- (df16s$A) / (df16s$A + df16s$T + df16s$C + df16s$G)

df16s$Tprop <- (df16s$T) / (df16s$A + df16s$T + df16s$C + df16s$G)

df16s$Gprop <- (df16s$G) / (df16s$A + df16s$T + df16s$C + df16s$G)

#Adding dinucleotide frequency
df16s <- cbind(df16s, as.data.frame(dinucleotideFrequency(df16s$Sequence2, as.prob = TRUE)))

#Adding oligonucleotide frequencies 3 and 4
df16s <- cbind(df16s, 
               as.data.frame(oligonucleotideFrequency(df16s$Sequence2, 
                                                      width = 3,
                                                      as.prob = TRUE)))

df16s <- cbind(df16s, 
               as.data.frame(oligonucleotideFrequency(df16s$Sequence2, 
                                                      width = 4,
                                                      as.prob = TRUE)))
df16s <- cbind(df16s, 
               as.data.frame(oligonucleotideFrequency(df16s$Sequence2, 
                                                      width = 5,
                                                      as.prob = TRUE)))
df16s <- cbind(df16s, 
               as.data.frame(oligonucleotideFrequency(df16s$Sequence2, 
                                                      width = 6,
                                                      as.prob = TRUE)))

view(df16s)

#converting string format back to character data
df16s$Sequence2 <- as.character(df16s$Sequence2)

#counts by subspecie
table(df16s$genus_name)

#set variable for sample size (based on the smaller Species_Name)
smaller_sample <- min(table(df16s$genus_name))
smaller_sample


#creating a Validation dataset
set.seed(28)

df16sValidate <- df16s %>%
  group_by(genus_name) %>%
  sample_n(floor(0.2 * smaller_sample))

#Checking
table(df16sValidate$genus_name)

#creating a training dataset
set.seed(765)

df16sTraining <- df16s %>%
  filter(!Title %in% df16sValidate$Title) %>%
  group_by(genus_name) %>%
  sample_n(ceiling(0.8 * smaller_sample))

names(df16s)


#Building classifiers

#first classifier using k-mer 1 nucleotide proportions

classifier1 <- randomForest(x = df16sTraining[, 9:11], 
                            y = as.factor(df16sTraining$genus_name), 
                            ntree = 40, importance = TRUE)

#check model performance on the validation data.
table(df16sValidate$genus_name)

predict16sValidate1 <- predict(classifier1, df16sValidate[, c(92, 9:11)])

#confusion matrix.
table(observed = df16sValidate$genus_name, predicted = predict16sValidate1)


#second classifier using k-mer 2 nucleotide proportions

classifier2 <- randomForest(x = df16sTraining[, 12:27], 
                            y = as.factor(df16sTraining$genus_name), 
                            ntree = 50, importance = TRUE)

predict16sValidate2 <- predict(classifier2, df16sValidate[, c(92, 12:27)])

#confusion matrix.
table(observed = df16sValidate$genus_name, predicted = predict16sValidate2)


#Third classifier using k-mer 3 nucleotide proportions

classifier3 <- randomForest(x = df16sTraining[, 28:91], 
                            y = as.factor(df16sTraining$genus_name), 
                            ntree = 50, importance = TRUE)

predict16sValidate3 <- predict(classifier3, df16sValidate[, c(92, 28:91)])

#confusion matrix.
table(observed = df16sValidate$genus_name, predicted = predict16sValidate3)


#Fourth classifier using k-mer 4 nucleotide proportions

classifier4 <- randomForest(x = df16sTraining[, 92:347], 
                            y = as.factor(df16sTraining$genus_name), 
                            ntree = 100, importance = TRUE)

predict16sValidate4 <- predict(classifier4, df16sValidate[, c(92, 92:347)])

#confusion matrix.
table(observed = df16sValidate$genus_name, predicted = predict16sValidate4)

#check accuracy by viewing OOB estimate of error rate

print(classifier1)
print(classifier4)


----####Edit3####---------------------------------------------------------------
# I would recommend defining functions to reduce code redundancy and increase
# code reusability
# Edit 3 replaces code from line 350 - 400.

# 1. Define function to build random forest classifiers using the training dataset
build_classifier <- function(dataset, predictor_indices, variable, num_trees) {
  classifier <- randomForest::randomForest(
    x = dataset[, predictor_indices],
    y = as.factor(dataset$variable),
    ntree = num_trees,
    importance = TRUE
  )
  return(classifier)
}


# 2. Define function to test the prediction accuracy of classifiers using 
# validation dataset

validate_classifier <- function(classifier, dataset, predictor_indices) {
  predictions <- predict(classifier, dataset[, predictor_indices])
  return(predictions)
}


# 3. Define function to create confusion matrix plots
confusion_matrix_plot <- function(observed, predicted) {
  # Create confusion matrix
  matrix_table <- table(observed = observed, predicted = predicted)
  # Store the table as a tibble
  matrix_tibble <- as_tibble(matrix_table)
  # Plot the tibble as a confusion matrix plot
  plot <- plot_confusion_matrix(matrix_tibble, target_col = "observed", 
                                prediction_col = "predicted", counts_col = "n", 
                                add_sums = TRUE, add_normalized = FALSE,
                                add_col_percentages = FALSE,
                                add_row_percentages = FALSE, palette = "Purples", 
                                sums_settings = sum_tile_settings(palette = "Oranges", 
                                                                  label = "Total"))
  return(plot)
}

# Build, validate, and plot random forest classifiers using the functions above
classifier1 <- build_classifier(df16sTraining, c(92, 9:11), "genus_name", 40)
valication1 <- validate_classifier(classifier1, df16sValidate, c(92, 9:11))
plot1 < confusion_matrix_plot(df16sValidate$genus_name, valication1)

  
  
