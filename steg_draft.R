setwd("/path_to_folder_containing_only_SRGA_files")

library(stringr)

files <- list.files(pattern = "^.*\\.txt$")     ## All .txt files in directory
seq_length <- 79                                ## Enter sequence length
bandwidth <- 5                                  ## Set the bandwidth for the density function
t1 <- 30                                        ## Set lower threshold one 
t2 <- 70                                        ## Set upper threshold two

## Define the data processing function
process <- function(x, y) {
      split <- strsplit(x$V1, split = "")                       ## Convert the reads into columns of bases
      cond <- lapply(split, function (x) length(x) == y)        ## Set condition to filter for reads of length argument "y"
      split2 <- split[unlist(cond)]                             ## Apply the filter condition to the "split" list
      df <- data.frame(matrix(unlist(split2), nrow=length(split2), byrow=T))   ## Convert the list to a data.frame
      dg <- matrix(NA, ncol=length(df[1,]), nrow=5)             ## Create an empty matrix with rows=5 and col=no of col in df
      dimnames(dg)[[1]] <- c("A","T","C","G","-")               ## Name the columns of the empty matrix
      for (i in 1:length(df[1,])){                              ## For loop to set factor levels for all columns
        for(j in  c("A","T","C","G","-")){
          dg[j, i] <- sum(df[,i]==j)   
        }
      }
      dg <- dg[-5, ]                                            ## Remove row containing counts of deletions
      dg_perc <- dg                                             ## Create dummy data.frame for accepting % base at position
      sapply(1:length(dg[1,]), FUN=function(i){ dg_perc[,i] <<- 100*dg[,i]/sum(dg[,i])})  ## Convert base counts matrix into % base at position
}

#########Process data#############################

data_list <- list()

data_list <- lapply(files, function(x) {
  data <- read.table(x, header = FALSE, stringsAsFactors=FALSE)   ## Read in the data from Dario Beraldi's Short Reference Global Aligner
  clean_data <- data[(data$V3 == max(data$V3)), ]                 ## Subset the sequences with top score

  plus_data <- clean_data[clean_data$V4 == "+", ]                 ## Contains all reads from plus strand
  minus_data <- clean_data[clean_data$V4 == "-", ]                ## Contains all reads from minus strand

  processed_plus <- process(plus_data, seq_length)                ## Use the process function on the plus strands
  processed_minus <- process(minus_data, seq_length)              ## Use the process funtion on the minus strands  
  averaged_data <- (processed_plus + processed_minus) / 2         ## Average output of plus and minus data (removes problems of strand imbalance)
})

result <- list()                   ## Create empty integer vector to hold binary output

## Assigning base by % at position. NB: Y = C/T (pyrimidine), R = A/G (purine)
for (j in seq(length(data_list))) {
  result[[j]] <- character()
  for (i in 1:length(data_list[[j]][1, ])){                              
    if (data_list[[j]]["A", i] >= t2) {result[[j]][i] <- "A"}
    else if (data_list[[j]]["T", i] >= t2) {result[[j]][i] <- "T"}
    else if (data_list[[j]]["C", i] >= t2) {result[[j]][i] <- "C"}
    else if (data_list[[j]]["G", i] >= t2) {result[[j]][i] <- "G"}
    else if ( (data_list[[j]]["C", i] > t1) & (data_list[[j]]["C", i] < t2) 
      & (data_list[[j]]["T", i] > t1) & (data_list[[j]]["T", i] < t2) ) {result[[j]][i] <- "Y"}
    else if ( (data_list[[j]]["G", i] > t1) & (data_list[[j]]["G", i] < t2) 
      & (data_list[[j]]["A", i] > t1) & (data_list[[j]]["A", i] < t2) ) {result[[j]][i] <- "R"}  
    else {result[[j]][i] <- "N"}
  }
}

############## Create fasta formatted file ################

export_draft_sequence <- lapply(result, paste0, collapse = "")  

drafts <- list()

for (i in seq(length(export_draft_sequence))) {
  drafts[[i]] <- paste0(">", i)
  drafts[[i]][2] <- export_draft_sequence[[i]]
  write.table(drafts[[i]], paste0("ref_", i, ".txt", collapse = ""), quote = FALSE, eol = "\n", row.names = FALSE, col.names = FALSE)  
}


############## Plots #######################################

for (i in seq(length(data_list))) {
    pdf(file = paste0("ref_", i, ".pdf", collapse = ""), width = 14, pointsize = 16)
    
    ## Set the plot layout
    layout(matrix(c(1,1,2), 3, 1, byrow=TRUE))
    
    ## Plot graphical representation of sequence (dg_perc matrix as stacked barplot)
    barplot(as.matrix(data_list[[i]]), col = c("blue2", "red", "chartreuse3","gold"), 
    xlab = "Sequence Position", ylab = "% base", axes = TRUE, xlim = c(4, ncol(data_list[[i]]) + 18))  
    
    ## Add legend
    legend(58, 145, horiz = TRUE, fill = c("blue2","red","green3","yellow2"), 
    legend = c("A", "T", "C","G"), bty = "n", xpd=TRUE)   
    
    ## Create an empty plot and paste the binary at the top (side 3) and the ASCII at the bottom (side 1)
    par(mar=c(3,3,3,3))
    plot.new()
    mtext(paste0(export_draft_sequence[[i]], collapse=""), side = 3, cex = 0.80)

    dev.off()
}