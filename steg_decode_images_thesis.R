## R script containing functions and library dependencies
suppressWarnings(source("STEG_functions.R"))

prefix <- format(Sys.time(), "%Y%m%d_%H%M")

parser <- ArgumentParser(description = "Settings for STEG_PORTRAITS")

parser$add_argument("--length", "-l",
                    type = "integer",
                    help = "Define the length of sequences")
parser$add_argument("--qualitycut", "-q",
                    type = "double",
                    default = 0.2,
                    help = "the quantile below which data is excluded")
parser$add_argument("--messagestart", "-s",
                    type = "integer",
                    default = 1,
                    help = "the first sequence position included in the message")
parser$add_argument("--messageend", "-e",
                    type = "integer",
                    help = "the final sequence position included in the message")
parser$add_argument("--lowerthreshold", "-t1",
                    type = "integer",
                    default = 30,
                    help = "the lower threshold for assigning bits")
parser$add_argument("--upperthreshold", "-t2",
                    type = "integer",
                    default = 70,
                    help = "the upper threshold for assigning bits")
parser$add_argument("--seq", 
                    action = "store_true",
                    help = "run script in standard mode")                  
parser$add_argument("--BS", 
                    action = "store_true",
                    help = "run script in BS mode") 
parser$add_argument("--oxBS", 
                    action = "store_true",
                    help = "run script in oxBS mode")
parser$add_argument("--mclist", "-mc", 
                    default = "empty",
                    type = "character",
                    help = "filename of list containing 5mC positions (generated from --oxBS)") 
parser$add_argument("--binaryoutput", "-bo", 
                    type = "character",
                    default = "bin",
                    help = "filename of the binary output") 
parser$add_argument("--asciioutput", "-ao", 
                    type = "character",
                    default = "ascii",
                    help = "filename of the ascii output") 
parser$add_argument("--graphicaloutput", "-go", 
                    type = "character",
                    default = "",
                    help = "filename of the graphical output") 
parser$add_argument("--portraitoutput", "-po", 
                    type = "character",
                    default = "portrait",
                    help = "filename of the graphical output")
parser$add_argument("--codestop", "-stop", 
                    default = "111",
                    type = "character",
                    help = "gamma elias stop code") 

args <- parser$parse_args()

if (is.null(args$length)) {
  cat("** No sequence length entered **\n")
  cat("** Exiting script **\n")
  quit(save = "no", status = 1)
}
if(args$qualitycut > 1) {
  cat("** Invalid quality cut chosen **\n")
  cat("** Exiting script **\n")
  quit(save = "no", status = 1)
}
if (is.null(args$messageend)) {
  cat("** No sequence position given for end of message **\n")
  cat("** Exiting script **\n")
  quit(save = "no", status = 1)
}
if (args$seq == FALSE & args$BS == FALSE & args$oxBS == FALSE) {
  cat("** No run mode selected **\n")
  cat("** Exiting script **\n")
  quit(save = "no", status = 1)
}
if ((args$BS + args$oxBS + args$seq) > 1) {
  cat("** Multiple run modes selected. Please select one run mode **\n")
  cat("** Exiting script **\n")
  quit(save = "no", status = 1)
}
if (args$BS == TRUE & args$mclist == "empty") {
  cat("** BS mode selected but no 5mC position list provided **\n")
  cat("** Unable to apply cipher **\n\n")
  cat("** Saving all modified cytosine positions**\n")
}

cat("\n** Please check run parameters **\n\n")
print(args)

### PARAMETERS ###

files <- list.files(pattern = "^.*\\.txt$")  ## Must be in working directory containing SRGA files 
seq_length <- args$length                    ## Enter sequence length
quality_cut <- args$qualitycut               ## Enter minimum quality quantile acceptable 
start_mes <- args$messagestart               ## Message start position in sequence
end_mes <- args$messageend                   ## Message end position in sequence
bandwidth <- 5                               ## Set the bandwidth for the density function
t1 <- args$lowerthreshold                    ## Set the lower threshold
t2 <- args$upperthreshold                    ## Set the lower threshold
bin_output <- paste0(prefix, "_", args$binaryoutput, ".out")    ## Enter name for binary output
ascii_output <- paste0(prefix, "_", args$asciioutput, ".out")   ## Enter name for ascii output 
pdf_output <- paste0(prefix, args$graphicaloutput, ".pdf")      ## Enter name for graphical output
portrait_output <- paste0(prefix, args$portraitoutput, ".pdf")  ## Enter name for portrait output


### PROCESS DATA ###

## Empty list to receive aligned data; filled with one oligo per item
data_list <- list()

cat("** PROCESSING DATA **\n")
data_list <- lapply(files, function(x) {
  data <- read.table(x, header = FALSE, stringsAsFactors=FALSE)  ## Read in the data from Dario Beraldi's Short Reference Global Aligner
  plus_data <- data[data$V4 == "+", ]                 ## Contains all reads from plus strand
  minus_data <- data[data$V4 == "-", ]                ## Contains all reads from minus strand
  
  clean_plus_data <- plus_data[(plus_data$V3 >= quantile(plus_data$V3, quality_cut)) & (plus_data$V5 == seq_length), ]      ## Set the alignment quality threshold
  clean_minus_data <- minus_data[(minus_data$V3 >= quantile(minus_data$V3, quality_cut)) & (minus_data$V5 == seq_length), ]      ## Set the alignment quality threshold
  
  processed_plus <- process(clean_plus_data, seq_length)       ## Process function on plus strands
  processed_minus <- process(clean_minus_data, seq_length)     ## Process funtion on minus strands  
  averaged_data <- (processed_plus + processed_minus) / 2      ## Average output of plus and minus data, removes problems from strand imbalance
})

## Isolate message portion, using user provided start and end points
message_portion <- lapply(data_list, `[`, i =, j = c(start_mes:end_mes))

### CONVERSION TO BINARY ###

cat("** CONVERTING TO BINARY **\n")
## Empty list to receive binary assignments; filled with one oligo per item
result <- list()                   
mC_list <- list()

for (j in seq(length(message_portion))) {
  result[[j]] <- numeric()
  mC_list[[j]] <- numeric(end_mes-start_mes+1)
  ## Assigning 0 or 1 based on rules (A & C = 0, T & G = 1, bisulfite G/A = 0, bisulfite C/T = 1)
  for (i in 1:length(message_portion[[j]][1,])){
    if (message_portion[[j]]["A", i] > t2) {result[[j]][i] <- 0}
    else if (message_portion[[j]]["T", i] >= t2) {result[[j]][i] <- 1}
    else if (message_portion[[j]]["C", i] >= t2) {result[[j]][i] <- 0; mC_list[[j]][i] <- 1}
    else if (message_portion[[j]]["G", i] >= t2) {result[[j]][i] <- 1; mC_list[[j]][i] <- 1}
    else if ( (message_portion[[j]]["C", i] > t1) & (message_portion[[j]]["C", i] < t2) & (message_portion[[j]]["T", i] > t1) & (message_portion[[j]]["T", i] < t2) ) {result[[j]][i] <- 1}
    else if ( (message_portion[[j]]["G", i] > t1) & (message_portion[[j]]["G", i] < t2) & (message_portion[[j]]["A", i] > t1) & (message_portion[[j]]["A", i] < t2) ) {result[[j]][i] <- 0}  
  }
}


## Save the mC_list as an R_object to be imported when run in mode == "BS"
## If mode != "oxBS" delete mC_list to ensure cipher cannot be applied
if (args$oxBS == TRUE) {
  save(mC_list, file = paste0(prefix, "_mC_list"))
} else if (args$BS == TRUE & args$mclist == "empty") {
  mods_list <- mC_list
  save(mods_list, file = paste0(prefix, "_modC_list"))
  } else if (args$seq == TRUE) {
    all_C_list <- mC_list
    save(all_C_list, file = paste0(prefix, "_allC_list"))
  } else {
  rm(mC_list)
}

## If mode == "BS" then load the R object 5mC_list created in mode == "oxBS"
## Apply cipher: at 5mC positions (as stored in the 5mC_list variable) the bit value is inverted
if (args$BS == TRUE & args$mclist != "empty") {
  load(args$mclist)
  for (j in seq(length(message_portion))) {
    for (i in 1:length(mC_list[[j]])) {
      if (mC_list[[j]][i] == 1) {
        if (result[[j]][i] == 0) {result[[j]][i] <- 1}
        else if (result[[j]][i] == 1) {result[[j]][i] <- 0}
      }  
    }
  }
}


### TEXT FILE OF BINARY OUTPUT ###

## Prep variable for printing to screen and output to file
export_bin_message <- lapply(result, paste0, collapse = "")   

## Export to .txt file
write.table(export_bin_message, file = bin_output, quote = FALSE, eol = "", 
            row.names = FALSE, col.names = FALSE)  

### DATA PLOTTING ###
cat("** PLOTTING **\n")


## Define plotting colors
colA <- rgb(0,255,255,max = 255)
colT <- rgb(255,0,0,max = 255)
colC <- rgb(255,255,0,max = 255)
colG <- "white"

## Start the graphics driver for producing pdf graphics
pdf(file = pdf_output, width = 12, pointsize = 16, colormodel = "cmyk")

## Set the plot layout
layout(matrix(c(1,1,2,2,3), 5, 1, byrow=TRUE))

## Produce the plots: one page per oligo
## Set up for Raven message
for (i in seq(length(data_list))) {      
  ## Plot graphical representation of sequence
  par(mar=c(5.1,4.1,4.1,2.1))
  barplot(message_portion[[i]], col = c(colA, colT, colC, colG), 
          xpd = TRUE, space = 0.2, border = c(colA, colT, colC, colG),
          xlab = "Position in message", ylab = "% base")
  box()
  
  legend("top", horiz = TRUE, fill = c(colA, colT, colC, colG), 
         legend = c("A", "T", "C","G"), bty = "n", xpd=TRUE,
         inset = -0.22)     ## Add a legend to the barplot
  
  ## Plot the four bases on top of each other in different colors
  
  unlisted_data <- as.vector(message_portion[[i]])
  columns_of_all_A_T_C_G <- matrix(unlisted_data, nrow = ncol(message_portion[[i]]), 
                                   ncol = 4, byrow = TRUE)
  
  #columns_of_all_A_T_C_G[1:5, ]
  vector_A <- columns_of_all_A_T_C_G[ ,1]
  vector_T <- columns_of_all_A_T_C_G[ ,2]
  vector_C <- columns_of_all_A_T_C_G[ ,3]
  vector_G <- columns_of_all_A_T_C_G[ ,4]
  
  plot(seq(length(vector_A)), (vector_A), 
       xlab = "Position in message", ylab = "% base", 
       col = colA, pch = 1, cex = 0.8)
  points(c(1:length(vector_G)), (vector_G), 
         col = "black", pch = 2, cex = 0.8, lwd = 0.75)
  points(c(1:length(vector_C)), (vector_C), 
         col = colC, pch = 3, cex = 0.8)
  points(c(1:length(vector_T)), (vector_T), 
         col = colT, pch = 4, cex = 0.8)
  
  legend("top", horiz = TRUE, fill = c(colA, colT, colC, "black"), 
         legend = c("A", "T", "C","G"), bty = "n", xpd=TRUE,
         inset = -0.22)     ## Add a legend to the barplot
  
  axis(4, at = t1, label = "T1", col = "lightgrey", col.axis = "lightgrey")
  axis(4, at = t2, label = "T2", col = "lightgrey", col.axis = "lightgrey")
  abline(h=c(t1,t2), lty = 2, col = "lightgrey") 
  
  ## Create empty plot and paste binary at the top (side 3) and the ASCII at the bottom (side 1)
  par(mar=c(3,3,3,3))
  plot.new()
  mtext(paste0(result[[i]], collapse=""), side = 3, cex = 0.80)
}

## Close the file connection
invisible(dev.off())

### Decode the elias gamma binary sequence ###
cat("** DECODING ELIAS**\n")

## Paste into a continuous sequence
decode1 <- paste0(export_bin_message, collapse = "")

## Decode Elias Gamma encoded integers
decode2 <- decodeGamma(decode1, args$codestop)

## Create empty list and apply inverse run length encoding to fill elements
decode3 <- list()
for (i in 1:length(decode2)) {
  if (i %% 2 == 0 ) {
    decode3[i] <- paste(rep("1", decode2[i]), collapse = "")
  } else {
    decode3[i] <- paste(rep("0", decode2[i]), collapse = "")
  }
}

## Split into segments of bits, of length = width of image
decode4 <- paste0(decode3, collapse = "")
decode5 <- strsplit(decode4, '(?<=.{60})', perl=TRUE)[[1]]

## Create a binary matrix (character)
decode6 <- matrix(decode5, nrow = length(decode5), ncol = 1)
decode7 <- sapply(decode6, strsplit, split = "")

## Create a binary matrix (numeric)
decode8 <- lapply(decode7, as.numeric)
decode9 <- NULL
for (i in 1:length(decode8)) {
  decode9 <- rbind(decode9, decode8[[i]])
}

### Visualise the binary matrix ###

pdf(file = portrait_output, width = 6, height = 9, pointsize = 16, colormodel = "cmyk")
image(rotate_mat(rotate_mat(rotate_mat(decode9))), col = grey(c(0,1)))
invisible(dev.off())

cat("** END **\n")
quit(save = "no", status = 0)