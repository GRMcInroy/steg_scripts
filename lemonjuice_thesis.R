## R script containing functions and library dependencies
suppressWarnings(source("STEG_functions.R"))

## Start the clock!
ptm <- proc.time()

input_file1 <- "asciitest1.txt"
input_file2 <- "asciitest2.txt"
input_file3 <- "asciitest3.txt"
output_file <- "DNA_sequence.txt"

binary_messages_vec <- character(length = 0)                                                              

## Create a matrix (not R matrix) of binarymessages
wrap_convert(c(input_file1, input_file2, input_file3))                                                   
binary_messages <- data.frame(messages = binary_messages_vec, stringsAsFactors = FALSE)                          

## Create character vector containing the encoded messages one position per element
split_bin <- strsplit(binary_messages$messages, split = "") ## list of 3 char vectors, 1 bit/element
split_binary_messages <- data.frame(matrix(unlist(split_bin), nrow=length(split_bin), byrow=T), stringsAsFactors = FALSE)    
code_vec <- apply(split_binary_messages, 2, paste0, collapse = "")                                        

## Convert the encoded messages into DNA sequence
base_code <- codetobase(code_vec)

## Split into two sequences for a top and a bottom strand
base_list <- NULL

for (i in 1:length(base_code)) {
  base_list[i] <- strsplit(base_code[i], split = "")
}

base_df <- data.frame(matrix(unlist(base_list), nrow = length(base_list), byrow = T), stringsAsFactors = FALSE)
colnames(base_df) <- c("top_strand", "bottom_strand")

## Reformat for exporting to .txt file
top_sequence <- paste0(base_df$top_strand, collapse = "")
bottom_sequence <- paste0(rev(base_df$bottom_strand), collapse = "") ## Rev to give 5' to 3' seq

## Export to .txt file
fileConn <- file(output_file)
writeLines(c(">top_sequence|5' to 3'", top_sequence, ">bottom_sequence|5' to 3'", bottom_sequence), fileConn)
close(fileConn)

## Stop the clock!
proc.time() - ptm