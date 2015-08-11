## Source file containing functions and library dependencies
suppressWarnings(source("STEG_functions.R"))

output_prefix <- format(Sys.time(), "%Y%m%d_%H%M")

## Read in command line arguments
args <- commandArgs(TRUE)

cat("\n***Options chosen***\n")
cat("Scanning direction:", args[1], "\n", sep = " ")
cat("Encoding method:", args[2], "\n\n", sep = " ")

## Choose scan method for rle: "horizontal", "vertical", "zigzag"
scan_method <- args[1]

## Choose encoding method for rle: "7bit", "gamma", "delta"
encoding <- args[2]

## File names of input images, which must be in .pbm format
file1 <- args[3]
file2 <- args[4]
file3 <- args[5]
splitlength <- as.numeric(args[6])

## Name of output text file containing the DNA sequences
output_name <- paste0(output_prefix, "_bakingsoda.out", collapse = "")

## Read in the binary image files
suppressWarnings(image1 <- read.pnm(file1))
suppressWarnings(image2 <- read.pnm(file2))
suppressWarnings(image3 <- read.pnm(file3))

## Obtain a list of binary matrices describing the files
img_mat_list <- list(img1 = getChannels(image1, colors = "grey"), 
                     img2 = getChannels(image2, colors = "grey"), 
                     img3 = getChannels(image3, colors = "grey"))

if (scan_method == "zigzag") {
      img_zigzagged_list <- lapply(img_mat_list, zigzag)
}

## Create empty lists ready to receive information from the binary oligos
oligos <- list(NULL, NULL, NULL)
bin_codes <- list(NULL, NULL, NULL)

for (i in 1:3) {    
      ## if zigzag then rle directly as already a vector
      ## if horizontal or vertical rle
      if (scan_method == "zigzag") {
            rle_num <- rle(img_zigzagged_list[[i]])[[1]]
      } else if (scan_method == "vertical") {
            rle_num <- rle(as.vector(img_mat_list[[i]]))[[1]]
      } else if (scan_method == "horizontal") {
            rle_num <- rle(as.vector(rotate_mat(img_mat_list[[i]])))[[1]]
      } else { stop("Please select a scanning method")
      }
      
      if (encoding == "7bit") {
            rle_num_bin <- NULL
            
            ## Replace decimal with binary
            ## Direct replacement if less than 64 (ie: can encode in six bits)
            ## Split and replace if greater than 63 (ie: cannot encode in six bits)
            for (z in 1:length(rle_num)) {
                  if (rle_num[z] < 64) {
                        rle_num_bin <- append(rle_num_bin, decToBin(rle_num[z], 6))
                        if (z %% 2 == 0) { rle_num_bin <- append(rle_num_bin, "0") }    ## 0=black
                        else if (z %% 2 != 0) {rle_num_bin <- append(rle_num_bin, "1")} ## 1=white
                  }
                  else if (rle_num[z] > 63) {
                        tmp_split <- bigNumSplit(rle_num[z])      ## Split into codeable numbers
                        tmp_bin <- sapply(tmp_split, decToBin, 6)             ## Convert to binary
                        if (z %% 2 == 0) {tmp_bin2 <- insert_black(tmp_bin)}                   
                        else if (z %% 2 != 0) {tmp_bin2 <- insert_white(tmp_bin)}                        
                        rle_num_bin <- append(rle_num_bin, tmp_bin2)   
                  }
            }
            
            bin_code <- paste0(rle_num_bin, collapse = "")
            bin_codes[[i]] <- bin_code     ## Captures output of intermediates in the for loop
            
            ## Splits the bin_code into oligo sized chunks of binary
            oligos[[i]] <- strsplit(bin_code, paste0("(?<=.{", splitlength, "})", collapse = ""), perl=TRUE)[[1]]

      } else if (encoding == "gamma") {
            gamma_coded <- NULL
            
            for (value in rle_num) {
                  gamma_coded <- append(gamma_coded, (eliasGamma(value)))    
            }
            
            pastedGamma <- paste0(gamma_coded, collapse = "")
            bin_codes[[i]] <- pastedGamma  ## Captures output of intermediates in the for loop                   
            oligos[[i]] <- strsplit(pastedGamma, paste0("(?<=.{", splitlength, "})", collapse = ""), perl=TRUE)[[1]]      
            
      } else if (encoding == "delta") {
            delta_coded <- NULL
            
            for (value in rle_num) {
                  delta_coded <- append(delta_coded, (eliasDelta(value)))    
            }
            
            pastedDelta <- paste0(delta_coded, collapse = "")
            bin_codes[[i]] <- pastedDelta  ## Captures output of intermediates in the for loop
            
            oligos[[i]] <- strsplit(pastedDelta, paste0("(?<=.{", splitlength, "})", collapse = ""), perl=TRUE)[[1]]    
      }
      
}

## Create an empty list ready to receive information from the binary oligos
max_length <- max(length(oligos[[1]]), length(oligos[[2]]), length(oligos[[3]]))
min_length <- min(length(oligos[[1]]), length(oligos[[2]]), length(oligos[[3]]))
mat_list_names <- c(1:min_length) 
mat_list <- sapply(mat_list_names, function(x) NULL)

## Append an "end of message" sequence to each picture
## If number of oligos is not equal, segregates excess oligos beyond minimum number into "run_over"

run_over <- list(NULL,NULL,NULL)
for (i in 1:3) {
      oligos[[i]][length(oligos[[i]])] <- paste0(append(oligos[[i]][length(oligos[[i]])], "111"), 
                                                 collapse = "")
      
      if (length(oligos[[i]]) > min_length) {
            run_over[[i]] <- oligos[[i]][(min_length+1):max_length]
      }
      oligos[[i]] <- oligos[[i]][1:min_length]
}

## Fill in final row of each oligo (with random bases) to generate oligos of equal length
for (i in 1:3) {
      tmp <- oligos[[i]][length(oligos[[i]])]
      if (nchar(oligos[[i]][length(oligos[[i]])]) != splitlength) {
            oligos[[i]][length(oligos[[i]])] <- capture.output(cat(tmp, 
                                                sample(c(0,1), (splitlength-nchar(tmp)), replace = TRUE),
                                                sep = ""))
      }
}

## Create a list where each element contains a row of all three images
for (i in 1:min_length) {
      mat_list[[i]] <- rbind(oligos[[1]][i], oligos[[2]][i], oligos[[3]][i])
}

mat_list_split <- lapply(mat_list, strsplit, split = "")

for (i in 1:length(mat_list_split)) {
      mat_list_split[[i]] <- data.frame(matrix(unlist(mat_list_split[[i]]), 
                                               nrow=length(mat_list_split[[i]]), 
                                               byrow=T), stringsAsFactors = FALSE)
}

## Paste together positions for all three images to form a binary code
for (i in 1:length(mat_list_split)) {
      mat_list_split[[i]] <- apply(mat_list_split[[i]], 2, paste0, collapse = "")  
}

## Convert from binary code to primary sequence
seq_code <- lapply(mat_list_split, codetobase)

## Create nested list; split into top and bottom strand
split_seq_code <- lapply(seq_code, strsplit, split = "")

## Unlist with recursive false to generate a single list
split_seq_list <- unlist(split_seq_code, recursive = FALSE)

## Convert to data frame, where each row is the sequence for oligos encoding itself
base_df <- data.frame(matrix(unlist(split_seq_list), nrow = length(split_seq_list), byrow = TRUE), 
                      stringsAsFactors = FALSE)
colnames(base_df) <- c("top_strand", "bottom_strand")

## Create character vector of all top sequences pasted together, then all bottom in the same manner
top <- paste0(base_df$top_strand, collapse = "")
bottom <- paste0(base_df$bottom_strand, collapse = "")   

## Chop the sequence character vectors into portions of oligo length (120))
top_sequences <- strsplit(top, paste0("(?<=.{", splitlength, "})", collapse = ""), perl=TRUE)[[1]]
bottom_sequences <- strsplit(bottom, paste0("(?<=.{", splitlength, "})", collapse = ""), perl=TRUE)[[1]]

## Use the strReverse function to obtain bottom_sequences in the 5' to 3' direction
five2three_bottom_sequences <- strReverse(bottom_sequences)

## Export to .txt file
fileConn <- file(output_name)
writeLines(c(">all top sequences|5' to 3'", top_sequences, ">all bottom sequences|5' to 3'", 
            five2three_bottom_sequences), fileConn)
close(fileConn)

## Image matrix (0 = black 1 = white)
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE))
image(rotate_mat(img_mat_list[[1]]), col = grey(c(0,1)))
image(rotate_mat(img_mat_list[[2]]), col = grey(c(0,1)))
image(rotate_mat(img_mat_list[[3]]), col = grey(c(0,1)))

cat("***Compression achieved***\n")
cat(args[3], "compressed from", dim(img_mat_list[[1]])[1]*dim(img_mat_list[[1]])[2], "to", nchar(bin_codes[[1]]), "\n", sep = " ")
cat(args[4], " compressed from", dim(img_mat_list[[2]])[1]*dim(img_mat_list[[2]])[2], "to", nchar(bin_codes[[2]]), "\n", sep = " ")
cat(args[5], " compressed from", dim(img_mat_list[[3]])[1]*dim(img_mat_list[[3]])[2], "to", nchar(bin_codes[[3]]), "\n\n", sep = " ")
cat("***Output created***\n")
cat(output_name, "\n\n")

if (!is.null(run_over[[1]]) | !is.null(run_over[[2]]) | !is.null(run_over[[3]])) {
      cat("***Run over binary strings***\n")
      print(run_over)
}
