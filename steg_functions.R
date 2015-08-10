## Code is dependent on the following libraries, please install via install.packages if missing
## eg: install.packages("stringr")
library(pixmap)
library(stringr)
library(plotrix)
library(argparse)

## STEG main process function
process <- function(x, y) {
      split <- strsplit(x$V1, split = "")                     ## Convert reads into columns of bases
      cond <- lapply(split, function (x) length(x) == y)      ## Set to filter for reads of length y
      split2 <- split[unlist(cond)]                           ## Apply filter condition to the list
      df <- data.frame(matrix(unlist(split2), nrow=length(split2), byrow=T))   
      dg <- matrix(NA, ncol=length(df[1,]), nrow=5)           ## Create empty matrix
      dimnames(dg)[[1]] <- c("A","T","C","G","-")             ## Name columns of empty matrix
      for (i in 1:length(df[1,])){                            ## Loop to set factor levels 
            for(j in  c("A","T","C","G","-")){
                  dg[j, i] <- sum(df[,i]==j)   
            }
      }
      dg <- dg[-5, ]                            ## Remove the row containing counts of deletions
      dg_perc <- dg                             ## Create dummy data.frame for base % at positions
      sapply(1:length(dg[1,]), FUN=function(i){ dg_perc[,i] <<- 100*dg[,i]/sum(dg[,i])})     
}

## Baking soda
## Function to convert from a binary code to primary DNA sequence
codetobase = function(code){
      code = str_replace_all(code, "000", "AT");  
      code = str_replace_all(code, "111", "TA");
      code = str_replace_all(code, "011", "CG");
      code = str_replace_all(code, "100", "GC");
      code = str_replace_all(code, "001", "5G");
      code = str_replace_all(code, "110", "G5");
      code = str_replace_all(code, "010", "6G");
      code = str_replace_all(code, "101", "G6");
      return (code);
}

## Baking soda
## Insertion function (adapted from ferdinanrd.kraft on stackoverflow)
insert.at <- function(a, pos, dots){
      stopifnot(length(dots)==length(pos))
      result <- vector("list",2*length(pos)+1)
      result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
      result[c(FALSE,TRUE)] <- dots
      unlist(result)
}

## Baking soda
## Function from strsplit help page, for reversing strings
strReverse <- function(x)
      sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")

## Baking soda
## Function to convert from a ascii numerals to binary
numtobin = function(num){
      num = str_replace_all(num, "0", "0000");  
      num = str_replace_all(num, "1", "0001"); 
      num = str_replace_all(num, "2", "0010"); 
      num = str_replace_all(num, "3", "0011"); 
      num = str_replace_all(num, "4", "0100"); 
      num = str_replace_all(num, "5", "0101"); 
      num = str_replace_all(num, "6", "0110"); 
      num = str_replace_all(num, "7", "0111"); 
      num = str_replace_all(num, "8", "1000"); 
      num = str_replace_all(num, "9", "1001"); 
      return (num);
}

## Baking soda
## Function to insert '0' character between all values of a vector
insert_black <- function(x) {
      if (length(x) > 1) {
            append(insert.at(x, seq(from = 1, to = (length(x)-1), by = 1), 
                             strsplit(rep("0", (length(x)-1)), split = "")), "0")
      }
      else {
            append(x, "0")
      }
}

## Baking soda
## Function to insert '1' character between all values of a vector
insert_white <- function(x) {
      if (length(x) > 1) {
            append(insert.at(x, seq(from = 1, to = (length(x)-1), by = 1), 
                             strsplit(rep("1", (length(x)-1)), split = "")), "1")
      }
      else {
            append(x, "1")
      }
}

## Baking soda
## Define a function to convert decimal to binary
## Arguments are x = number, num_bits = number of bits to display
decToBin <- function(x, num_bits) {
      print_from <- (32 - num_bits + 1)
      temp1 <- strsplit(paste(rev(intToBits(x))), split = "")
      paste0(sapply(temp1[print_from:32], `[[`, 2), collapse = "")
}

## Baking soda
## Function to split large numbers into components
bigNumSplit <- function(x) {
      mod_x <- NULL 
      i <- 1
      while (i <= nchar(x)) {
            if (substr(x, start = i, stop = i+1) < 64) {
                  mod_x <- append(mod_x, substr(x, start = i, stop = i+1))
                  i <- i+2
            }
            else {
                  mod_x <- append(mod_x, substr(x, start = i, stop = i))
                  i <- i+1
            }
      }
      return (mod_x)
}

## Baking soda
## Function that outputs a vector of values read from a matrix in a zigzag fashion
zigzag <- function(input) {
      
      h <- 1   ## horizontal index
      v <- 1   ## vertical index
      
      vmin <- 1
      hmin <- 1
      
      vmax <- nrow(input)
      hmax <- ncol(input)
      
      i <- 1
      
      output <- NULL
      
      while ((v <= vmax) & (h <= hmax)) {
            if (((h+v) %% 2) == 0) {
                  if (v == vmin) {
                        output[i] <- input[v,h]   
                        if (h == hmax) {
                              v <- v + 1
                        }
                        else {
                              h <- h + 1
                        }
                        i <- i + 1
                  }
                  else if ((h == hmax) & (v < vmax)) {
                        output[i] <- input[v,h]
                        v <- v + 1
                        i <- i + 1
                  }
                  else if ((v > vmin) & (h < hmax)) {
                        output[i] <- input[v,h]
                        v <- v - 1
                        h <- h + 1
                        i <- i + 1
                  }
            }
            else {
                  if ((v == vmax) & (h <= hmax)) {   
                        output[i] <- input[v, h]
                        h <- h + 1
                        i <- i + 1
                  }
                  else if (h == hmin) {                  
                        output[i] <- input[v, h]
                        
                        if (v == vmax) {
                              h <- h + 1
                        }
                        else {
                              v <- v + 1
                        }
                        
                        i = i + 1
                  }
                  else if ((v < vmax) & (h > hmin)) {
                        output[i] <- input[v,h]
                        v <- v + 1
                        h <- h - 1
                        i <- i + 1
                  }
            }
            if ((v == vmax) & (h == hmax)) {
                  output[i] <- input[v,h]
                  break
            }
      }
      return(output)
}

## Baking soda / decode_images
## Convert a decimal number to uanry
unary <- function(x) {
      N <- as.integer(log2(x))
      return(append(rep("0", N), "1"))
}

## Baking soda / decode_images
## Convert a decimal number to binary
## Optionally, can set the number of bits to display
decToBin <- function(x, num_bits = NULL) {
      temp1 <- strsplit(paste(rev(intToBits(x))), split = "")
      if (is.null(num_bits)) {
            temp2 <- paste0(sapply(temp1, `[[`, 2), collapse = "")
            return(substr(temp2, regexpr("[^0]", temp2), nchar(temp2)) )
      }
      else {
            print_from <- (32 - num_bits + 1)
            return(paste0(sapply(temp1[print_from:32], `[[`, 2), collapse = ""))
      }
}

## Baking soda / steg_decode_images
## Generic Elias encoding function
eliasGamma <- function(x) {
      if (x == 0) {
            return("0")
      }
      else {
            bin <- decToBin(x)
            encoded <- append(paste0(unary(x), collapse = ""), 
                              substr(bin, start = 2, stop = nchar(bin)))
            return(paste0(encoded, collapse = ""))
      }
}

## Baking soda
## Elias delta encoding
eliasDelta <- function(x) {
      N <- as.integer(log2(x))   ## Let N=log2 x be the highest power of 2, so 2N <= x < 2N+1
      Bin <- decToBin(x)         ## Convert x to binary
      gamma <- eliasGamma(N + 1)  ## Gamma encode (N+1)
      remNbin <- substr(Bin, start = 2, stop = nchar(Bin))  ## Obtain remaining N binary digits of x 
      paste0(append(gamma, remNbin), collapse = "")
}

## Baking soda
## Elias Omega encoding (recursive encoding)
eliasOmega <- function(x) {
      code <- "0"
      if (x == 1) {
            return(code)
      } else {
            
            while (nchar(decToBin(x)) >= 2) {
                  code <- append(code,decToBin(x), after = 0)
                  x <- nchar(code[1]) - 1
            }
            paste0(code, collapse = "")
      }
}

## Lemonjuice
## Define function: convert strings of asciiary into into strings of ascii chaarcters
asciitobin = function(ascii){
      ascii = str_replace_all(ascii, "A", "01000001");  ## Each line gives asciiary for ASCII char
      ascii = str_replace_all(ascii, "B", "01000010");
      ascii = str_replace_all(ascii, "C", "01000011");
      ascii = str_replace_all(ascii, "D", "01000100");
      ascii = str_replace_all(ascii, "E", "01000101");
      ascii = str_replace_all(ascii, "F", "01000110");
      ascii = str_replace_all(ascii, "G", "01000111");
      ascii = str_replace_all(ascii, "H", "01001000");
      ascii = str_replace_all(ascii, "I", "01001001");
      ascii = str_replace_all(ascii, "J", "01001010");
      ascii = str_replace_all(ascii, "K", "01001011");
      ascii = str_replace_all(ascii, "L", "01001100");
      ascii = str_replace_all(ascii, "M", "01001101");
      ascii = str_replace_all(ascii, "N", "01001110");
      ascii = str_replace_all(ascii, "O", "01001111");
      ascii = str_replace_all(ascii, "P", "01010000");
      ascii = str_replace_all(ascii, "Q", "01010001");
      ascii = str_replace_all(ascii, "R", "01010010");
      ascii = str_replace_all(ascii, "S", "01010011");
      ascii = str_replace_all(ascii, "T", "01010100");
      ascii = str_replace_all(ascii, "U", "01010101");
      ascii = str_replace_all(ascii, "V", "01010110");
      ascii = str_replace_all(ascii, "W", "01010111");
      ascii = str_replace_all(ascii, "X", "01011000");
      ascii = str_replace_all(ascii, "Y", "01011001");
      ascii = str_replace_all(ascii, "Z", "01011010");
      ascii = str_replace_all(ascii, "a", "01100001"); 
      ascii = str_replace_all(ascii, "b", "01100010"); 
      ascii = str_replace_all(ascii, "c", "01100011"); 
      ascii = str_replace_all(ascii, "d", "01100100"); 
      ascii = str_replace_all(ascii, "e", "01100101"); 
      ascii = str_replace_all(ascii, "f", "01100110"); 
      ascii = str_replace_all(ascii, "g", "01100111"); 
      ascii = str_replace_all(ascii, "h", "01101000"); 
      ascii = str_replace_all(ascii, "i", "01101001"); 
      ascii = str_replace_all(ascii, "j", "01101010"); 
      ascii = str_replace_all(ascii, "k", "01101011"); 
      ascii = str_replace_all(ascii, "l", "01101100"); 
      ascii = str_replace_all(ascii, "m", "01101101"); 
      ascii = str_replace_all(ascii, "n", "01101110"); 
      ascii = str_replace_all(ascii, "o", "01101111"); 
      ascii = str_replace_all(ascii, "p", "01110000"); 
      ascii = str_replace_all(ascii, "q", "01110001"); 
      ascii = str_replace_all(ascii, "r", "01110010"); 
      ascii = str_replace_all(ascii, "s", "01110011"); 
      ascii = str_replace_all(ascii, "t", "01110100"); 
      ascii = str_replace_all(ascii, "u", "01110101"); 
      ascii = str_replace_all(ascii, "v", "01110110"); 
      ascii = str_replace_all(ascii, "w", "01110111"); 
      ascii = str_replace_all(ascii, "x", "01111000"); 
      ascii = str_replace_all(ascii, "y", "01111001"); 
      ascii = str_replace_all(ascii, "z", "01111010"); 
      ascii = str_replace_all(ascii, "\n", "00001001");   ## newline
      ascii = str_replace_all(ascii, " ", "00100000");    ## space
      ascii = str_replace_all(ascii, ",", "00101100");
      ascii = str_replace_all(ascii, fixed("."), "00101110");  ## Prevents str_replace_all using regexp
      ascii = str_replace_all(ascii, "!", "00100001");
      ascii = str_replace_all(ascii, "-", "00101101");
      ascii = str_replace_all(ascii, "'", "00100111");
      ascii = str_replace_all(ascii, ";", "00111011");
      return (ascii);
}

## Lemon juice
## Define function: creates matrix of binary messages (note: not actually an R matrix)
wrap_convert <- function(files) {
      
      for (file_name in files) {
            data <- readChar(file_name, nchars = file.info(file_name)$size)
            
            split_data_list <- strsplit(data, split = "")
            
            split_data <- unlist(split_data_list)
            
            ## Create empty character vector to accept converted message
            binary_conv <- character(length = 0)   
            
            for (i in 1:length(split_data)) {
                  binary_conv[i] <- (asciitobin(split_data[[i]]))
            }
            
            ## Prep variable for printing to screen and output to file
            binary_messages_vec <<- c(binary_messages_vec, paste0(binary_conv, collapse = ""))
            
      }
}

## Function to rotate a matrix (required as image() seems to rotate when plotting)
rotate_mat <- function(x) t(apply(x, 2, rev))

## steg_decode_text
## Define function that converts strings of binary into into strings of ascii chaarcters
bintoascii = function(bin){
      bin = str_replace_all(bin, "01000001", "A");  ## Each line gives binary for an ASCII character
      bin = str_replace_all(bin, "01000010", "B");
      bin = str_replace_all(bin, "01000011", "C");
      bin = str_replace_all(bin, "01000100", "D");
      bin = str_replace_all(bin, "01000101", "E");
      bin = str_replace_all(bin, "01000110", "F");
      bin = str_replace_all(bin, "01000111", "G");
      bin = str_replace_all(bin, "01001000", "H");
      bin = str_replace_all(bin, "01001001", "I");
      bin = str_replace_all(bin, "01001010", "J");
      bin = str_replace_all(bin, "01001011", "K");
      bin = str_replace_all(bin, "01001100", "L");
      bin = str_replace_all(bin, "01001101", "M");
      bin = str_replace_all(bin, "01001110", "N");
      bin = str_replace_all(bin, "01001111", "O");
      bin = str_replace_all(bin, "01010000", "P");
      bin = str_replace_all(bin, "01010001", "Q");
      bin = str_replace_all(bin, "01010010", "R");
      bin = str_replace_all(bin, "01010011", "S");
      bin = str_replace_all(bin, "01010100", "T");
      bin = str_replace_all(bin, "01010101", "U");
      bin = str_replace_all(bin, "01010110", "V");
      bin = str_replace_all(bin, "01010111", "W");
      bin = str_replace_all(bin, "01011000", "X");
      bin = str_replace_all(bin, "01011001", "Y");
      bin = str_replace_all(bin, "01011010", "Z");
      bin = str_replace_all(bin, "01100001", "a"); 
      bin = str_replace_all(bin, "01100010", "b"); 
      bin = str_replace_all(bin, "01100011", "c"); 
      bin = str_replace_all(bin, "01100100", "d"); 
      bin = str_replace_all(bin, "01100101", "e"); 
      bin = str_replace_all(bin, "01100110", "f"); 
      bin = str_replace_all(bin, "01100111", "g"); 
      bin = str_replace_all(bin, "01101000", "h"); 
      bin = str_replace_all(bin, "01101001", "i"); 
      bin = str_replace_all(bin, "01101010", "j"); 
      bin = str_replace_all(bin, "01101011", "k"); 
      bin = str_replace_all(bin, "01101100", "l"); 
      bin = str_replace_all(bin, "01101101", "m"); 
      bin = str_replace_all(bin, "01101110", "n"); 
      bin = str_replace_all(bin, "01101111", "o"); 
      bin = str_replace_all(bin, "01110000", "p"); 
      bin = str_replace_all(bin, "01110001", "q"); 
      bin = str_replace_all(bin, "01110010", "r"); 
      bin = str_replace_all(bin, "01110011", "s"); 
      bin = str_replace_all(bin, "01110100", "t"); 
      bin = str_replace_all(bin, "01110101", "u"); 
      bin = str_replace_all(bin, "01110110", "v"); 
      bin = str_replace_all(bin, "01110111", "w"); 
      bin = str_replace_all(bin, "01111000", "x"); 
      bin = str_replace_all(bin, "01111001", "y"); 
      bin = str_replace_all(bin, "01111010", "z"); 
      bin = str_replace_all(bin, "00001001", "\n");   ## newline
      bin = str_replace_all(bin, "00100000", " ");    ## space
      bin = str_replace_all(bin, "00101100", ",");
      bin = str_replace_all(bin, "00101110", ".");
      bin = str_replace_all(bin, "00100001", "!");
      bin = str_replace_all(bin, "00101101", "-");
      bin = str_replace_all(bin, "00100010", "\"");
      bin = str_replace_all(bin, "00100111", "'");
      bin = str_replace_all(bin, "00001010", "\n");   ## linefeed (newline)
      bin = str_replace_all(bin, "00111011", ";");
      return (bin);
}

## Lemonjuice (inverse of baking soda)
## Function to convert from a primary DNA sequence to a binary matrix
basetocode = function(base){
      base = str_replace_all(base, "AT", "000");  
      base = str_replace_all(base, "TA", "111");
      base = str_replace_all(base, "CG", "011");
      base = str_replace_all(base, "GC", "100");
      base = str_replace_all(base, "5G", "001");
      base = str_replace_all(base, "G5", "110");
      base = str_replace_all(base, "6G", "010");
      base = str_replace_all(base, "G6", "101");
      return (base);
}

## steg_decode_images
## Function to decode a series of Elias Gamma numbers
## Input a single character string of numbers
## Termination code of "111"
decodeGamma <- function(code, finish = "111") {
      split_code <- unlist(strsplit(code, split = ""))
      decoded <- numeric()
      element_index <- 1
      
      while (element_index < nchar(code)) {
            N <- 0                                     ## Set the "zero counter" to zero
            
            ## For an element in code, if 0 then increment the "zero counter"
            for (element in split_code[element_index:length(split_code)]) {       
                  if (finish == "111" & N == 0 & element == "1" & split_code[element_index + 1] == "1" & split_code[element_index + 2] == "1") {
                        print("End of code found (code element starts '111')")
                        return(decoded)
                  } else if (finish == "1111" & N == 0 & element == "1" & split_code[element_index + 1] == "1" & split_code[element_index + 2] == "1" & split_code[element_index + 3] == "1") {
                        print("End of code found (code element starts '1111')")
                        return(decoded)
                  } else if (finish == "11111" & N == 0 & element == "1" & split_code[element_index + 1] == "1" & split_code[element_index + 2] == "1" & split_code[element_index + 3] == "1" & split_code[element_index + 4] == "1") {
                        print("End of code found (code element starts '11111')")
                        return(decoded)
                  } else if (element == "0") {         ## If one then break out of this loop
                        N = N + 1 
                        element_index <- element_index + 1
                  } else if (element == "1") {
                        break
                  }
            }
            
            ## Read the next N elements of split_code
            bin_value <- paste0(split_code[element_index:(element_index + N)], collapse = "")   
            ## Convert the pasted elements into decimal
            dec_value <- strtoi(bin_value, base = 2L)   
            ## Append this value to the decoded list
            decoded <- append(decoded, dec_value)     
            ## Increment the element counter by elements read
            element_index <- element_index + (N+1)                                                 
      }
      
      print(element_index)
      return(decoded)                         ## Return the vector of decoded decimal values
}
