#!/bin/bash
set -e

# Have a directory containing the reference files
# Have a directory containing the barcode files
# Use "cat input.file" and pipe (|) into this script
# Select the adapter to trim as first argument

# cat data.fastq | autoscript_split_DFT_new.sh 6
# (to trim for adapter 6 during fastx_clipper)

# Define the date variable
DATE="$(date +"%Y%m%d_%k%M%S")"
# Define the current working directory
myDIR="$(pwd -P)"

# Create a temp_split folder
mkdir -p "${DATE}"_temp_split
# Create an output folder
mkdir -p "${DATE}"_autoscript_output

# Define the trimming sequence
trimSEQ="$1"
# Define the forward barcodes file for fastx_barcode_splitter.pl
fwd_bcode="DFT_barcodes.txt"
# Define the reverse barcodes file for fastx_barcode_splitter.pl
rev_bcode="DFT_barcodes_revcomp.txt"
# Default adapter sequence
adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTC"


# Trim adapter based upon index given in argument 1
if [ "${trimSEQ}" = "1" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"
elif [ "${trimSEQ}" = "2" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "3" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "4" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "5" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "6" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "7" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "8" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "9" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "10" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "11" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "12" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "13" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "14" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "15" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "16" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "18" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "19" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "20" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "21" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "22" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "23" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "25" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
elif [ "${trimSEQ}" = "27" ]; then
	adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG"
	echo "$adapter"	
fi

# Trimming lines
fastx_clipper -Q 33 -v -a "${adapter}" -l 70 | \
fastq_quality_filter -Q 33 -v -q 30 -p 75 -o temp_trimmed.fastq

# Split file by inline barcodes (credit to Hannon lab for fastxtools)
# First by barcode at beginning of line, then rev comp at end of line
cat temp_trimmed.fastq | fastx_barcode_splitter.pl --bcfile "${myDIR}"/barcodes/"${fwd_bcode}" \
--bol --exact --prefix "${myDIR}"/"${DATE}"_temp_split/fwd_ --suffix ".fastq"

cat temp_trimmed.fastq | fastx_barcode_splitter.pl --bcfile "${myDIR}"/barcodes/"${rev_bcode}" \
--eol --exact --prefix "${myDIR}"/"${DATE}"_temp_split/rev_ --suffix ".fastq"

# Loop to merge the pairs of split files
for j in {1..28}; do
	i=$(printf "%02d" "$j")
	cat "${myDIR}"/"${DATE}"_temp_split/fwd_"${i}".fastq "${myDIR}"/"${DATE}"_temp_split/rev_"${i}".fastq \
	> "${myDIR}"/"${DATE}"_autoscript_output/"${i}".fastq
done
	
# Delete the temporary trimmed fastq output
rm -f temp_trimmed.fastq

# Delete the temp_split folder
rm -dfR "${DATE}"_temp_split

# Run the aligner

for j in {1..28}; do
	i=$(printf "%02d" "$j")
	echo FASTQ FILE "${myDIR}"/"${DATE}"_autoscript_output/"${i}".fastq
	echo REFERENCE "${myDIR}"/reference/ref_"${i}".txt
	java -jar /Users/Gordon/prog/ShortRefGlobalAligner/ShortRefGlobalAligner.jar \
	-f "${myDIR}"/"${DATE}"_autoscript_output/"${i}".fastq -r "${myDIR}"/reference/ref_"${i}".txt \
	-n "${i}" \
	> "${myDIR}"/"${DATE}"_autoscript_output/SRGA_"${i}".txt
done