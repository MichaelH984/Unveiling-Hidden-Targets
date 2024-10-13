#!/bin/bash

# Create bin and output directories if they do not exist
mkdir -p ../bin
mkdir -p ../output

# Path to the input files in the Raw data directory
GWIPS_VIZ_BED="../Raw data/gwipsvizRiboseq.bed"
CCDS_BED="../Raw data/CCDS.bed"

# Step 1: Intersect and filter unique entries
bedtools intersect -a "$GWIPS_VIZ_BED" -b "$CCDS_BED" -wa | sort | uniq > "../bin/unique_filtered_output.bed"

# Step 2: Move the third column from the first row to the rest
awk 'NR==1 {prev=$3; $3=""} NR>1 {temp=$3; $3=prev; prev=temp} 1' "../bin/unique_filtered_output.bed" > "../bin/unique_moved.bed"

# Step 3: Remove specific columns using awk
awk '{$1=$4=$5=""; $0=$0}1' OFS="\t" "../bin/unique_moved.bed" > "../bin/unique_deleted.bed"

# Step 4: Remove the first and last line using sed
sed '1d;$d' "../bin/unique_deleted.bed" > "../bin/unique_sorted.bed"

# Step 5: Calculate the third column based on the first two
awk '{ $3 = $1 - $2 }1' "../bin/unique_sorted.bed" > "../bin/unique_third.bed"

# Step 6: Create frequency table from the third column
awk '{counts[$3]++} END {for (score in counts) print score, counts[score]}' "../bin/unique_third.bed" | sort -n > "../bin/frequency_table.txt"

# Step 7: Filter frequency table for counts greater than or equal to 0
awk '$2 >= 0' "../bin/frequency_table.txt" > "../output/FQtable.txt"

# Step 8: Calculate the total count of scores
total_count=$(awk '{sum += $2} END {print sum}' "../output/FQtable.txt")

# Step 9: Calculate cumulative frequency and find scores for each 1% increment from 90% to 100%
awk -v total_count="$total_count" '
BEGIN { 
    print "Percentage\tScore";
    print "-----------\t-----";
}
{
    cumulative += $2; 
    for (i = 90; i <= 100; i++) {
        if (cumulative / total_count >= i / 100 && !scores[i]) {
            scores[i] = $1;  # Capture the score at each 1% increment
            print i "%\t" scores[i]; 
            # Removed the blank line to print each entry continuously
        }
    }
}' "../output/FQtable.txt" > "../output/percentile_scores.txt"