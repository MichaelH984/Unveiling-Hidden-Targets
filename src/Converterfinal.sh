#!/bin/bash

# Define input and output directories for processing .bw files and storing results.
input_dir="../Raw data"  # Input directory containing .bw files
output_dir="../bin"      # Output directory for generated files

# Checks if the output and bin directory exists; create it if it doesn't.
mkdir -p "$output_dir"
mkdir -p "$bin_dir"

# Make the bigWigToWig script executable to convert .bw files to .wig format.
chmod +x bigWigToWig

# Loop through each .bw file in the input directory for processing.
for bw_file in "$input_dir"/*.bw  
do
    # Extract the base filename (without the .bw extension) for naming output files.
    base_name=$(basename "$bw_file" .bw)

    # Define output file paths for the .wig and .bed files.
    wig_file="$output_dir/$base_name.wig"  # Output .wig file path
    bed_file_prefix="$output_dir/$base_name"  # Base prefix for subsequent .bed files

    # Convert the .bw file to .wig format using the bigWigToWig tool.
    ../src/bigWigToWig "$bw_file" "$wig_file"
    echo "Converted $bw_file to $wig_file"  # Notify that conversion to .wig is complete

    # Convert the generated .wig file to .bed format.
    wig2bed --zero-indexed <"$wig_file" > "$bed_file_prefix.bed"
    echo "Converted $wig_file to $bed_file_prefix.bed"  # Notify that conversion to .bed is complete

    # Perform a series of transformations on the .bed file using various tools.
    awk '$5 >= 10' "$bed_file_prefix.bed" > "$bed_file_prefix+v1.bed"  # Filter entries based on a score threshold
    
    bedtools merge -i "$bed_file_prefix+v1.bed" -d 44 > "$bed_file_prefix+v2.bed"  # Merge overlapping intervals
    
    bedtools subtract -A -a "$bed_file_prefix+v2.bed" -b "$input_dir/CCDS.bed" > "$bed_file_prefix+v3.bed"  # Subtract CCDS regions
    
    awk '$3 - $2 != 1' "$bed_file_prefix+v3.bed" > "$bed_file_prefix+v4.bed"  # Remove intervals of length 1
    
    bedtools intersect -wa -a "$bed_file_prefix+v4.bed" -b "$input_dir/sno_miRNA.bed" "$input_dir/lincRNA_TUCP.bed" "$input_dir/lincRNA_RNA-Seq.bed" "$input_dir/tRNA_Genes.bed" > "$bed_file_prefix+vneg.bed"  # Intersect with multiple RNA feature files
    
    bedtools subtract -A -a "$bed_file_prefix+v4.bed" -b "$bed_file_prefix+vneg.bed" > "$bed_file_prefix+v5.bed"  # Subtract intersected features
    
    bedtools subtract -A -a "$bed_file_prefix+v5.bed" -b "$input_dir/GRCh38_SimpleRepeat_homopolymer_4to6_slop5.bed" "$input_dir/GRCh38_SimpleRepeat_homopolymer_7to11_slop5.bed" "$input_dir/GRCh38_SimpleRepeat_homopolymer_gt11_slop5.bed" > "$bed_file_prefix+v6.bed"  # Subtract homopolymer repeat regions
    
    bedtools intersect -wb -sorted -a "$bed_file_prefix+v6.bed" -b "$bed_file_prefix.bed" > "$bed_file_prefix+v7.bed"  # Intersect with original .bed file
    
    awk '{ $4=""; $5=""; $6=""; $7=""; print $0 }' "$bed_file_prefix+v7.bed" | awk '{$1=$1; print}' > "$bed_file_prefix+v8.bed"  # Clean up unnecessary columns
    
    sed 's/ \{1,\}/\t/g' "$bed_file_prefix+v8.bed" > "$bed_file_prefix+v9.bed"  # Replace multiple spaces with tabs
    
    sort -k1,1 -k2,2n "$bed_file_prefix+v9.bed" > "$bed_file_prefix+v10.bed"  # Sort the .bed file by chromosome and start position
    
    bedtools merge -d 50 -c 4 -o mean -i "$bed_file_prefix+v10.bed" > "$bed_file_prefix+v11.bed"  # Merge nearby intervals and calculate mean
    
    bedtools subtract -A -a "$bed_file_prefix+v11.bed" -b "$input_dir/hg38.knownGene.gtf" > "$bed_file_prefix+v12.bed"  # Subtract known gene regions
    
    bedtools subtract -A -a "$bed_file_prefix+v12.bed" -b "$input_dir/generalhg38.bed" > "$bed_file_prefix+v13.bed"  # Subtract general feature regions
    
    bedtools subtract -A -a "$bed_file_prefix+v13.bed" -b "$input_dir/mt_hg38max.bed" > "$bed_file_prefix+v14.bed"  # Subtract mitochondrial regions

    # Add a new value in the fifth column for identification purposes.
    
    awk -v new_val="$base_name" 'BEGIN {FS=OFS="\t"} {$5=new_val; print}'  $bed_file_prefix+v14.bed > $bed_file_prefix+v15.bed
    
    echo "$base_name finished"  # Notify that processing for the current file is complete
done

# Combine all processed .bed files into a single file.
cat "$output_dir"/*+v15.bed > "$output_dir"/combined_file.bed

# Sort the combined .bed file.
bedtools sort -i "$output_dir"/combined_file.bed > "$output_dir"/combined_file+1.bed

# Merge entries in the sorted combined .bed file while calculating mean and distinct values.
bedtools merge -c 4,5 -o mean,distinct -i "$output_dir"/combined_file+1.bed > "$output_dir"/combined_final.bed
