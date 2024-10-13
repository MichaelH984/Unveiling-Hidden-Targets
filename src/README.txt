Scripts README


Converterfinal.sh Script

Overview
The Converterfinal script is designed to automate the conversion of .bw files to .wig and subsequently to .bed files, while performing multiple transformations throughout the process.

Dependencies
- bigtowig is a command-line utility that converts BigBed files into Wiggle format. This format is often used in genomic data visualization and allows for easy interpretation of coverage and signal data (This should be installed in the src folder).
- The only external tool required for this script is bedtools, which can be installed using the following command:
  
	$sudo apt-get install bedtools

Usage
To run the script, you must start from the src directory where the Converterfinal script is located. Ensure that the working directory contains the necessary input files and directories. The .bw files should be placed in the ../Raw data directory, while output files will be saved in ../bin.

Steps:
1. Navigate to the src directory:
   cd src
2. Run the script:
   bash Converterfinal.sh

Input:
- The script reads .bw (BigWig) files from the ../Raw data directory.

Output:
- The resulting .wig and .bed files will be stored in the ../bin directory.

Notes:
- Ensure that all necessary files are placed in the correct directories before running the script (Nothing should be changed from downloading from GitHub). The script will handle conversions and any transformations automatically, making it easier to process multiple .bw files in one go.


###########################################################################################################################################################


Distance Calculation and Frequency Analysis Script

## Overview
This script automates the process of analyzing the GWIPS-viz Global aggregate Ribo-seq data by intersecting with CCDS regions, filtering unique entries, and calculating the frequency distributions of scores.

Dependencies
- Ensure bedtools is installed and accessible in your system's PATH.

Steps:
1. Navigate to the src directory:
   cd src
2. Run the script:
   bash Distance_calculator.sh

Output: Two files will be available in the output directory:

	FQtable.txt contains a frequency table of all the distances between peaks.

	percentile_scores.txt contains the scores from the 90th to the 100th percentile.


###########################################################################################################################################################


Rplot.txt
## Overview 
This txt file contains the R script used to construct figure 1.

Dependencies
- ggplot2
- dplyr

Notes: Used in R studio, make sure to change lines 4 and 5 to suit the path to your files.
