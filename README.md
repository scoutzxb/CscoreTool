# CscoreTool


1. Purpose

This program is to do A/B compartment analysis for Hi-C data. Test data are also provided.
Several other files are also available for download, including the bedgraph files involved in the CscoreTool paper, the generateEqualLengthBed.cpp file, and the chromosome size files used for generating 

2. Installation

The program is simple one-file executable for linux system. Just download CscoreTool1.1 Then run the following command:
  
  chmod +x CscoreTool1.1

Then you can type ./CscoreTool1.1 to run it.

The executable is compiled on a RedHat linux machine. If it doesn't work, download the CscoreTool1.1.cpp, twister.h and twister.cpp files and put them in the same folder. Then run

  g++ CscoreTool1.1.cpp twister.cpp -fopenmp -O3 -o CscoreTool1.1

You'll get an execuable file CscoreTool1.1. 

3. Usage: CscoreTool1.1 < windows.bed> < input.summary> < outputPrefix> < session> < minDis> [chrName]

Input parameters

a. windows.bed 
This file is to specify the genomic windows to analyze. It should be equal-length bed files covering the region of interest, presumably a whole chromosome or whole genome. An example hg19_chr1_10k.bed can be downloaded.These files can also be generatd using the generateEqualLengthBed.cpp program provided here. The chromosome size files used for generating windows.bed are also available for download.

b. input.summary
This file is the main input file for Hi-C interactions. We accept the same format as the HiCsummary file format for HOMER runHiCpca.pl. See http://homer.ucsd.edu/homer/interactions/HiCtagDirectory.html
An example file test.summary.gz can be downloaded. This is 0.5% randomly selected reads in chr1 from the High-resolution GM12878 cell Hi-C dataset (Rao, 2014).  

c. OutputPrefix
This is the prefix for output files.

d. session
This the number of sessions to use. The number of choide depends on the resource available. The program is not perfectly parallelized, so using more sessions just partly improve the speed. 

e. minDis
This is the criteria of minimum interaction distance to be considered. It should be no less than the analysis resolution. We suggest using 1000000 (1M), because interactions of shorter distance may be more affected by TAD structure than A/B compartments.

f. chrName
This is a parameter by choice. Users can speficy a certain chromosome to analyze using it.

Example run:

CscoreTool1.1 hg19_chr1_10k.bed Test.summary Test_chr1_10k 12 1000000

4. Output

There are 4 output files.

XXX_bias.txt 
This is the bias factor estimated for each genomic window. 
XXX_hh.txt 
This is the estimated distance curve. The corresponding distance is 10^(0.04*n), where n is the number in the first column.
XXX_cscore.txt 
This is the Cscore estimated for each genomic window. 
XXX_cscore.bedgraph
This is the bedgraph file made for visualization. Low-mappability windows (bias<0.2) or high copy-number windows (bias>5) are already filtered. 
All bedgraph files involved in the analysis in the CscoreTool paper are available in the Bedgraphs.zip file.
