# GPXplore
Gene Promoter Extraction Tool
GPXplore is an advanced computational tool designed to efficiently and precisely extract gene promoters from genomic data. Leveraging a user-friendly interface and robust algorithms, GPXplore facilitates the retrieval of upstream and downstream sequences relative to gene loci, which are critical for understanding gene regulation

### Input Files required to process GPXplore
1. Gene list (fatsa)
2. Annotation file (GFF)
3. Reference Genome File (Fatsa)

## Dual Interface Application 
GPXplore has been developed to cater to users with varying levels of technical expertise, offering both a command-line interface and a graphical user interface (GUI)

## Download the software of your choice

### GPXplore GUI for Windows 10 and Later: 
      https://drive.google.com/file/d/1iaM40kyXVcZGYCa5ysvV-qKgJ7vPqfKc/view?usp=drive_link
      
Download the file using the provided link, open the application by double-clicking, browse for the input files and output directory, set the length for upstream and downstream sequence length and click the "Extract Promoter" button.

### GPXplore GUI application for ubuntu 21.10 and above: 
      https://drive.google.com/file/d/1laAXVWQ29nFGOZ5a6Bc-XeU-pXQBg08l/view?usp=drive_link
      
 Download the file using the provided link, open the application by double-clicking, browse for the input files and output directory, set the length for upstream and downstream sequence length and click the "Extract Promoter" button.  
 
### GPXplore GUI python code:
     https://drive.google.com/file/d/1bZ6hugE0JrACsgDBGB9MivPZv3XQ9lDr/view?usp=drive_link
     
Download the python script using the provided link, open the terminal and use the command to open the file 

      python GPXplore_GUI.py 
      
This command open the GUI of GPXplore
 
### GPXplore for Command line interface
     https://drive.google.com/file/d/1zKIogkPLg4tE0GLj4SIbf2GriUbrYYeN/view?usp=drive_link
     
Download the script using the provided link, open the terminal and use the command to run the GPXplore

      python GPXplore_terminal.py --genome_path <reference genome.fasta> --gene_path <gene.fasta> --gff_path <annotation gff file> --output_path output file --upstream_len 100 --downstream_len 100

--genome_path : reference genome file in fasta format

--gene_path : gene sequences for which you want to extract promoter  (this file can have one or more gene sequnce)

--gff_path : Annotation file in gff format

--output_path : output file location where you want to store your results 

--upstream_len 100 : length of upstream sequence 

--downstream_len : length of downstream sequence

      
## Developed by:
### Dr. Shbana Begam, Scientist
#### ICAR-National Institute for Plant Biotechnology, New Delhi-110012
#### Contact: shbana.begam@icar.gov.in

### Dr. Samarth Godara, Scientist
#### Division of Computer Applications, ICAR-Indian Agricultural Statistics Research Institute (IASRI), New Delhi-110012 
#### Contact: samarth.godara@icar.gov.in
