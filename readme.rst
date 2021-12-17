===============================
C-WAP
CFSAN Wastewater Analysis Pipeline
===============================

The C-WAP is a bash-based system for the analysis of whole genome sequencing
data of DNA samples extracted from wastewater. While it is possible to adjust
for other organisms of interest, it primarily targets SARS-CoV2 and its variants.

The C-WAP was developed by the United States Food
and Drug Administration, Center for Food Safety and Applied Nutrition.

* Free software: See license below.
* Documentation: ???
* Source Code: ??
* PyPI Distribution: ??


Introduction
------------

The CFSAN Wastewater Analysis Pipeline uses a reference-based alignments to create a matrix of
SNPs for a given set of samples. The process generally starts off by finding
a reference that is appropriate for the samples of interest, and collecting
the sample sequence data into an appropriate directory structure. The SNP
pipeline can then be used to perform the alignment of the samples to the
reference. Once the sample sequences are aligned, a list of SNP positions is
generated. The list of SNP positions is then used in combination with
alignments of the samples to the reference sequence to call SNPs. The SNP
calls are organized into a matrix containing (only) the SNP calls for all
of the sequences.

This software was developed with the objective of creating high quality
SNP matrices for sequences from closely-related pathogens, e.g., different
samples of Salmonella enteriditis from an outbreak investigation. The
focus on closely related sequences means that this code is not suited for
the analysis of relatively distantly related organisms, where there is not
a single reference sequence appropriate for all the organisms for which an
analysis is desired.

The CFSAN Wastewater Analysis Pipeline was written in predominantly in Bash 
with some functionalities implemented in Python with Numpy and Matplotlib. The
code is designed to be straightforward to deploy and run
from the command line. The configuration file (prepareEnvironment.sh) 
provided aims to enumerate the 
dependencies on other tools and the other scripts assume that all referred
executables are available in the search path. URL contained therein 
refer to the repositories that they can be acquired from.
In situations where additional customization is desired, the
code is not highly complex and should be easy to modify as necessary.



Example usage
-------------------
Examples of using the code are provided. These examples serve as both
unit tests, and as examples that can be modified to work on other data
sets of interest.



Citing C-WAP
-------------------
This work is currently unpublished. If you are making use of this package, 
we would appreciate if you gave credit to our repository. 

The html->pdf conversion tool (wkhtmltopdf) is not our work and has been 
provided for user convenience. For more information and updates, please visit 
https://github.com/wkhtmltopdf. 



License
-------

See the LICENSE.txt file included in the C-WAP Pipeline distribution.

