#!/bin/bash

########################################################################
# Hope4Gene                                   
#   Hopefield-like class prediction algorithm for gene expression data (Hope4Gene) for >=2 classes
#   Cantini Laura (Univeristà degli studi di Torino)
#
#    input: (1) Gene expression matrix in log2ratio
#                  1st column: Gene name
#                  1st row: samples
#           (2) Marker genes matrix
#                  1st column: Gene name
#                  2nd column: Class (1,2,...)
#           (3) res.number number of signature resamplings only if the number of classes is greater than 2
#           (4) gene.res number of signature genes to be resampled only if the number of classes is greater than 2
#           (5) permutations number number of permutations for the FDR computation
#                  
#  
#    output: classification report
#            
########################################################################
#########################################################

echo -n "Enter dataset file name > "
read text1
echo "You entered: $text1"
echo -n "Enter signature file name > "
read text2
echo "You entered: $text2"
echo -n "Enter total class number > "
read text3
echo "You entered: $text3"
if [ $text3 -ge 3 ]; then
	echo -n "Enter res.number > "
	read text4
	echo "You entered: $text4"
	echo -n "Enter gene.res > "
	read text5
	echo "You entered: $text5"
fi
echo -n "Enter permutations number  > "
read text6
echo "You entered: $text6"
if [ $text3 -ge 3 ]; then
	Rscript Hope4Gene.R $text1 $text2 $text3 $text4 $text5 >& output1 &
	Rscript permutations.R $text1 $text2 $text3 $text5 $text6 >& output2 &
	
else
	Rscript Hope4Gene.R $text1 $text2 $text3 >& output1 &
	Rscript permutations.R $text1 $text2 $text3 $text6 >& output2 &
	
fi	


wait
Rscript FDR.R $text3 > output3 

