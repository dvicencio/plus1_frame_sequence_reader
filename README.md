# +1 frame sequence reader
+1 frame sequence reader is a Perl program to read through a nucleotide sequence as a [ribosome](https://www.genome.gov/genetics-glossary/Ribosome) would do, in 3-nuclotide ([codon](https://www.genome.gov/genetics-glossary/Codon#:~:text=A%20codon%20is%20a%20trinucleotide,in%20groups%20of%20three%20bases.)) steps.
# Description
+1 frame sequence reader is a Perl program consisting of three [loops](https://www.perltutorial.org/perl-for-loop/) (subsequent loops inside the previous one), two '[if statements](https://www.perltutorial.org/perl-if/),' and two '[elsif statements](https://www.perltutorial.org/perl-if/)' to spot sequences in the +1-frame of a genome file.
* Loop 1: places all gene sequences into an [array](https://www.perltutorial.org/perl-array/).

* Loop 2: inside loop 1, lists all nucleotides of the gene in groups of three also known as codons. This way the loop is also commanded to read the sequence in `codon + codon` steps just as a ribosome would move through an mRNA. 
   * If statement 1: inside loop2, to search for sequences with the `6` nucleotides (di-codon) of interest (this could be modified to more or less nucleotides depending of the length of sequence searched) and to retrieve the +1-frame sequence downstream of the di-codon.
   
* Loop 3: inside loop 2, keeps reading through the gene squence, although in the +1-frame now, listing the nucleotides in codon groups once again. 
   * If satement 2: inside Loop 3, to search for the stop codon `TAA` and , if this is the first stop codon identified downstream the +1 frame after the di-codon, push the sequence in between to a new file. 
   * Elsif 1: inside loop 3, to search for the stop codon `TAG` and, if this is the first stop codon identified dowsntream the +1 frame after the di-codon, push the sequence in between to a new file.
   * Elsif 2: inside loop 3, to search for the stop codon `TGA` and, if this is the first stop codon identified downstream the +1 frame after the di-codon, push the sequence in between to a new file.
   
Each loop can be manipulated to retrieve and quantify specific information like: 
* total number of codons 
* total number of stop codons 
* total number of di-codons 
* total number of genes 
* total number of nucleotides in a gene 
* the numerical position of any sequence within the genome    
