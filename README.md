# +1 frame sequence reader
+1 frame sequence reader is a Perl program to read through a nucleotide sequence as a [ribosome](https://www.genome.gov/genetics-glossary/Ribosome) would do, in 3-nuclotide ([codon](https://www.genome.gov/genetics-glossary/Codon#:~:text=A%20codon%20is%20a%20trinucleotide,in%20groups%20of%20three%20bases.)) steps.
# Description
+1 frame sequence reader is a Perl program consisting of three [loops](https://www.perltutorial.org/perl-for-loop/) (subsequent loops inside the previous one), two '[if statements](https://www.perltutorial.org/perl-if/),' and two '[elsif statements](https://www.perltutorial.org/perl-if/)' to spot sequences in the +1-frame of a genome file.
* Loop 1: places all gene sequences into a separate [array](https://www.perltutorial.org/perl-array/).

* Loop 2: inside loop 1, lists all nucleotides of the gene in groups of three also known as codons. This way the loop is also commanded to read the sequence in `codon + codon` steps just as a ribosome would move through an mRNA. 
   * If statement 1: inside loop2, to search for sequences with the `6` nucleotides (di-codon) of interest (this could be modified to more or less nucleotides depending of the length of sequence searched) and to retrieve the +1-frame sequence downstream of the di-codon.
   
* Loop 3: inside loop 2 and "if statement 1", keeps reading through the gene squence, although in the +1-frame now, listing the nucleotides in codon groups once again. 
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
# Usage
## Initiation
Selects the file with the genome sequence of interest (e.g. s. cerevisiae ORFs in this case) and defines it as an array containing genes as elements
```
#!/usr/bin/perl -w
open(SEQFILE, "file_with_genome.txt")||die "opening file $!";
@ORFarray = <SEQFILE>; # first, we define our sequence file as an array
close (SEQFILE);
# this segment of code reads each line of the file, defining genes details, into an array

@NEWDATA=();
open (RESULTS, ">>file_with_results.txt") ||die "cannot open results.txt: $!";
@NEWDATA = <RESULTS>;
# the 3 previous lines create a new file to deliver results

push (@NEWDATA, "frame-codon2.pl\n");
push (@NEWDATA, "Gene name\t");
push (@NEWDATA, "Gene length\t");
push (@NEWDATA, "Nucleotide position\t");
push (@NEWDATA, "Potential frameshift bi-codon\t");
push (@NEWDATA, "Hypothetical +1 sequence\t");
push (@NEWDATA, "Stop codon position and seq downstream\n");
# the 8 lines create labels in the new file created
```

## Loop 1
For each gene information, list their letters (nucleotides) into an array
```
for($index=0; $index<@ORFarray; $index++){
# this line defines the length of the array, or gene, in unit "nucleotide" elements             
    $gene = $ORFarray [$index];
    # this line defines each gene (ORF) as the sequence read until a line break is found
    
    $findtext = index ($gene, ">" , 0);
    # finds out where the symbol ">" is, defining the beginning of the S. cerevisiae gene name
    
    $scername= substr ($gene,$findtext,8);
    # extracts eight letters of the S. cerevisiae gene name
    
    $ATGregion = index ($gene, "???", 0);
    # dentifies the beginning of the ORF by indicating the "???" characters situated before ATG start codons
    
    $ORFseq = substr ($gene, $ATGregion+3);
    # extracts the ORF sequence from the start codon to the stop codon
    
    $genelen = length ($ORFseq) -2;
    # measures the number of nucleotides in the ORF
   
    print "$scername \t  $genelen \n";
    
my $len = 3;
# defines 3 nucleotides length used in loop 2
    $gene = $ORFseq; # Specifies that the gene is equivalent to the ORF
```
## Loop 2
For each ORF, list their nucleotides in sets of three reading the sequence from the start codon to the end of the ORF in codon steps
```
for (my $ORFcod = 1; $ORFcod <= length $ORFseq; $ORFcod += ($len)) {
# this line defines the length of the array beggining from "ATG" to the end of the ORF in 3 nucleotide steps 

my $codon = substr ($ORFseq, $ORFcod - 1, $len);
# extracts the codons from the ORF in order until the end

     $sixnt = substr ($ORFseq, $ORFcod - 1, $len +3);
     # extracts the codons plus the following 3 nucleotides (di-codon) for each codon from the beginning of the ORF to the end. 
     
     $position = ($ORFcod + (length $sixnt) - 1);
     # extracts the range position for each di-codon in the ORF (e.g. position 206)
     
     %pos = ($sixnt => $position);
     # generates key-value pairs for each di-codon => position to retrieve them when needed
```
### If Statement 1
Looks for specific di-codon sequences and pushes all the related information into the new file
```
 if ($sixnt eq $fsitectt1) {
 # As Loop 2 reads through each codon + codon in the sequence if a di-codon is equal to $fsitectt1 (a string which is CTTACG in this case) the following code is applied:
 
   $newstart = index ($gene, "???", 0);
   # finds out the beginning of the ORF to map out the positions of each di-codon
    
   $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
   # finds the di-codon and its exact position in the ORF using the key-value previously defined in Loop 2 and retrieves the +1 frame downstream sequence by skipping one nucleotide. For instance, when the program finds the di-codon CTTACG, it will retrieve the new sequence strating from CGX.
  
        push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
        # the 5 previous lines retrieve all the information related to the di-codon of interest which includes: gene name, gene length, di-codon position within ORF, and the +1-frame sequence downstream
        
  $newseqlen = length($newseq);
  # defines the new frame sequence length for future reference
  ```
## Loop 3
for each +1 sequence, list the codons until the end of the ORF
```
for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
# this line, once again, defines the length of the array; however, it starts reading from the +1 frame of the sequence downstream the di-codon starting from the fourth nucleotide from left to right.  

    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    # the 3 previous lines define the stop codons as a string
    
     $newcodon = substr($newseq, $stop, $len);
     # extracts the codons from the +1 frame sequences in order
     
     $newposition = ($stop + (length $newcodon));
     #defines the position of each codon in the +1 frame sequence
     
    %newpos = ($newcodon => $newposition);
    # generates key-value pairs for each codon => +1 frame position to retrieve them when needed
```
### If statement 2
Looks for the first immediate stop codon ( in case TAA is the first) after the di-codon +1 frame initiates and pushes the in between sequence into the new file
 ```
if ($newcodon eq $stopsite1 ){
    
    $stopsite = index ($newcodon, $stopsite1);
    # identifies the first stop codon in the +1 frame sequence
    
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    # extracts the new sequence with the identified stop codon at the end
    
    $stop1 = $newposition;
    # identifies position of the stop codon in the +1 frame sequence
    
    push (@a,"$stop1\t $stopseq \n");
    # pushes the new +1 frame sequence with the first stop codon encountered downstream and its position into a new array
  ```
 ### Elsif statement 1
 

Looks for the first immediate stop codon ( in case TAG is the first) after the di-codon +1 frame initiates and pushes the in between sequence into the new file
    ```
    $stopsite = index ($newcodon, $stopsite2);
    # identifies the first stop codon in the +1 frame sequence
  
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    # extracts the new sequence with the identified stop codon at the end
    
    $stop1 = $newposition;
    # identifies position of the stop codon in the +1 frame sequence
    
    push (@a, "$stop1\t $stopseq \n");
    # pushes the new +1 frame sequence with the first stop codon encountered downstream and its position into a new array
    ```

