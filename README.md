# +1 frame sequence reader
+1 frame sequence reader is a Perl program to read through a nucleotide sequence as a [ribosome](https://www.genome.gov/genetics-glossary/Ribosome) would do, in 3-nuclotide ([codon](https://www.genome.gov/genetics-glossary/Codon#:~:text=A%20codon%20is%20a%20trinucleotide,in%20groups%20of%20three%20bases.)) steps.
# Description
+1 frame sequence reader is a Perl program consisting of three [loops](https://www.perltutorial.org/perl-for-loop/) (subsequent loops inside the previous one), two '[if statements](https://www.perltutorial.org/perl-if/),' and two '[elsif statements](https://www.perltutorial.org/perl-if/)' to spot sequences in the +1-frame of a genome file.
* Loop 1: places all gene sequences into a separate [array](https://www.perltutorial.org/perl-array/).

* Loop 2: inside loop 1, lists all nucleotides of the gene in groups of three (also known as codons). This way the loop is also commanded to read the sequence in `codon + codon` steps just as a ribosome would move through an mRNA. 
   * If statement 1: inside loop2, to search for sequences with the `6` nucleotides (di-codon) of interest (this could be modified for more or less nucleotides depending of the length and specifications of sequence of interest) and to retrieve the +1-frame sequence downstream of the di-codon.
   
* Loop 3: inside loop 2 and "if statement 1", keeps reading through the gene squence, although in the +1-frame now, listing the nucleotides in codon groups once again. 
   * If satement 2: inside Loop 3, to search for the stop codon `TAA` and , if this is the first stop codon identified downstream the +1 frame after the di-codon, pushes the sequence in between to a new file. 
   * Elsif 1: inside loop 3, to search for the stop codon `TAG` and, if this is the first stop codon identified dowsntream the +1 frame after the di-codon, pushes the sequence in between to a new file.
   * Elsif 2: inside loop 3, to search for the stop codon `TGA` and, if this is the first stop codon identified downstream the +1 frame after the di-codon, pushes the sequence in between to a new file.
   
Each loop can be manipulated to retrieve and quantify specific information like: 
* total number of codons 
* total number of stop codons 
* total number of di-codons 
* total number of genes 
* total number of nucleotides in a gene 
* the numerical position of any sequence within the genome    
# Usage
## Initiation
Selects the file with the genome sequence of interest (e.g. *S. cerevisiae* ORFs in this case) and defines it as an array containing genes as elements
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
# the 8 previous lines create labels in the new file created
```

## Loop 1
For each gene sequence information, list their letters (nucleotides) into an array to quantify them and identify their locartion in the ORF.
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
    # identifies the beginning of the ORF by indicating the "???" characters situated before the "ATG" start codons
    
    $ORFseq = substr ($gene, $ATGregion+3);
    # extracts the ORF sequence from the start codon to the end
    
    $genelen = length ($ORFseq) -2;
    # measures the number of nucleotides in the ORF
   
    print "$scername \t  $genelen \n";
    
my $len = 3;
# defines 3 nucleotides length used in loop 2
    $gene = $ORFseq; # Specifies that the gene is equivalent to the ORF
```
Example of the information extracted by Loop 1 from the genome file: gene name and length (a shorter genome file was used to exemplify the execution of the code)

![](loop1.gif)
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
     # extracts the position for each di-codon in the ORF (e.g. position 206)
     
     %pos = ($sixnt => $position);
     # generates key-value pairs for each di-codon => position to retrieve them when needed
```
Example of the information extracted by Loop 2 from the genome file: codons and di-codons with total number
![](loop2.gif)
### If Statement 1
Looks for specific di-codon sequences of interest and pushes all the downstream +1 frame sequences into the new file
```
 if ($sixnt eq $fsitectt1) {
 # As Loop 2 reads through each codon + codon in the sequence, if a di-codon is equal to $fsitectt1 (a string which is CTTACG in this case) the following code is applied:
 
   $newstart = index ($gene, "???", 0);
   # finds out the beginning of the ORF to map out the positions of each di-codon
    
   $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
   # finds the di-codon and its exact position in the ORF using the key-value previously defined in Loop 2 and retrieves the +1 frame downstream sequence by skipping one nucleotide. For instance, when the program finds the di-codon CTTACG, it will retrieve the new sequence starting from CGX.
  
        push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
        # the 5 previous lines retrieve all the information related to the di-codon of interest which includes: gene name, gene length, di-codon position within ORF, and the +1-frame sequence downstream
        
  $newseqlen = length($newseq);
  # defines the new frame sequence length for future reference
  ```
  Example of the information extracted by the If statement 1 from the genome file: Specific di-codons of interest and their downstream +1 sequence
  
  ![](ifstate.gif)
## Loop 3
For each +1 sequence, list the codons until the end of the ORF
```
for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
# this line, once again, defines the length of the array; however, it starts reading from the +1 frame of the sequence downstream the di-codon starting from the fourth nucleotide from left to right.  

    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    # the 3 previous lines define the stop codons as strings
    
     $newcodon = substr($newseq, $stop, $len);
     # extracts the codons from the +1 frame sequences in order
     
     $newposition = ($stop + (length $newcodon));
     # defines the position of each codon in the +1 frame sequence
     
    %newpos = ($newcodon => $newposition);
    # generates key-value pairs for each codon => +1 frame position to retrieve them when needed
```
Example of the information extracted by Loop 3 from the genome file: codons in the downstream +1 frame

![](loop3.gif)
### If statement 2
Looks for the first immediate stop codon ( in case that TAA is the first) after the di-codon +1 frame initiates and pushes the "in between sequence" into the new file
 ```
if ($newcodon eq $stopsite1 ){
# As Loop 3 reads through each codon + codon in the + 1 frame sequence if a di-codon is equal to $stopsite1 (TAA) the following code is applied: 

    $stopsite = index ($newcodon, $stopsite1);
    # identifies the first stop codon in the +1 frame sequence
    
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    # extracts the new sequence with the identified stop codon at the end
    
    $stop1 = $newposition;
    # identifies position of the stop codon in the +1 frame sequence
    
    push (@a,"$stop1\t $stopseq \n");
    # pushes the new +1 frame sequence with the first stop codon encountered downstream and its position into a new array
    
    } # end of if statement 2
  ```
### Elsif statement 1
Looks for the first immediate stop codon ( in case that TAG is the first) after the di-codon +1 frame initiates and pushes the "in between sequence" into the new file
    
```
elsif ($newcodon eq $stopsite2) {
# As Loop 3 reads through each codon + codon in the + 1 frame sequence if a di-codon is equal to $stopsite1 (TAG) the following code is applied:
    
    $stopsite = index ($newcodon, $stopsite2);
    # identifies the first stop codon in the +1 frame sequence
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    # extracts the new sequence with the identified stop codon at the end
    
    $stop1 = $newposition;
    # identifies position of the stop codon in the +1 frame sequence
    
    push (@a, "$stop1\t $stopseq \n");
    # pushes the new +1 frame sequence with the first stop codon encountered downstream and its position into a new array
    
    } # end of elsif statement 1
```
### Elsif statement 2
Looks for the first immediate stop codon ( in case that TGA is the first) after the di-codon +1 frame initiates and pushes the "in between sequence" into the new file
```
elsif ($newcodon eq $stopsite3) {
# As Loop 3 reads through each codon + codon in the + 1 frame sequence if a di-codon is equal to $stopsite1 (TAG) the following code is applied:   
      
    $stopsite = index ($newcodon, $stopsite3);
    # identifies the first stop codon in the +1 frame sequence
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    # extracts the new sequence with the identified stop codon at the end
    
    $stop1 = $newposition;
    # identifies position of the stop codon in the +1 frame sequence
    
   push (@a,"$stop1\t $stopseq \n");
   # pushes the new +1 frame sequence with the first stop codon encountered downstream and its position into a new array  
     
     } # end of elsif statement 2
     
 } # end of loop 3

    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
    # Once the code has read through all the +1 sequence, it pushes the first sequence (the +1 frame sequnece with the first stop codon identified) in the @a array into the new file (@NEWDATA;results)
     
    } # end of if statement 1
     
   shift @a; # after finish reading through the di-codon indicated in statement 1, the code deletes the @a array to avoid overstacking of sequences and allow the next +1 sequence with the first stop codon to take its place
 ```
 Example of the information extracted by If statement 2 and Elsif statements 1 and 2 from the genome file: +1 sequeces with stop codons to retrieve ORF like sequences
 
 ![](elsifs.gif)
 ## Big Picture
The code starts reading through a gene with loops cycling through the sequence looking for specific di-codons. 

When a di-codon is found, subsequent loops take place to retrieve more information.

Further details are obtained and indicated by the code like the +1 sequences and the immediate stop codons.

When the code has ran through the first gene and the first di-codon, if there is another di-codon of interst in the gene the code applies the same loops and if statements to retrieve its information.

On the other hand, if there are no more di-codons of interest in the gene, the program jumps to the next ORF.

This method is applied to each gene in the Saccharomyces cerevisiae's genome file until the last gene in the file has been read by the program. (For a broad and further understanding of the program please read the [research paper](https://github.com/dvicencio/plus1_frame_sequence_reader/blob/main/Dissertation%20final.docx))

# Genome File Specifications and Edition
Genome files are usually represented in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format). However, slight changes to the file have to be done in order to run the code as intended. Such changes are (the following changes could be done using a text editor or the command line):
* Remove all line breaks in the file
* Add the characters "???" before every start codon
* Add a line break before every ">" symbol to separate each gene

note: some symbols (e.g.">") are found in between the gene information which were edited/removed individually (approx 7)
# Additional Notes
The program was sligthly edited for different circumstances were other informations was required like:
* Total number of stop codons in a +1 (achieved by skipping a nucleotide at the beginning of each gene)
* Total number of codons, di-codons, and ORFs
* Concatenation of all +1 sequences

All codes retrieve the information into a new "Results file"

