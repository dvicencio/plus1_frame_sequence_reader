#!/usr/bin/perl -w
open(SEQFILE, "Queryfile.txt")||die "opening file $!";
@ORFarray = <SEQFILE>;
close (SEQFILE);
#this segment of code reads each line of the file, defining ORFs details, into an array

@NEWDATA=();
open (RESULTS, ">>Resultsfile.txt") ||die "cannot open results.txt: $!";
@NEWDATA = <RESULTS>;



push (@NEWDATA, "frame-codon2.pl\n");

push (@NEWDATA, "Codons\n");



for($index=0; $index<@ORFarray; $index++){
                 
    $gene = $ORFarray [$index];
    #this 'for' loop takes each line of sequence details in turn out of the array ready for processing
    
    $findtext = index ($gene, ">" , 0);
    # finds out where the letter Y is, defining the beginning of the Scer yeast name
    
    $scername= substr ($gene,$findtext,8);
    # extracts the Scer gene name 
    
    $ATGregion = index ($gene, "???", 0);
    
    $ORFseq = substr ($gene, $ATGregion+3);
    
    $genelen = length ($ORFseq) -2;
    #measures length of gene sequence
   
   print "$scername \t $genelen\n";
    
   
my $len = 3;

    $gene = $ORFseq;

for (my $ORFcod = 1; $ORFcod <= length $ORFseq; $ORFcod += ($len)) {
    
     $codon = substr ($ORFseq, $ORFcod - 1, $len);
    
    my $sixnt = substr ($ORFseq, $ORFcod - 1, $len +3);
    
    my $position = ($ORFcod + (length $sixnt) - 1);
    
    push (@NEWDATA, "$codon");
   
   


   }

   


}


print RESULTS @NEWDATA;
print (scalar(@NEWDATA)-scalar(@ORFarray));
close(RESULTS);

