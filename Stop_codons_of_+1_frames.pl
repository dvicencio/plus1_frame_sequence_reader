#!/usr/bin/perl -w
open(SEQFILE, "ScORFs.txt")||die "opening file $!";
@ORFarray = <SEQFILE>;
close (SEQFILE);
#this segment of code reads each line of the file, defining ORFs details, into an array


@NEWDATA=();
open (RESULTS, ">>practice1.txt") ||die "cannot open results.txt: $!";
@NEWDATA = <RESULTS>;




push (@NEWDATA, "Normal Stopcodons in ORFs\n");



for($index=0; $index<@ORFarray; $index++){
                
    $gene = $ORFarray [$index];
    #this 'for' loop takes each line of sequence details in turn out of the array ready for processing
    
    $findtext = index ($gene, ">" , 0);
    # finds out where the letter Y is, defining the beginning of the Scer yeast name
    
    $scername= substr ($gene,$findtext,8);
    # extracts the Scer gene name 
    
    $ATGregion = index ($gene, "ATG", 0);
    
    $ORFseq = substr ($gene, $ATGregion+1);
    
    $genelen = length ($ORFseq) -2;
    #measures length of gene sequence
   
    print "$scername \t  $genelen \n";
    
my $len = 3;

    $gene = $ORFseq;

for (my $ORFcod = 1; $ORFcod <= length $ORFseq; $ORFcod += ($len)) {
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
  
  
  
  
   $codon = substr ($ORFseq, $ORFcod - 1, $len);
    
    if ($codon eq $stopsite1 ){   
        push (@NEWDATA,"$codon");
        print "$codon\n";
     }
     elsif ($codon eq $stopsite2 ){   
        push (@NEWDATA,"$codon");
        print "$codon\n";
     }
     elsif ($codon eq $stopsite3 ){   
        push (@NEWDATA,"$codon");
        print "$codon\n";
     }
     
       }
 

 }

     

print RESULTS @NEWDATA;
print scalar(@NEWDATA)-1;
close(RESULTS);

