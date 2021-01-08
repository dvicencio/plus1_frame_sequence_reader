#!/usr/bin/perl -w
open(SEQFILE, "file_with_genome.txt")||die "opening file $!";
@ORFarray = <SEQFILE>; # first, we define our sequence file as an array
close (SEQFILE);
# this segment of code reads each line of the file, defining ORFs details, into an array


@NEWDATA=();
open (RESULTS, ">>file_with_results.txt") ||die "cannot open results.txt: $!";
@NEWDATA = <RESULTS>;
# lines 8 to 11 create a new file to deliver results


push (@NEWDATA, "frame-codon2.pl\n");
push (@NEWDATA, "Gene name\t");
push (@NEWDATA, "Gene length\t");
push (@NEWDATA, "Nucleotide position\t");
push (@NEWDATA, "Potential frameshift bi-codon\t");
push (@NEWDATA, "Hypothetical +1 sequence\t");
push (@NEWDATA, "Stop codon position and seq downstream\n");
# lines 14 to 21 create labels in the new file created

for($index=0; $index<@ORFarray; $index++){
# this line define the length of the array, or ORFs, in unit "nucleotide" elements                
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
    # measures the number of nucleotides in the 
   
    print "$scername \t  $genelen \n";
    
my $len = 3;

    $gene = $ORFseq; # Specifies that the gene is equivalent to the ORF

for (my $ORFcod = 1; $ORFcod <= length $ORFseq; $ORFcod += ($len)) {
# this line defines the length of the array beggining from "ATG" to the end of the ORF in 3 nucleotide steps

   $fsitectt1 = "CTTACG";
   $fsitectt2 = "CTTTCG";
   $fsitectt3 = "CTTCGG";
   $fsitectt4 = "CTTAGG";
   $fsitectt5 = "CTTCTC";
   $fsitectt6 = "CTTCAG";
   $fsitectt7 = "CTTAGT";
   
   $fsiteggg1 = "GGGACG";
   $fsiteggg2 = "GGGTCG";
   $fsiteggg3 = "GGGCGG";
   $fsiteggg4 = "GGGAGG";
   $fsiteggg5 = "GGGCTC";
   $fsiteggg6 = "GGGCAG";
   $fsiteggg7 = "GGGAGT";
   
   $fsitegcg1 = "GCGACG";
   $fsitegcg2 = "GCGTCG";
   $fsitegcg3 = "GCGCGG";
   $fsitegcg4 = "GCGAGG";
   $fsitegcg5 = "GCGCTC";
   $fsitegcg6 = "GCGCAG";
   $fsitegcg7 = "GCGAGT";
   
   $fsiteccg1 = "CCGACG";
   $fsiteccg2 = "CCGTCG";
   $fsiteccg3 = "CCGCGG";
   $fsiteccg4 = "CCGAGG";
   $fsiteccg5 = "CCGCTC";
   $fsiteccg6 = "CCGCAG";
   $fsiteccg7 = "CCGAGT";
   
   $fsitecuc1 = "CTCACG";
   $fsitecuc2 = "CTCTCG";
   $fsitecuc3 = "CTCCGG";
   $fsitecuc4 = "CTCAGG";
   $fsitecuc5 = "CTCCTC";
   $fsitecuc6 = "CTCCAG";
   $fsitecuc7 = "CTCAGT";
   
   $fsiteagg1 = "AGGACG";
   $fsiteagg2 = "AGGTCG";
   $fsiteagg3 = "AGGCGG";
   $fsiteagg4 = "AGGAGG";
   $fsiteagg5 = "AGGCTC";
   $fsiteagg6 = "AGGCAG";
   $fsiteagg7 = "AGGAGT";
   
   
   $fsitecug1 = "CTGACG";
   $fsitecug2 = "CTGTCG";
   $fsitecug3 = "CTGCGG";
   $fsitecug4 = "CTGAGG";
   $fsitecug5 = "CTGCTC";
   $fsitecug6 = "CTGCAG";
   $fsitecug7 = "CTGAGT";
   
    my $codon = substr ($ORFseq, $ORFcod - 1, $len);
    # extracts the codons from the ORF in order until the end
    
     $sixnt = substr ($ORFseq, $ORFcod - 1, $len +3);
     # extracts the codons plus the following 3 nucleotides (di-codon) for each codon from the beginning of the ORF to the end
     
     $position = ($ORFcod + (length $sixnt) - 1);
     # extracts the range position for each di-codon in the ORF
     
     %pos = ($sixnt => $position);
     # generates key-value pairs for each di-codon => position to retrieve them when needed
      
 if ($sixnt eq $fsitectt1) {
  # As Loop 2 reads through each codon + codon in the sequence if a di-codons is equal to $fsitectt1 (which is CTTACG in this case) the following code is applied: 
  
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
 
    if ($newcodon eq $stopsite1 ){
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
  
 
 
 
 if ($sixnt eq $fsitectt2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 if ($sixnt eq $fsitectt3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
  
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 if ($sixnt eq $fsitectt4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 if ($sixnt eq $fsitectt5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
  
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 if ($sixnt eq $fsitectt6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
   
   
   
   if ($sixnt eq $fsitectt7) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 if ($sixnt eq $fsiteggg1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
  
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 if ($sixnt eq $fsiteggg2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 if ($sixnt eq $fsiteggg3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 if ($sixnt eq $fsiteggg4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 if ($sixnt eq $fsiteggg5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 if ($sixnt eq $fsiteggg6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
   
   
   if ($sixnt eq $fsiteggg7) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 if ($sixnt eq $fsitegcg1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 if ($sixnt eq $fsitegcg2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 if ($sixnt eq $fsitegcg3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 if ($sixnt eq $fsitegcg4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 if ($sixnt eq $fsitegcg5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 if ($sixnt eq $fsitegcg6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 if ($sixnt eq $fsitegcg7) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 if ($sixnt eq $fsiteccg1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 if ($sixnt eq $fsiteccg2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 if ($sixnt eq $fsiteccg3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 if ($sixnt eq $fsiteccg4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 if ($sixnt eq $fsiteccg5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 if ($sixnt eq $fsiteccg6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 if ($sixnt eq $fsiteccg7) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 if ($sixnt eq $fsitecuc1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 if ($sixnt eq $fsitecuc2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecuc3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecuc4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecuc5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
  
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecuc6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 if ($sixnt eq $fsitecuc7) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 if ($sixnt eq $fsiteagg1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 if ($sixnt eq $fsiteagg2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 
 if ($sixnt eq $fsiteagg3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 
 if ($sixnt eq $fsiteagg4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 
 if ($sixnt eq $fsiteagg5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 
 if ($sixnt eq $fsiteagg6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 if ($sixnt eq $fsiteagg7) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 if ($sixnt eq $fsitecug1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecug2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 if ($sixnt eq $fsitecug3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecug4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecug5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
 $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecug6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  
    $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;

 
 
 
 if ($sixnt eq $fsitecug7) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  push (@NEWDATA, "$scername\t");
        push (@NEWDATA, "$genelen\t");
        push (@NEWDATA, "($ORFcod-$position)\t");
        push (@NEWDATA, "$sixnt\t");    
        push (@NEWDATA, "$newseq\t");
  
    $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     $newposition = ($stop + (length $newcodon));
     
    %newpos = ($newcodon => $newposition);
 
 
    if ($newcodon eq $stopsite1 ){   
    
     
     $stopsite = index ($newcodon, $stopsite1);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    push (@a,"$stop1\t $stopseq \n");
    
    
    
     
       
    }
    
    
    
     elsif ($newcodon eq $stopsite2) {
     
      
      $stopsite = index ($newcodon, $stopsite2);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
    $stop1 = $newposition;
    
    
    push (@a, "$stop1\t $stopseq \n");
     
     
     
    }
    
     
     elsif ($newcodon eq $stopsite3) {
      
      
      
      $stopsite = index ($newcodon, $stopsite3);
     
    $stopseq = substr ($newseq,$stopsite,$newpos{$newcodon} );
     $stop1 = $newposition;
     
   push (@a,"$stop1\t $stopseq \n");
     
     
     }
 
 
 
    
 }


    push (@NEWDATA,"Stop sequence\t", $a[0], "\n");
   
  
    
    }
     
   shift @a;

       }
    
}

push (@NEWDATA, "\n");
print RESULTS @NEWDATA;
close(RESULTS);
