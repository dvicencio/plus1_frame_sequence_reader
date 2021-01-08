#!/usr/bin/perl -w
open(SEQFILE, "ScORFs.txt")||die "opening file $!";
@ORFarray = <SEQFILE>;
close (SEQFILE);
#this segment of code reads each line of the file, defining ORFs details, into an array


@NEWDATA=();
open (RESULTS, ">>numstopcodonsafter+1site.txt") ||die "cannot open results.txt: $!";
@NEWDATA = <RESULTS>;




push (@NEWDATA, "stopcodons\n");


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
   
    print "$scername \t  $genelen \n";
    
my $len = 3;

    $gene = $ORFseq;

for (my $ORFcod = 1; $ORFcod <= length $ORFseq; $ORFcod += ($len)) {
  
   $fsitectt1 = "CTTACG";
   $fsitectt2 = "CTTTCG";
   $fsitectt3 = "CTTCGG";
   $fsitectt4 = "CTTAGG";
   $fsitectt5 = "CTTCTC";
   $fsitectt6 = "CTTCAG";
   
   $fsiteggg1 = "GGGACG";
   $fsiteggg2 = "GGGTCG";
   $fsiteggg3 = "GGGCGG";
   $fsiteggg4 = "GGGAGG";
   $fsiteggg5 = "GGGCTC";
   $fsiteggg6 = "GGGCAG";
   
   $fsitegcg1 = "GCGACG";
   $fsitegcg2 = "GCGTCG";
   $fsitegcg3 = "GCGCGG";
   $fsitegcg4 = "GCGAGG";
   $fsitegcg5 = "GCGCTC";
   $fsitegcg6 = "GCGCAG";
   
   $fsiteccg1 = "CCGACG";
   $fsiteccg2 = "CCGTCG";
   $fsiteccg3 = "CCGCGG";
   $fsiteccg4 = "CCGAGG";
   $fsiteccg5 = "CCGCTC";
   $fsiteccg6 = "CCGCAG";
   
   $fsitecuc1 = "CTCACG";
   $fsitecuc2 = "CTCTCG";
   $fsitecuc3 = "CTCCGG";
   $fsitecuc4 = "CTCAGG";
   $fsitecuc5 = "CTCCTC";
   $fsitecuc6 = "CTCCAG";
   
   $fsiteagg1 = "AGGACG";
   $fsiteagg2 = "AGGTCG";
   $fsiteagg3 = "AGGCGG";
   $fsiteagg4 = "AGGAGG";
   $fsiteagg5 = "AGGCTC";
   $fsiteagg6 = "AGGCAG";
   
   
   $fsitecug1 = "CTGACG";
   $fsitecug2 = "CTGTCG";
   $fsitecug3 = "CTGCGG";
   $fsitecug4 = "CTGAGG";
   $fsitecug5 = "CTGCTC";
   $fsitecug6 = "CTGCAG";
   
   

    my $codon = substr ($ORFseq, $ORFcod - 1, $len);
    
     $sixnt = substr ($ORFseq, $ORFcod - 1, $len +3);
    
     $position = ($ORFcod + (length $sixnt) - 1);
     %pos = ($sixnt => $position);
     
      
 if ($sixnt eq $fsitectt1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
 }
 if ($sixnt eq $fsitectt2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
  $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 
 push (@NEWDATA,"\n")
 }
 
 if ($sixnt eq $fsitectt3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
  
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 if ($sixnt eq $fsitectt4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 if ($sixnt eq $fsitectt5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
  
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 if ($sixnt eq $fsitectt6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 if ($sixnt eq $fsiteggg1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
  
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 if ($sixnt eq $fsiteggg2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 if ($sixnt eq $fsiteggg3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 if ($sixnt eq $fsiteggg4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 if ($sixnt eq $fsiteggg5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 if ($sixnt eq $fsiteggg6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 if ($sixnt eq $fsitegcg1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 if ($sixnt eq $fsitegcg2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 if ($sixnt eq $fsitegcg3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 if ($sixnt eq $fsitegcg4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 if ($sixnt eq $fsitegcg5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 if ($sixnt eq $fsitegcg6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 if ($sixnt eq $fsiteccg1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 if ($sixnt eq $fsiteccg2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 if ($sixnt eq $fsiteccg3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 if ($sixnt eq $fsiteccg4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 if ($sixnt eq $fsiteccg5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 if ($sixnt eq $fsiteccg6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 if ($sixnt eq $fsitecuc1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 if ($sixnt eq $fsitecuc2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecuc3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecuc4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecuc5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
  
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecuc6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 if ($sixnt eq $fsiteagg1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 if ($sixnt eq $fsiteagg2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 if ($sixnt eq $fsiteagg3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 if ($sixnt eq $fsiteagg4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 if ($sixnt eq $fsiteagg5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 if ($sixnt eq $fsiteagg6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecug1) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecug2) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 if ($sixnt eq $fsitecug3) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecug4) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecug5) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 push (@NEWDATA,"\n");
  
 }
 
 
 
 
 
 
 
 
 if ($sixnt eq $fsitecug6) {
  
   $newstart = index ($gene, "???", 0);
    
    $newseq = substr ($gene, $newstart+($pos{$sixnt})-1);
  
   
  
        $newseqlen = length($newseq);
 
 
 for ($stop =0; $stop <= $newseqlen; $stop = $stop += ($len)){
  
    $stopsite1 = "TAA";
    $stopsite2 = "TAG";
    $stopsite3 = "TGA";
    
     $newcodon = substr($newseq, $stop, $len);
    
     
 
    if ($newcodon eq $stopsite1 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite2 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     elsif ($newcodon eq $stopsite3 ){   
        push (@NEWDATA,"$newcodon\t")
     }
     
       }
 
  
    }

 
}
   
         
 
}

print scalar @NEWDATA;
print RESULTS @NEWDATA;
close(RESULTS);


