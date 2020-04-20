#! /usr/bin/perl -w
use strict;
#----------------------------------------------------------
#statistics of Multiple sequence alignment of two sequences
#----------------------------------------------------------
die "Usage : perl <VS_SEQ1_SEQ2.ms.fa.xls> <outprefix> \n" unless @ARGV == 2;

my %hash;
open IN,"$ARGV[0]" or die $!;
open POS,">$ARGV[1].pos" or die $!;
open IDL,">$ARGV[1].indel" or die $!;
<IN>;
while(<IN>){
    chomp;
    my @temp = split /\s+/,$_;
    my $pos = (split /_/,$temp[0])[0];
    if(($temp[1] =~ /[ACGTacgt]/) && ($temp[2] =~ /[ACGTacgt]/)){
        print POS "$pos\n";
    }
    if(($temp[1] =~ /[ACGTacgt]/) && ($temp[2] =~ /\-/) || ($temp[1] =~ /\-/) && ($temp[2] =~ /[ACGTacgt]/)){
        if($pos > 265 && $pos < 29675 ){ #five_prime_UTR  1 265   #three_prime_UTR 29675   29903
            print IDL "$_\n";  
        }
    }
    $hash{"$temp[1]$temp[2]"} += 1;
}
close IN;
close POS;
close IDL;

open LOG,">$ARGV[1].log" or die $!;  
foreach my $key (sort keys %hash){
    print LOG "$key\t$hash{$key}\n";
}
close LOG;

