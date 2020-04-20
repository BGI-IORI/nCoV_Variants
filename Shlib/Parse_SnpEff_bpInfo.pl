#!/usr/bin/perl -w
use strict;
#------------------------------------------------------------------------
#Input file are snpEff annotated vcf and output of bpInfo from pysamstats
#------------------------------------------------------------------------

die "Usage : perl <SnpEff.vcf> <bpInfo> <outFile>\n" unless @ARGV==3;

###Variation SnpEff vcf file
my %hash;
open IN,"$ARGV[0]" or die $!;
while(<IN>){
	next if (/^#/);
    chomp;
    my @temp=split(/ANN\=/,$_);
    my @eff = split /\|/,$temp[1];
    my @allsnp=split(/\s+/,$temp[0]);
    my $pos=$allsnp[1];
    my $ref=$allsnp[3];
    my $snp=$allsnp[4];
    my $qual = $allsnp[5];
    $eff[1] = ($eff[1] ne "") ? $eff[1] : "NA";
    $eff[2] = ($eff[2] ne "") ? $eff[2] : "NA";
    $eff[3] = ($eff[3] ne "") ? $eff[3] : "NA";
    $eff[8] = ($eff[8] ne "") ? $eff[8] : "NA";
    $eff[9] = ($eff[9] ne "") ? $eff[9] : "NA";
    $eff[10] = ($eff[10] ne "") ? $eff[10] : "NA";
    $hash{$pos} = "$ref\t$snp\t$qual\t$eff[1]\t$eff[2]\t$eff[3]\t$eff[8]\t$eff[9]\t$eff[10]";
}
close IN;

###pysamstats results
open INN,"$ARGV[1]" or die $!;
open OUT,">$ARGV[2]" or die $!;
print OUT "POSITION\tREF\tSNP\tQUAL\tSNP_type\tIMPACT\tGENE\tGENOTYPE\tNT_muts\tAA_muts\tDEPTH\tDEL\tINS\tA\tC\tT\tG\n";
while(<INN>){
    chomp;
    my @file = split("\t",$_);  
    if(exists $hash{$file[1]}){
        print OUT "$file[1]\t$hash{$file[1]}\t";
        print OUT "$file[6]\t$file[24]\t$file[30]\t$file[36]\t$file[42]\t$file[48]\t$file[54]\n";
    }
}
close OUT;
close INN;
