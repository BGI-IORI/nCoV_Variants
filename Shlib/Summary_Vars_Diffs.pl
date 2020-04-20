#!/usr/bin/perl -w
use strict;
#---------------------------------------------------------------
#Summarize SNV from freebayes, consensus and assembly
#---------------------------------------------------------------

die "Usage : <raw.eff.vcf.txt> <flt.eff.vcf.txt> <CNS.xls> <ASS.xls> <bpInfo> <name> <SumVars.xls>\n" unless @ARGV==7;

my $sample = $ARGV[5];
my %allPos;

my %anno;
open RVCF,"$ARGV[0]" or die;
<RVCF>;
while(<RVCF>){
    chomp;
    my @temp=split(/\s+/,$_);
    my $pos=$temp[0];
    $anno{$pos}= "$temp[8]\t$temp[9]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]";
}
close RVCF;

my %fltPos;
open FVCF,"$ARGV[1]" or die;
<FVCF>;
while(<FVCF>){
    chomp;
    my @temp=split(/\s+/,$_);
    my $pos=$temp[0];
    $fltPos{$pos}{'ref'}= "$temp[1]";
    $fltPos{$pos}{'snp'}= "$temp[2]";
    $allPos{$pos} += 1;
}
close FVCF;

my %CNSPos;
open CNS,"$ARGV[2]" or die;
<CNS>;
while(<CNS>){
    chomp;   
    my @temp=split(/\s+/,$_);
    my $pos=(split /_/,$temp[0])[0];
    next if($pos>=29675 || $pos<=265);
    next if($temp[2] eq "-");
    next if($temp[2] =~ /[Nn]/);
    $CNSPos{$pos}{'ref'}= "$temp[1]";
    $CNSPos{$pos}{'snp'}= "$temp[2]";
    $allPos{$pos} += 1;
}
close CNS;

my %assPos;
if(!-e $ARGV[3]){`touch $ARGV[3]`;}
open ASS,"$ARGV[3]" or die;
<ASS>;
while(<ASS>){
    chomp;   
    my @temp=split(/\s+/,$_);
    my $pos=(split /_/,$temp[0])[0];
    next if($pos>=29675 && $pos<=265);
    next if($temp[2] eq "-");
    next if($temp[2] =~ /[Nn]/);
    $assPos{$pos}{'ref'}= "$temp[1]";
    $assPos{$pos}{'snp'}= "$temp[2]";
    $allPos{$pos} += 1;
}
close ASS;

open IN,"$ARGV[4]" or die;#INFO
open OUT,">$ARGV[6]" or die;#result
my $header = <IN>;
chomp($header);
$header =~ s/$//;
print OUT "Sample\tpos\tref\tfreebayes\tCNS\tASS\treads_pp\tdeletions_pp\tinsertions_pp\tA_pp\tC_pp\tT_pp\tG_pp\tN_pp\t";
print OUT "NT_muts\tAA_muts\tQUAL\tSNP_type\tIMPACT\tGENE\n";
while(<IN>){
    chomp;
    $_ =~ s/$//;
    my @file=split("\t",$_);
    next unless ($file[6] >=10); #depth >=10X
    if(exists $allPos{$file[1]}){
        print OUT "$sample\t$file[1]\t";
        $anno{$file[1]} = (defined $anno{$file[1]}) ? $anno{$file[1]} : "NA";
        $fltPos{$file[1]}{'snp'} = (defined $fltPos{$file[1]}{'snp'}) ? $fltPos{$file[1]}{'snp'} : "NA";
        $CNSPos{$file[1]}{'snp'} = (defined $CNSPos{$file[1]}{'snp'}) ? $CNSPos{$file[1]}{'snp'} : "NA";
        $assPos{$file[1]}{'snp'} = (defined $assPos{$file[1]}{'snp'}) ? $assPos{$file[1]}{'snp'} : "NA";
        print OUT "$file[2]\t";
        print OUT "$fltPos{$file[1]}{'snp'}\t$CNSPos{$file[1]}{'snp'}\t$assPos{$file[1]}{'snp'}\t";
        print OUT "$file[6]\t$file[24]\t$file[30]\t$file[36]\t$file[42]\t$file[48]\t$file[54]\t$file[60]\t";
        print OUT "$anno{$file[1]}\n";
    }
}
close IN;
close OUT;
