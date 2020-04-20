#! /usr/bin/perl -w
use strict;
#------------------------------------------------------------
#Summarize annotated vcf of snpEff
#------------------------------------------------------------

die "Usage: perl <eff_vcf_xls.np> <outfile>\n" unless @ARGV==2;

my @gene = qw/orf1a orf1b  S   ORF3a   E   M   ORF6    ORF7a  ORF7b  ORF8    N   ORF10/;

my %amino_acid = (
    Ala=>"A",
    Arg=>"R",
    Asn=>"N",
    Asp=>"D",
    Cys=>"C",
    Glu=>"E",
    Gln=>"Q",
    Gly=>"G",
    His=>"H",
    Ile=>"I",
    Leu=>"L",
    Lys=>"K",
    Met=>"M",
    Phe=>"F",
    Pro=>"P",
    Ser=>"S",
    Thr=>"T",
    Trp=>"W",
    Tyr=>"Y",
    Val=>"V",
);

open IN,"$ARGV[0]" or die $!;
open OUT,">$ARGV[1]" or die $!;
print OUT "Sample\tSNP_ns\tNT_muts\t",join("\t",@gene),"\n";
while(<IN>){
    chomp;
    next if (/^#/);
    my($name,$file) = (split /\s+/,$_);
    #my %nt;
    my $nt = "";
    my %aa;
    my ($SNP,$n,$s) = (0,0,0);
    open FILE,"$file" or die $!;
    <FILE>;
    while(<FILE>){
        chomp;
        my @temp = split /\s+/,$_;
        $SNP += 1;
        if($temp[4] =~ /synonymous_variant/){
            $s += 1;
        }elsif($temp[4] =~ /missense_variant/){
            $n += 1;
        }
        $nt .= "$temp[1]$temp[0]$temp[2],";
        $temp[9] =~ s/p\.//;
        $temp[9] =~ /(\D+)(\d+)(\D+)/;
        my ($ori,$pos,$tar) = ($1,$2,$3);
        if(!defined $ori){($ori,$pos,$tar) = ("NA","NA","NA");}
        next if ($ori eq $tar);
        $ori = (exists $amino_acid{$ori}) ? $amino_acid{$ori} : $ori;
        $tar = (exists $amino_acid{$tar}) ? $amino_acid{$tar} : $tar;
        $aa{$temp[6]} .= "$ori$pos$tar,";
    }
    close FILE;
    if($nt eq ""){$nt = "NA";}
    print OUT "$name\t$SNP\($n\/$s\)\t$nt";
    if(%aa){
        foreach my $gene (@gene) {
            if (exists $aa{$gene}){
                print OUT "\t$aa{$gene}";
                delete $aa{$gene};
            }else{
                print OUT "\tNA";
            }
        }
    }else{
        print OUT "\tNA" x 12;
    }
    if(%aa){
        print OUT "\t#",%aa,"\n";
    }else{
        print OUT "\n";
    }
}
close OUT;
close IN;
