#! /usr/bin/perl -w
use strict;
#-------------------------------------------------
#Get diff pos from Multiple sequence alignment
#-------------------------------------------------

die "Usage : perl <MSA> <Ref_ID> <OUTmatrix>\n" unless @ARGV==3;

my $Refid = $ARGV[1] or die $!;
#######################################Read Postions and bp
###Read Ref position
my %Refpos;
$/ = "\>";
open IN,"$ARGV[0]" or die $!;
<IN>;
while(<IN>){
    chomp;
    my($name,$seq) = (split /\n/,$_,2)[0,1];
    $name = (split /\s+/,$name)[0];
    next unless ($name eq $Refid);
    $seq =~ s/\n//g;
    my @seq = split //,$seq;
    for(my $i=0;$i<=$#seq;$i++){
        $Refpos{$i} = $seq[$i];
    }
}
close IN;

###Read All seq position
$/ = "\>";
my %Allpos;
my %overlap;
my @allid;
open INN,"$ARGV[0]" or die $!;
<INN>;
while(<INN>){
    chomp;
    my($name,$seq) = (split /\n/,$_,2)[0,1];
    next if ($name eq $Refid);
    push @allid,$name; 
    $name = (split /\s+/,$name)[0];
    $seq =~ s/\n//g;
    my @seq = split //,$seq;
    for(my $i=0;$i<=$#seq;$i++){
        $Allpos{$name}{$i} = $seq[$i];
        if($Refpos{$i} ne $Allpos{$name}{$i}){
            $overlap{$i} += 1;
        }
    }
}
close INN;
##########################################output diff matrix
$/ = "\n";
open OUT,">$ARGV[2]" or die $!;
print OUT "Pos\t$Refid\t";
print OUT join "\t",@allid;
print OUT "\n";
foreach my $i (sort {$a <=> $b} keys %overlap){
    print OUT $i+1,"_",$overlap{$i},"\t",$Refpos{$i};
    foreach my $name (@allid){
        print OUT "\t$Allpos{$name}{$i}";  
    }
    print OUT "\n";
}
close OUT;
