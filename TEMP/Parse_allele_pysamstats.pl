#!/usr/bin/perl -w
use strict;
#--------------------------------------------------------
#bpInfo.txt include two cols: name output_of_pysamstats
#---------------------------------------------------------
die "Usage : perl <bpInfo.txt> <outprefix>\n" unless @ARGV==2;

my $heteRatio = 0.05;
my $DEP = 10;
open IN,"$ARGV[0]" or die;
open OUT,">$ARGV[1].hete.xls" or die;#result
open SUM,">$ARGV[1].hete.Summary.xls" or die;#count
open FIG,">$ARGV[1].hete.List.xls" or die $!;
print OUT "Sample\tpos\tmutates\treads_pp\tdeletions_pp\tinsertions_pp\tA_pp\tC_pp\tT_pp\tG_pp\n";
print SUM "Sample\tHete_count\tCoverbp\tHete_rate\n";

my @stran;
my @sort_stran;
my @stran_minor;
my @sort_stran_minor;

while(<IN>){
    chomp;
    my($name,$file) = (split /\s+/,$_)[0,1];
    my ($all,$n) = (0,0);
    open FILE,"$file" or die $!;
    <FILE>;
    while(<FILE>){
        chomp;
        my @file=split("\t",$_);  
        my @number=($file[24],$file[30],$file[36],$file[42],$file[48],$file[54]);
        if($file[2] =~ m/A|T|C|G/){
            $all=$all+1;
        }next if ($file[1]<=265 || $file[1]>=29675);
        	my $depth=$file[6];
        	if($depth >=$DEP ){
            	my @sorted_numbers = sort {$b <=> $a} @number;
            	my $maxindex= $#number;
            	my $ratio=$depth*$heteRatio;
	   	my %count = ($file[24]=>"del",$file[30]=>"ist",$file[36]=>"A",$file[42]=>"C",$file[48]=>"T",$file[54]=>"G");
 	   	my $Major = $count{$sorted_numbers[0]};#$count{(sort ({$count{$b} <=> $count{$a}}keys %count))[0]};
		my $Minor = $count{$sorted_numbers[1]}; 
		@stran=""; @sort_stran="";
		if($Major eq "del"){
			@stran = ($file[25],$file[26]);
		  	@sort_stran = sort {$b <=> $a} @stran;	
			} elsif($Major eq "ist"){
		  	@stran = ($file[31],$file[32]);
                   	@sort_stran = sort {$b <=> $a} @stran;
			}elsif($Major eq "A"){
                  	@stran = ($file[37],$file[38]);
                 	@sort_stran = sort {$b <=> $a} @stran;
             		}elsif($Major eq "C"){
                  	@stran = ($file[43],$file[44]);
                 	@sort_stran = sort {$b <=> $a } @stran;
			}elsif($Major eq "T"){
                  	@stran = ($file[49],$file[50]);
                 	@sort_stran = sort {$b <=> $a} @stran;
			}else{
                  	@stran = ($file[55],$file[56]);
                 	@sort_stran = sort {$b <=> $a} @stran;
			}	

	 	@stran_minor=""; @sort_stran_minor="";
                if($Minor eq "del"){
                  @stran_minor = ($file[25],$file[26]);
                  @sort_stran_minor = sort {$b <=> $a} @stran_minor;
                  } elsif($Minor eq "ist"){
                  @stran_minor = ($file[31],$file[32]);
                  @sort_stran_minor = sort {$b <=> $a} @stran_minor;
                }elsif($Minor eq "A"){
                  @stran_minor = ($file[37],$file[38]);
                  @sort_stran_minor = sort {$b <=> $a} @stran_minor;
             	}elsif($Minor eq "C"){
                  @stran_minor = ($file[43],$file[44]);
                  @sort_stran_minor = sort {$b <=> $a } @stran_minor;
                }elsif($Minor eq "T"){
                  @stran_minor = ($file[49],$file[50]);
                  @sort_stran_minor = sort {$b <=> $a} @stran_minor;
                }else{
                  @stran_minor = ($file[55],$file[56]);
                  @sort_stran_minor = sort {$b <=> $a} @stran_minor;
		}


	if($sorted_numbers[1] >= $ratio  && $sort_stran_minor[1]*10>$sort_stran_minor[0] && $sort_stran_minor[1]>4){

		my $ref_call = $stran[1]/($stran[1]+$stran[0]);
		my $var_call = $stran_minor[1]/($stran_minor[1]+$stran_minor[0]);
		my $ratio2 = $ref_call/$var_call;
		if ($ratio2 <5 && $ratio2 >0.2 ){
	                $n=$n+1;
               		print OUT "$name\t$file[1]\t$Minor\t$depth\t";
                	print OUT "$file[24]\t$file[30]\t$file[36]\t$file[42]\t$file[48]\t$file[54]\n";
                	print FIG "$name\t",$sorted_numbers[1]/$depth*100,"\n";
            	}
        }
    }
}
    close FILE;
    print SUM "$name\t$n\t$all\t",$n/$all,"\n";
}
close IN;
close OUT;
close SUM;
close FIG;
