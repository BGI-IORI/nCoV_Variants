#! /usr/bin/perl -w
use Getopt::Long; 
use File::Spec;
use FindBin qw($Bin);

#---------------------------------------------------
#Version 1
#NOTE: Pipeline for analyzing SNV of nCoV genome
#---------------------------------------------------

my $usage=" NOTE: Pipeline for analyzing SNV of nCoV
    Usage:
    Options:
        -s <software>  input software input.config
        -b <bam.txt>   input bam.txt contain 2 col (name sample.bam)
        -c <CNS.txt>   input CNS.txt contain 2 col (name consensus.fasta)
        -a <ASS.txt>   input ASS.txt contain 2 col (name assembly.fasta)
        -o <outpath>   output path [./result]
        -h|?Help!
    Example:perl $0 -i data.txt -c input.config
";
my ($cfg,$bam,$cns,$ass,$outpath,$help);
GetOptions(
    '-s=s' => \$cfg,
    '-b=s' => \$bam,
    '-c=s' => \$cns,
    '-a=s' => \$ass,
    '-o=s' => \$outpath,
    'h|?'  => \$help,
);
if($help or !$bam){die "$usage\n";}
$outpath ||= "./result";
mkdir("$outpath") unless (-d $outpath);
$outpath = File::Spec->rel2abs($outpath);
my $shellall = "$outpath/shellall";
mkdir($shellall) unless (-d $shellall);
my %config = &readConf($cfg);

###Shlib
my $diff = "$Bin/Shlib/Show_diff_pos.pl";
my $diff_stat = "$Bin/Shlib/Show_diff_pos_stat.pl";
my $Parse_SnpEff = "$Bin/Shlib/Parse_SnpEff_bpInfo.pl";
my $Summary = "$Bin/Shlib/Summary_eff_vcf_xls.pl";
my $SumVars = "$Bin/Shlib/Summary_Vars_Diffs.pl";

my %CNS = &LIST($cns);
my %ASS = &LIST($ass);

open IN,"$bam" or die $!;
open RDS,">$outpath/shellall/A1_Var_Reads.sh" or die $!;
open CNS,">$outpath/shellall/A2_Var_CNS.sh" or die $!;
open ASS,">$outpath/shellall/A3_Var_ASS.sh" or die $!;
open SV,">$outpath/shellall/A4_SumVars.sh" or die $!;
open SUM,">$outpath/all_freebayesVar.np" or die $!;
while(<IN>){
    chomp;
    my($name,$bampath) = (split /\s+/,$_)[0,1];
    mkdir ("$outpath/$name") unless (-d "$outpath/$name");
    ###Variation: freebayes calling and SnpEff annotation
    open ORDS,">$outpath/$name/Var_Reads_$name.sh" or die $!; 
    print ORDS "$config{'freebayes'} $config{'freebayes_Parameters'} -f $config{'Ref'} $bampath > $outpath/$name/$name.raw.vcf 2>$outpath/$name/$name.raw.vcf.log && \\\n";
    print ORDS "$config{'snippy_Vcf_filter'} $config{'snippy_Parameters'} $outpath/$name/$name.raw.vcf > $outpath/$name/$name.flt.vcf 2>$outpath/$name/$name.flt.vcf.log && \\\n"; 
    print ORDS "$config{'java'} -Xmx14g -jar $config{'snpEff'} $config{'snpEff_Parameters'} $outpath/$name/$name.raw.vcf -c $config{'snpEff_cfg'} >$outpath/$name/$name.raw.eff.vcf -stats $outpath/$name/$name.raw.eff.html 2>$outpath/$name/$name.raw.eff.vcf.log && \\\n";
    print ORDS "$config{'java'} -Xmx14g -jar $config{'snpEff'} $config{'snpEff_Parameters'} $outpath/$name/$name.flt.vcf -c $config{'snpEff_cfg'} >$outpath/$name/$name.flt.eff.vcf -stats $outpath/$name/$name.flt.eff.html 2>$outpath/$name/$name.flt.eff.vcf.log && \\\n";
    print ORDS "$config{'pysamstats'} $config{'pysamstats_Parameters'} --fasta $config{'Ref'} $bampath > $outpath/$name/$name.bam.bpInfo && \\\n";
    print ORDS "$config{'perl'} $Parse_SnpEff $outpath/$name/$name.raw.eff.vcf $outpath/$name/$name.bam.bpInfo $outpath/$name/$name.raw.eff.vcf.xls && \\\n";
    print ORDS "$config{'perl'} $Parse_SnpEff $outpath/$name/$name.flt.eff.vcf $outpath/$name/$name.bam.bpInfo $outpath/$name/$name.flt.eff.vcf.xls && \\\n";
    print ORDS "echo ==========end at : `date` ========== && \\\necho Still_waters_run_deep 1>&2 && \\\n";  
    print ORDS "echo Still_waters_run_deep > $outpath/$name/Var_Reads_$name.sh.sign\n";
    close ORDS;
    print RDS "sh $outpath/$name/Var_Reads_$name.sh\n";
    print SUM "$name $outpath/$name/$name.flt.eff.vcf.xls\n"; 

    ###Variation: consensus-vs-reference
    if(exists $CNS{$name}){
        open OCNS,">$outpath/$name/Var_CNS_$name.sh" or die $!;
        print OCNS "cat $config{'Ref'} $CNS{$name} > $outpath/$name/VS_REF_CNS.$name.fa && \\\n";
        print OCNS "$config{'mafft'} $config{'mafft_Parameters'} $outpath/$name/VS_REF_CNS.$name.fa >$outpath/$name/VS_REF_CNS.$name.ms.fa 2>$outpath/$name/VS_REF_CNS.$name.ms.fa.log && \\\n";
        print OCNS "$config{'perl'} $diff $outpath/$name/VS_REF_CNS.$name.ms.fa '$config{'RefID'}' $outpath/$name/VS_REF_CNS.$name.ms.fa.xls && \\\n";
        print OCNS "$config{'perl'} $diff_stat $outpath/$name/VS_REF_CNS.$name.ms.fa.xls $outpath/$name/VS_REF_CNS.$name.ms.Diff && \\\n";
        print OCNS "echo ==========end at : `date` ========== && \\\necho Still_waters_run_deep 1>&2 && \\\n";
        print OCNS "echo Still_waters_run_deep > $outpath/$name/Var_CNS_$name.sh.sign\n";
        close OCNS;
        print CNS "sh $outpath/$name/Var_CNS_$name.sh\n";
    }
    
    ###Variation: assembly-vs-reference
    if(exists $ASS{$name}){
        open OASS,">$outpath/$name/Var_ASS_$name.sh" or die $!;
        print OASS "cat $config{'Ref'} $ASS{$name} > $outpath/$name/VS_REF_ASS.$name.fa && \\\n";
        print OASS "$config{'mafft'} $config{'mafft_Parameters'} $outpath/$name/VS_REF_ASS.$name.fa >$outpath/$name/VS_REF_ASS.$name.ms.fa 2>$outpath/$name/VS_REF_ASS.$name.ms.fa.log && \\\n";
        print OASS "$config{'perl'} $diff $outpath/$name/VS_REF_ASS.$name.ms.fa '$config{'RefID'}' $outpath/$name/VS_REF_ASS.$name.ms.fa.xls && \\\n";
        print OASS "$config{'perl'} $diff_stat $outpath/$name/VS_REF_ASS.$name.ms.fa.xls $outpath/$name/VS_REF_ASS.$name.ms.Diff && \\\n";
        print OASS "echo ==========end at : `date` ========== && \\\necho Still_waters_run_deep 1>&2 && \\\n";
        print OASS "echo Still_waters_run_deep > $outpath/$name/Var_ASS_$name.sh.sign\n";
        close OASS;
        print ASS "sh $outpath/$name/Var_ASS_$name.sh\n";
    }

    ###SumVars: freebayes-CNS-ass
    open OSV,">$outpath/$name/SumVars_$name.sh" or die $!;
    print OSV "perl $SumVars $outpath/$name/$name.raw.eff.vcf.xls $outpath/$name/$name.flt.eff.vcf.xls $outpath/$name/VS_REF_CNS.$name.ms.fa.xls $outpath/$name/VS_REF_ASS.$name.ms.fa.xls $outpath/$name/$name.bam.bpInfo $name $outpath/$name/SumVars.$name.xls\n";
    close OSV;
    print SV "sh $outpath/$name/SumVars_$name.sh\n";
}
close RDS;
close ASS;
close CNS;
close SV;
close IN;

open ST,">$outpath/all_stat.sh" or die $!;
print ST "perl $Summary $outpath/all_freebayesVar.np $outpath/all_freebayesVar.xls\n" or die $!;
print ST "cat $outpath/*/SumVars.*.xls > $outpath/all_Var_Summary.xls\n";
close ST;

###################all shells
`cat $outpath/shellall/A1_Var_Reads.sh $outpath/shellall/A2_Var_CNS.sh $outpath/shellall/A3_Var_ASS.sh $outpath/shellall/A4_SumVars.sh $outpath/all_stat.sh > $outpath/shellall/allDependent.sh`;

sub LIST{
    my $file = shift;
    my %hash;
    open LST,"$file" or die $!;
    while(<LST>){
        chomp;
        my($name,$path) = split /\s+/,$_;
        $hash{$name} = $path;
    }
    close LST;
    return %hash;
}

sub readConf{
    my $confFile = shift @_;
    my %hash;
    open IN, $confFile or die "Cannot open file $confFile:$!\n";
    while (<IN>){
        chomp;
        next if(/^\s*$/ || /^\s*\#/);
        $_ =~ s/^\s*//;
        $_ =~ s/#(.)*//;
        $_ =~ s/\s*$//;
        if (/^(\w+)\s*=\s*(.*)$/xms){
            next if ($2 =~ /^\s*$/);
            my $key = $1;
            my $value = $2;
            $value =~ s/\s*$//;
            $hash{$key} = $value;
        }
    }
    return %hash;
}
