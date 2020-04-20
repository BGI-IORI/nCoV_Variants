# nCoV_Variants

## Introduction
nCoV_Variants is a pipeline for variants calling of HCoV-19 virus. The pipeline could  call SNVs (Single-Nucleotide Variants) from HCoV-19 sequencing data.

![Image](https://github.com/BGI-IORI/nCoV_Variants/blob/master/nCoV_Variants.png)

## Requirements:
perl: v5.22.0  
python: v2.7.16   
java: v1.8.0  

Software for SNV calling:  
* Freebayes v1.3.1 (https://github.com/ekg/freebayes)  

Softwares for SNV filtering:  
* Snippy v3.2 (https://github.com/tseemann/snippy)

Software for SNV annotation:
* SnpEff v4.3 (http://snpeff.sourceforge.net/)  

Software for multiple sequence alignment
* mafft v7.222 (https://mafft.cbrc.jp/alignment/software/)  

Software for base calls
* pysamstats v1.1.2 (https://github.com/alimanfoo/pysamstats)

## Installation
```
git clone 
```
Notes: The above dependent software needs to be installed separately according to their instructions. After installing, the users should edit the input.config file, and change the software path to your own path.


## Usage
### 1.Prepare the input.config file
Edit the input.config file, and change the database path to your own path. We provide the HCoV-19 reference  (EPI_ISL_402119)  in the folder “database”, and you can also change it to your own reference.

### 2.Run the pipeline.
```
perl nCoV_Variants.pl -s input.config -b bam.txt -c CNS.txt -a ASS.txt -o ./outpath/
cd ./outpath/shellall/
sh allDependent.sh
#Notes: bam.txt file is required. CNS.txt and ASS.txt is optional
#bam.txt includes two columns, name and bam path (the bam should be sorted)
#CNS.txt includes two columns, name and consensus.fa (the orientation should be same with reference) 
#ASS.txt includes two columns, name and assembly.fa (the orientation should be same with reference, and it’s better to cut off the UTR region)
```

## Output
1.SNV from fFreebayes
```
./outpath/all_freebayesVar.xls
```
2.SNV summary from freebayes/consensus/assembly
```
./outpath/all_Var_Summary.xls
#The results will be generated only when CNS.txt or ASS.txt was as input.
```

