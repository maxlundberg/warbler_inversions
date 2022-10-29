# Demographic and selection analyses


## Demographic modelling with gIMble

Use `freebayes` to call variants in the four high coverage resequenced samples

```
freebayes -f ww_southern_hifi_bionano.filt.new.fasta -b 1A05.sorted.nodup.bam -b UK06.sorted.nodup.bam -b 1L19.sorted.dedup.bam -b 1M13.sorted.nodup.bam -b DW83.sorted.nodup.bam -v pacbio.vcf
gzip -c pacbio.vcf > pacbio.vcf.gz
```


### Run gIMble

Preprocess the data. The file `bam_symlinks` provides a list of the bam files used for variant calling

```
gIMble preprocess -f ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -v pacbio.vcf.gz -b bam_symlinks -g 2 -m 8 -M 0.75 -t 8
```


Next, set up zarr stores

```
gIMble setup -g chromosome_1.gimble.genomefile -b chromosome_1.gimble.bed -s chromosome_1.gimble.samples -v gimble.vcf.gz -o chromosome_1.64
gIMble setup -g chromosome_3.gimble.genomefile -b chromosome_3.gimble.bed -s chromosome_3.gimble.samples -v gimble.vcf.gz -o chromosome_3.64
gIMble setup -g chromosome_5.gimble.genomefile -b chromosome_5.gimble.bed -s chromosome_5.gimble.samples -v gimble.vcf.gz -o chromosome_5.64
```


Get blocks

```
gIMble blocks -z chromosome_1.64.z -l 64 --force
gIMble blocks -z chromosome_3.64.z -l 64 --force
gIMble blocks -z chromosome_5.64.z -l 64 --force
```


Set up models

```
gIMble model -p 2 -s A,B -n 1,1 -j A,B
gIMble model -p 2 -s A,B -n 1,1 -j A,B -m 'A>B'
gIMble model -p 2 -s A,B -n 1,1 -j A,B -m 'B>A'
```


For each divergent region, optimise each model

```
gIMble optimize -z chromosome_1.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.ini -m gimble.model.s_A_B.n_1_1.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.ini -m gimble.model.s_A_B.n_1_1.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_5.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.ini -m gimble.model.s_A_B.n_1_1.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_1.64.z -c gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_5.64.z -c gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_1.64.z -c gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.tsv -b -t 16,1 -i 1000
gIMble optimize -z chromosome_5.64.z -c gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.tsv -b -t 16,1 -i 1000
```


Simulate SI data for each divergent region

```
gIMble simulate -z chromosome_1.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.sims_1_SI.ini -b 72432 -r 100 -t 8 -l optimised_SI
gIMble simulate -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.sims_3_SI.ini -b 83537 -r 100 -t 8 -l optimised_SI
gIMble simulate -z chromosome_5.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.sims_5_SI.ini -b 28351 -r 100 -t 8 -l optimised_SI
```


For each simulated SI dataset, optimise under an SI and the best fitting (IM) model

```
gIMble optimize -z chromosome_1.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.ini -m gimble.model.s_A_B.n_1_1.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_SI
gIMble optimize -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.ini -m gimble.model.s_A_B.n_1_1.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_SI
gIMble optimize -z chromosome_5.64.z -c gimble.model.s_A_B.n_1_1.J_A_B.ini -m gimble.model.s_A_B.n_1_1.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_SI
gIMble optimize -z chromosome_1.64.z -c gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_SI
gIMble optimize -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_SI
gIMble optimize -z chromosome_5.64.z -c gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_A_B.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_SI
```


For the chromosome 3 region, simulate and optimise IM2 data to obtain 95% CIs

```
gIMble simulate -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.sims_3_SI.ini -b 83537 -r 100 -t 8 -l optimised_IM
gIMble optimize -z chromosome_3.64.z -c gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.ini -m gimble.model.s_A_B.n_1_1.M_B_A.J_A_B.tsv -s -t 7,7 -i 1000 -l optimised_IM
```


## Calculate divergence and RND in the divergent regions

Use `freebayes` to call variants for high coverage willow warblers alongside the dusky warbler in each of the three divergent regions 

```
freebayes -f ww_southern_hifi_bionano.filt.new.fasta --region Scaffold19:58159-11678265 --bam 1A05.sorted.nodup.bam --bam 1L19.sorted.dedup.bam --bam 1M13.sorted.nodup.bam --bam DW83.sorted.nodup.bam --bam UK06.sorted.nodup.bam | gzip -c > freebayes_DW83_southern_hifi_bionano.raw.chromosome_1.vcf.gz
freebayes -f ww_southern_hifi_bionano.filt.new.fasta --region Scaffold61:56078010-69204987 --bam 1A05.sorted.nodup.bam --bam 1L19.sorted.dedup.bam --bam 1M13.sorted.nodup.bam --bam DW83.sorted.nodup.bam --bam UK06.sorted.nodup.bam | gzip -c > freebayes_DW83_southern_hifi_bionano.raw.chromosome_3.vcf.gz
freebayes -f ww_southern_hifi_bionano.filt.new.fasta --region Scaffold0:327560-4361789 --bam 1A05.sorted.nodup.bam --bam 1L19.sorted.dedup.bam --bam 1M13.sorted.nodup.bam --bam DW83.sorted.nodup.bam --bam UK06.sorted.nodup.bam | gzip -c > freebayes_DW83_southern_hifi_bionano.raw.chromosome_5.vcf.gz
```

Use gIMble preprocess to filter the vcfs

```
gIMble preprocess -f data_2022/ww_southern_hifi_bionano.filt.new.fasta -v data_2022/freebayes_DW83_southern_hifi_bionano.raw.chromosome_1.vcf.gz -b preprocess_DW83_2022/bams -g 2 -q 1 -m 8 -M 3 -t 20 -o preprocess_DW83_2022/phylloscopus.DW83.chromosome_1 -k
gIMble preprocess -f data_2022/ww_southern_hifi_bionano.filt.new.fasta -v data_2022/freebayes_DW83_southern_hifi_bionano.raw.chromosome_3.vcf.gz -b preprocess_DW83_2022/bams -g 2 -q 1 -m 8 -M 3 -t 20 -o preprocess_DW83_2022/phylloscopus.DW83.chromosome_3 -k
gIMble preprocess -f data_2022/ww_southern_hifi_bionano.filt.new.fasta -v data_2022/freebayes_DW83_southern_hifi_bionano.raw.chromosome_5.vcf.gz -b preprocess_DW83_2022/bams -g 2 -q 1 -m 8 -M 3 -t 20 -o preprocess_DW83_2022/phylloscopus.DW83.chromosome_5 -k
```

Filter bed files for each divergent region so that they contain sequences that are intergenic and non-repeat

```
chroms="1 3 5"
for chrom in $chroms
do
bedtools subtract -a phylloscopus.DW83.chromosome_$chrom.bed -b ww_southern_hifi_augustus_webapollo.combined.genes.bed > phylloscopus.DW83.chromosome_$chrom.intergenic.bed
bedtools subtract -a phylloscopus.DW83.chromosome_$chrom.intergenic.bed -b ww_southern_hifi_bionano.filt.new.repeats.bed > phylloscopus.DW83.chromosome_$chrom.intergenic.no_repeats.bed
done 
```

These files will also be trimmed to exclude any non-divergent sequence from the bed files and have the suffix "subset" added. From these files we will
use gIMble to setup stores, make blocks and finally calculate genetic summary statistics

```
chroms="1 3 5"
for chrom in $chroms
do
gIMble setup -g preprocess_DW83_2022/phylloscopus.DW83.chromosome_$chrom.genomefile -b preprocess_DW83_2022/phylloscopus.DW83.chromosome_$chrom.intergenic.no_repeats.subset.bed -s preprocess_DW83_2022/phylloscopus.DW83.chromosome_$chrom.samples.csv -v preprocess_DW83_2022/phylloscopus.DW83.chromosome_$chrom.vcf.gz -o chromosome_$chrom/phylloscopus.DW83.chromosome_$chrom
gIMble blocks -z chromosome_$chrom/phylloscopus.DW83.chromosome_$chrom.z -l 64 --force
gIMble info -z chromosome_$chrom/phylloscopus.DW83.chromosome_$chrom.z
done
```


## MSMC2


Split the vcf file containing variants called in high coverage willow warbler samples by samples and by scaffold. For each of these vcf files,
filter them to contain only heterozygous variants, split the callable sites file the same way and generate input for MSMC2.

```
#First get positions of intergenic sites not overlapping with repeats, where phylloscopus.bed is all the scaffold intervals
bedtools subtract -a phylloscopus.bed -b ww_southern_hifi_augustus_webapollo.combined.genes.bed > phylloscopus.intergenic.bed
bedtools subtract -a phylloscopus.bed -b ww_southern_hifi_bionano.filt.new.repeats.bed > phylloscopus.intergenic.no_repeats.bed

samples="1A05 1L19 1M13 UK06"
for sample in $samples
do
for i in {0..545}
do

bcftools filter -r Scaffold$i pacbio.vcf.gz | bcftools view -s $sample > phylloscopus.Scaffold$i.$sample.vcf

bcftools filter -i 'AC == 1 && AN == 2' phylloscopus.Scaffold$i.$sample.vcf > phylloscopus.Scaffold$i.$sample.clean.vcf
gzip phylloscopus.Scaffold$i.$sample.clean.vcf
rm -f phylloscopus.Scaffold$i.$sample.vcf

cat phylloscopus.intergenic.no_repeats.bed | grep -P "Scaffold$i\t" | grep $sample | bedtools merge > phylloscopus.Scaffold$i.$sample.bed

~/software/msmc2/msmc-tools/generate_multihetsep.py --mask phylloscopus.Scaffold$i.$sample.bed scaffold$i.$sample.clean.vcf.gz > phylloscopus.Scaffold$i.$sample.MSMC2.input
done
done                                                                
```

Remove Scaffold0, Scaffold19 and Scaffold61 input files (divergent regions) as well as scaffolds shorter than 500kb and empty files

```
falen ww_southern_hifi_bionano.filt.new.fasta | sort -k2 -nr | tail -436 | cut -f 1 > short_sequences.txt

cat short_sequences.txt | while read line
do rm -f phylloscopus.$line.*
done

find -empty | while read file; do rm -f $file; done
```

Run MSMC2

```
for sample in $samples
do
msmc2_linux64bit -i 30 -o $sample -r 1 -t 24 -I 0,1 ../2_input/phylloscopus.Scaffold*.$sample.MSMC2.input
done
```



## Get phased variants

Phase the filtered dataset with `beagle`

```
#There are some scaffolds with only one variant. If not removing these variants, the software will stop
echo "Scaffold155:13151" > exclude_markers.list
echo "Scaffold208:106810" >> exclude_markers.list
echo "Scaffold268:2233" >> exclude_markers.list
echo "Scaffold284:22680" >> exclude_markers.list
echo "Scaffold468:27408" >> exclude_markers.list
echo "Scaffold473:101" >> exclude_markers.list

java -Xmx240G -jar ~/beagle.22Jul22.46e.jar gt=freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz out=freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased excludemarkers=exclude_markers.list
```

Extract bi-allelic SNPs and create a file for pure northern and pure southern samples

```
vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.vcf.gz --recode --remove-indels --max-alleles 2 --stdout | bgzip -c > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.vcf.gz

vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.vcf.gz --keep ~/ww_samples_pure_southern.txt  --recode --stdout | bgzip -c > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.vcf.gz
vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.vcf.gz --keep ~/ww_samples_pure_northern.txt  --recode --stdout | bgzip -c > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.vcf.gz
```

Calculate nucleotide diversity and Tajima's D in 10 kb windows for pure southern and northern samples using `vcftools`

```
vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.vcf.gzfreebayes_reseq_southern_hifi_bionano.raw.vcf.gz --stdout > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.tajd_10kb.out
vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.vcf.gz --TajimaD 10000 --stdout > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.tajd_10kb.out
vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.vcf.gz --window-pi 10000 --stdout > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pi_10kb.out
vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.vcf.gz --window-pi 10000 --stdout > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pi_10kb.out
```

## Infer ancestral and derived variants 

Use `bcftools` and customized scripts to get genotypes from the aligned chiffchaff and dusky warbler reads for the phased bi-allelic willow warbler SNPs

```
zcat freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.vcf.gz | grep -v "#" | cut -f1-2 > positions.phased_snps.list

bcftools mpileup DW83.sorted.nodup.bam --fasta-ref ww_southern_hifi_bionano.filt.new.fasta -T positions.phased_snps.list --config illumina --annotate INFO/AD >  DW.bcftools.phased_pos.vcf
bcftools mpileup ww_southern_hifi_bionano.filt.fasta.chiffchaff.hifi_reads.sorted.bam --fasta-ref ww_southern_hifi_bionano.filt.new.fasta -T positions.phased_snps.list --config pacbio-ccs --annotate INFO/AD > CC.bcftools.phased_pos.vcf

cat DW.bcftools.phased_pos.vcf | extract_allele_counts_bcftools_mpileup.pl > DW.bcftools.phased_pos.ac.out
cat CC.bcftools.phased_pos.vcf | extract_allele_counts_bcftools_mpileup.pl > CC.bcftools.phased_pos.ac.out

#Calculate the mean coverage for each species across the genotypes
cat DW.bcftools.phased_pos.ac.out | perl -ne 'use List::Util qw(sum);$tot_sum=0;$n_sites;while(<STDIN>){my @input=split("\t",$_);$tot_sum+=sum(@input[2..$#input]);$n_sites++} printf("%.0f",$tot_sum/$n_sites); print "\n"' #33
cat CC.bcftools.phased_pos.ac.out | perl -ne 'use List::Util qw(sum);$tot_sum=0;$n_sites;while(<STDIN>){my @input=split("\t",$_);$tot_sum+=sum(@input[2..$#input]);$n_sites++} printf("%.0f",$tot_sum/$n_sites); print "\n"' #43

#Summarize the outgroup data
find_ancestral_allele_outgroup_species.2.pl --ref_species_vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.vcf.gz --outgroup1_file CC.bcftools.phased_pos.ac.out --outgroup2_file DW.bcftools.phased_pos.ac.out --outgroup1_cov 43 --outgroup2_cov 33 > outgroup_allele_info.out

#Change the vcf files 
polarize_alleles_vcf.pl --ref_species_vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.vcf.gz --outgroup_info outgroup_allele_info.out | bgzip -c > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pol.vcf.gz
tabix -p vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pol.vcf.gz
polarize_alleles_vcf.pl --ref_species_vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.vcf.gz --outgroup_info outgroup_allele_info.out | bgzip -c > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pol.vcf.gz
tabix -p vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pol.vcf.gz
```


## Calculate XP-nsl

First divide the northern and southern phased and polarized vcf files into scaffold-specific files 

```
mkdir vcf_files_pol_selscan

zcat freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pol.vcf.gz | grep -v "#" | cut -f1 > scaffolds.phased.vcf.list

cat scaffolds.phased.vcf.list | while read scaffold
do
bcftools view freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pol.vcf.gz $scaffold | bgzip -c > vcf_files_pol_selscan/$scaffold.northern.vcf.gz
bcftools view freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pol.vcf.gz $scaffold | bgzip -c > vcf_files_pol_selscan/$scaffold.southern.vcf.gz
done
```

Next, use `selscan 1.3.0` to compute XP-nsl between pure southern and northern willow warblers 

```
mkdir selscan_xpnsl_pol

cat scaffolds.phased.vcf.list | while read scaffold
do
selscan --vcf vcf_files_pol_selscan/$scaffold.southern.vcf.gz --vcf-ref vcf_files_pol_selscan/$scaffold.northern.vcf.gz --threads 20 --out selscan_xpnsl_pol/$scaffold --xpnsl
mv selscan_xpnsl_pol/$scaffold.xpnsl.out selscan_xpnsl_pol/$scaffold.xpnsl.out.old
cat selscan_xpnsl_pol/$scaffold.xpnsl.out.old | sed "s/^./$scaffold/" > selscan_xpnsl_pol/$scaffold.xpnsl.out
rm selscan_xpnsl_pol/$scaffold.xpnsl.out.old
done
```

Normalize the raw values in 10 kb windows using `norm` from the selscan package

```
norm --xpnsl --files selscan_xpnsl_pol/*.out --bp-win --winsize 10000

echo -e "scaffold\tstart\tend\tN_scores_win\tfrac_scores_gt_threshold\tfrac_scores_lt_threshold\tapprox_perc_gt_win\tapprox_perc_lt_win\tmax_score\tmin_score" > xpnsl.10kb.out

cat scaffolds.phased.vcf.list | while read scaffold
do
export scaffold
cat selscan_xpnsl_pol/$scaffold.xpnsl.out.norm.10kb.windows | perl -ne 'while(<STDIN>){print "$ENV{scaffold}\t$_"}' >> xpnsl.10kb.out
done
```


## Calculate Fay and Wu's H

Add the outgroup consensus genotype to the phased vcf files

```
./add_outgroup_geno_vcf.pl --ref_species_vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.vcf.gz --outgroup_info outgroup_allele_info.out | bgzip -c > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.with_anc.vcf.gz
tabix -p vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.with_anc.vcf.gz

./add_outgroup_geno_vcf.pl --ref_species_vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.vcf.gz --outgroup_info outgroup_allele_info.out | bgzip -c > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.with_anc.vcf.gz
tabix -p vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.with_anc.vcf.gz
```

Load the data into `PopGenome` in R. First, SNP intervals for each scaffold will be calculated to improve the data loading.   

```
zcat  freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.with_anc.vcf.gz | grep -v "#" | cut -f1-2 | sort -k1,1 -k2,2nr | sort -u -k1,1 | perl -lane 'print "$F[0]\t1\t$F[1]"' > positional_info.vcf_files_with_anc.txt

#In R
library(PopGenome)

scaffold_info=read.delim("positional_info.vcf_files_with_anc.txt",header=F)

for(i in c(1:nrow(scaffold_info))){
	scaffold=scaffold_info[i,1]
	start=scaffold_info[i,2]
	end=scaffold_info[i,3]

input=readVCF("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.with_anc.vcf.gz",frompos=start,topos=end,tid=scaffold,numcols=100)
input<-set.outgroup(input,c("OUTGROUP"),diploid =TRUE)

sliding_windows_data<-sliding.window.transform(input,width=10000,jump=10000,type=2)

test=try(sliding_windows_data<-neutrality.stats(sliding_windows_data),silent=T)

if(class(test)=="try-error"){next}
else{
scaff_data=as.data.frame(array(0,c(length(sliding_windows_data@region.names),4)))
scaff_data[,1]=rep(scaffold,nrow(scaff_data))
scaff_data[,2]=sliding_windows_data@region.names
scaff_data[,3]=sliding_windows_data@Tajima.D
scaff_data[,4]=sliding_windows_data@Fay.Wu.H
write.table(scaff_data,file="southern_with_anc_popgenome.out",append=T,quote=F,row.names=F,col.names=F,sep="\t")
}

#
input=readVCF("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.with_anc.vcf.gz",frompos=start,topos=end,tid=scaffold,numcols=100)
input<-set.outgroup(input,c("OUTGROUP"),diploid =TRUE)

sliding_windows_data<-sliding.window.transform(input,width=10000,jump=10000,type=2)

test=try(sliding_windows_data<-neutrality.stats(sliding_windows_data),silent=T)

if(class(test)=="try-error"){next}
else{
scaff_data=as.data.frame(array(0,c(length(sliding_windows_data@region.names),4)))
scaff_data[,1]=rep(scaffold,nrow(scaff_data))
scaff_data[,2]=sliding_windows_data@region.names
scaff_data[,3]=sliding_windows_data@Tajima.D
scaff_data[,4]=sliding_windows_data@Fay.Wu.H
write.table(scaff_data,file="northern_with_anc_popgenome.out",append=T,quote=F,row.names=F,col.names=F,sep="\t")
}
```

## Sweepfinder2

Use `Sweepfinder2` to detect signals of recent positive selection. To get suitable input data we will first extract derived genotype counts
using `vcftools` 

```
echo -e "position\tx\tn\tfolded" > header.sweepfinder.txt

vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pol.vcf.gz --counts2 --stdout | sed '1d' | cut -f1,2,6 > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pol.counts.out
vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pol.vcf.gz --counts2 --stdout | sed '1d' | cut -f1,2,6 > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pol.counts.out

(cat header.sweepfinder.txt && cat freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pol.counts.out | perl -lane 'if($F[2]>0){$,="\t";print $F[1],$F[2],20,1}') > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.sweepfinder.freq
(cat header.sweepfinder.txt && cat freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pol.counts.out | perl -lane 'if($F[2]>0){$,="\t";print $F[1],$F[2],14,1}') > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.sweepfinder.freq
```


Get a global allele frequency file from both the southern and the northern samples

```
SweepFinder2 -f freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.sweepfinder.freq Spectsouth
SweepFinder2 -f freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.sweepfinder.freq Spectnorth
```

Create input files for each scaffold containing the divergent regions 

```
(cat header.sweepfinder.txt && cat freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pol.counts.out | perl -lane 'if($F[0] eq "Scaffold0" && $F[2]>0){$,="\t";print $F[1],$F[2],20,1}') > Scaffold0.southern.sweepfinder.freq
(cat header.sweepfinder.txt && cat freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pol.counts.out | perl -lane 'if($F[0] eq "Scaffold19" && $F[2]>0){$,="\t";print $F[1],$F[2],20,1}') > Scaffold19.southern.sweepfinder.freq
(cat header.sweepfinder.txt && cat freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pol.counts.out | perl -lane 'if($F[0] eq "Scaffold61" && $F[2]>0){$,="\t";print $F[1],$F[2],20,1}') > Scaffold61.southern.sweepfinder.freq

(cat header.sweepfinder.txt && cat freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pol.counts.out | perl -lane 'if($F[0] eq "Scaffold0" && $F[2]>0){$,="\t";print $F[1],$F[2],14,1}') > Scaffold0.northern.sweepfinder.freq
(cat header.sweepfinder.txt && cat freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pol.counts.out | perl -lane 'if($F[0] eq "Scaffold19" && $F[2]>0){$,="\t";print $F[1],$F[2],14,1}') > Scaffold19.northern.sweepfinder.freq
(cat header.sweepfinder.txt && cat freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pol.counts.out | perl -lane 'if($F[0] eq "Scaffold61" && $F[2]>0){$,="\t";print $F[1],$F[2],14,1}') > Scaffold61.northern.sweepfinder.freq
```


Run `sweepfinder2` on each of the input files and specify 10000 between each grid point. To speed things up we will use `GNU parallel`.

```
echo -e "Scaffold0.southern Spectsouth\nScaffold0.northern Spectnorth\nScaffold61.southern Spectsouth\nScaffold61.northern Spectnorth\nScaffold19.southern Spectsouth\nScaffold19.northern Spectnorth" > info.txt

parallel -j 6 --col-sep " " '~/SF2/SweepFinder2 -lg 10000 {1}.sweepfinder.freq {2} {1}.sweepfinder.out' :::: info.txt
```


## Calculate LD

LD will be calculated in 10 kb windows using `vcftools` and data extraction using `bcftools`. For each window, we will obtain a mean value for r^2 and D'
among SNPs that have a minor allele frequency of at least 0.2 and are at least 1 kb from each other. As LD measurements are sensitive to sample size
differences, we will randomly remove three southern willow warblers (0G03, 0G04 and UK06) to get even sample sizes. 

```
#Southern samples
echo -e "0G03\n0G04\nUK06" > samples_southern_remove.txt

vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.vcf.gz --remove samples_southern_remove.txt --recode --recode-INFO-all --stdout | bgzip -c > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.7_samples.vcf.gz
tabix -p vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.7_samples.vcf.gz 

cat ww_southern_hifi_bionano.filt.new.10kb_windows.bed | while read chr start end
do
start2=$(echo "$start+1" | bc)
ld_data=$(bcftools view freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.7_samples.vcf.gz $chr:$start-$end | vcftools --vcf - --ld-window-bp-min 1000  --hap-r2 --maf 0.2 --stdout | sed '1d' | cut -f5,7 | grep -v "nan" | perl -e '$sum_r2=0;$sum_dprime=0;$count=0;while($input=<STDIN>){($r2,$dprime)=split("\t",$input);$sum_dprime+=abs($dprime);$sum_r2+=$r2;$count++}if($count>0){print $sum_r2/$count,"\t",$sum_dprime/$count,"\n"}else{print "0\t0\n"}')
echo -e "$chr\t$start\t$end\t$ld_data"
done > ww_southern_hifi_bionano.filt.southern.7_samples.ld.out

#Northern samples

cat ww_southern_hifi_bionano.filt.new.10kb_windows.bed | while read chr start end
do
start2=$(echo "$start+1" | bc)
ld_data=$(bcftools view freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.vcf.gz $chr:$start-$end | vcftools --vcf - --ld-window-bp-min 1000  --hap-r2 --maf 0.2 --stdout | sed '1d' | cut -f5,7 | grep -v "nan" | perl -e '$sum_r2=0;$sum_dprime=0;$count=0;while($input=<STDIN>){($r2,$dprime)=split("\t",$input);$sum_dprime+=abs($dprime);$sum_r2+=$r2;$count++}if($count>0){print $sum_r2/$count,"\t",$sum_dprime/$count,"\n"}else{print "0\t0\n"}')
echo -e "$chr\t$start\t$end\t$ld_data"
done > ww_southern_hifi_bionano.filt.northern.ld.out
```


