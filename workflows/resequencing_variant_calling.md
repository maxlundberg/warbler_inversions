# Resequencing data


## Trim reads

Trim raw reads using `trimmomatic`

```
java -jar trimmomatic.jar PE -threads 16 -basein ${sample}_R1_001.fastq.gz -baseout ${sample}_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
```

## Map reads to assembly

Map reads for each sample using `bwa mem` with a list of samples and read prefixes. Read duplicates will be removed using `picardtools`.

```
genome="ww_southern_hifi_bionano.filt.fasta"
sample_list="sample_info_phyllocopus_new_genome.txt"

bwa index $genome

cat $sample_list | while read name prefix
do
bwa mem $genome -M -t 19 -R "@RG\tID:${name}\tLB:$name\tSM:$name\tPL:Illumina" ${prefix}_[12]P.fastq.gz | samtools view -bS - > $name.bam
samtools sort -@ 20 $name.bam > $name.sorted.bam
java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$name.sorted.bam O=$name.sorted.nodup.bam M=$name.sorted.bam.duplicatedata.txt
samtools index $name.sorted.nodup.bam
cp $name.sorted.nodup.bam* $out_dir
cp $name.sorted.bam.duplicatedata.txt $out_dir
done
```


## Call variants

Call variants using `freebayes`. Use `GNU parallel` to make the variant calling more efficient

```
#create a size-sorted list of scaffolds and remove the mitochondrial scaffold from the list
samtools faidx ww_southern_hifi_bionano.filt.new.fasta
cat ww_southern_hifi_bionano.filt.new.fasta.fai | grep -v mito | sort -k2,2nr | cut -f1 > ww_southern_hifi_bionano.filt.new.scaffolds_excl_mito.list

ls *.bam > bamfiles.list

parallel --jobs 20 "freebayes --bam-list bamfiles.list -f ww_southern_hifi_bionano.filt.new.fasta --region {} | gzip -c > ./vcf_files/{}.vcf.gz; echo {} processed" :::: ww_southern_hifi_bionano.filt.new.scaffolds_excl_mito.list

vcf-concat vcf_files/*.vcf.gz | bgzip -c > freebayes_reseq_southern_hifi_bionano.raw.vcf.gz
tabix -p vcf freebayes_reseq_southern_hifi_bionano.raw.vcf.gz
```


## Filter variants

Use a combination of `vcflib` and `vcftools` to filter the raw set of variants


First, filter variants based on quality, strand support, read placement and genotype coverage

```
vcffilter -f 'QUAL > 30 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0' -g 'DP > 4' freebayes_reseq_southern_hifi_bionano.raw.vcf.gz > freebayes_reseq_southern_hifi_bionano.filt1.vcf
```

Next, remove sites that are missing a maximum of four genotypes in each population

```
vcftools --vcf freebayes_reseq_southern_hifi_bionano.filt1.vcf --keep ~/southern_samples_reseq_new.list --max-missing-count 4 --removed-sites --stdout > removed_sites_southern.txt
vcftools --vcf freebayes_reseq_southern_hifi_bionano.filt1.vcf --keep ~/northern_samples_reseq_new.list --max-missing-count 4 --removed-sites --stdout > removed_sites_northern.txt

cat removed_sites_southern.txt removed_sites_northern.txt | grep -v CHROM | sort -k1,1 -k2,2n | uniq > removed_sites_combined.txt

vcftools --vcf freebayes_reseq_southern_hifi_bionano.filt1.vcf --exclude-positions removed_sites_combined.txt --recode --recode-INFO-all --stdout > freebayes_reseq_southern_hifi_bionano.filt2.vcf
```

Remove sites with a mean coverage twice that of the median mean coverage of all sites (30x) and monomorphic sites

```
vcftools --vcf freebayes_reseq_southern_hifi_bionano.filt2.vcf --max-meanDP 30 --removed-sites --stdout > excessive_coverage_positions.out

vcftools --vcf freebayes_reseq_southern_hifi_bionano.filt2.vcf --freq --stdout | sed '1d' | grep ":1" > monomorphic_pos.out

cat excessive_coverage_positions.out monomorphic_pos.out | grep -v CHROM | cut -f1-2 | sort -k1,1 -k2,2n | uniq > excessive_coverage_monomorphic_positions.out 

vcftools --vcf freebayes_reseq_southern_hifi_bionano.filt2.vcf --exclude-positions excessive_coverage_monomorphic_positions.out --recode --recode-INFO-all --stdout | bgzip -c > freebayes_reseq_southern_hifi_bionano.filt3.vcf.gz
tabix -p vcf freebayes_reseq_southern_hifi_bionano.filt3.vcf.gz
```

Decompose variants and remove variants overlapping annotated repeats

```
vcfallelicprimitives -kg freebayes_reseq_southern_hifi_bionano.filt3.vcf.gz | bgzip -c > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.vcf.gz
tabix -p vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.vcf.gz

echo "#Scaffold start end" > header.txt
cat header.txt ww_southern_hifi_bionano.filt.fasta.bed > ww_southern_hifi_bionano.filt.fasta.new.bed

vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.vcf.gz --recode --recode-INFO-all --exclude-bed  ww_southern_hifi_bionano.filt.fasta.new.bed --stdout | bgzip -c > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.vcf.gz
tabix -p vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.vcf.gz
```


For several analyses later on, such as fst calculations, we will need to fix how a small number of missing genotypes have been recoded in vcftools

```
zcat freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.vcf.gz | perl -ne 'if($_=~/\t\.:/){$_=~s/\t\.:/\t\.\/\.:/g} print $_' | bgzip -c > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz
tabix -p vcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz
```


## Calculate FST

To calculate FST and extract highly differentiated positions we will use `vcftools`

```
fst_thresholds="0.5 0.6 0.7 0.8 0.9 1"
for fst_threshold in $fst_thresholds
do
export fst_threshold
zcat freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz | vcftools --vcf - --weir-fst-pop ~/southern_samples_reseq_new.list  --weir-fst-pop ~/northern_samples_reseq_new.list --stdout | sed '1d' | perl -lane 'if($F[2]>=$ENV{fst_threshold}){print $F[0],"\t",$F[1]-1,"\t",$F[1],"\t",$F[2]}' > fst_${fst_threshold}_variants.bed
done
```

Also calculate FST for bi-allelic SNPs in non-overlapping 10 kb windows and also do the same for variants with a minor allele frequency (MAF) of at least 0.1

```
vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz --weir-fst-pop ~/southern_samples_reseq_new.list --weir-fst-pop ~/northern_samples_reseq_new.list --fst-window-size 10000 --stdout --max-alleles 2 --remove-indels > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.fst_10kb.out
vcftools --gzvcf freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz --weir-fst-pop ~/southern_samples_reseq_new.list --weir-fst-pop ~/northern_samples_reseq_new.list --fst-window-size 10000 --maf 0.1 --stdout --max-alleles 2 --remove-indels > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.fst_10kb.maf0.1.out

```

Finally, get FST estimates for bi-allelic SNPs with a MAF of at least 0.1 for scaffolds containing the divergent regions and predicted adjacent scaffolds

```
scaffolds="Scaffold11 Scaffold19 Scaffold12 Scaffold61 Scaffold38 Scaffold0"
for scaffold in $scaffolds
do
bcftools view freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz $scaffold | vcftools --vcf - --weir-fst-pop ~/southern_samples_reseq_new.list --weir-fst-pop ~/northern_samples_reseq_new.list --maf 0.1 --stdout 
done > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.fst_div_region_adj_scaffolds.fst.out

(head -n1 freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.fst_div_region_adj_scaffolds.fst.out && grep -v POS freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.fst_div_region_adj_scaffolds.fst.out) > freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.fst_div_region_adj_scaffolds.fst.2.out 
```



## Annotate highly differentiated variants in the divergent regions

Extract variants that have a FST of at least 0.7 between homozygotes of northern and southern haplotypes in each of the divergent regions
For this purpose use `vcftools` and `bcftools` to efficiently retrieve variants from each region 

```
mkdir highly_diff_variants_div_region


#Chr1

bcftools view freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz  Scaffold19 | vcftools --vcf - --weir-fst-pop ~/chr1_southern.txt  --weir-fst-pop ~/chr1_northern.txt --stdout | perl -lane 'if($F[2]>=0.7){print $F[0],"\t",$F[1],"\t",$F[2]}' > highly_diff_variants_div_region/chr1_highly_diff_positions.fst.out
cat highly_diff_variants_div_region/chr1_highly_diff_positions.fst.out | cut -f1-2 > highly_diff_variants_div_region/chr1_highly_diff_positions.list

bcftools view freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz  Scaffold19 | vcftools --vcf - --positions highly_diff_variants_div_region/chr1_highly_diff_positions.list --recode --recode-INFO-all --stdout > highly_diff_variants_div_region/chr1_highly_differentiated_variants.vcf



#Chr3
bcftools view freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz  Scaffold61:56000000- | vcftools --vcf - --weir-fst-pop ~/chr3_southern.txt  --weir-fst-pop ~/chr3_northern.txt --stdout | perl -lane 'if($F[2]>=0.7){print $F[0],"\t",$F[1],"\t",$F[2]}' > highly_diff_variants_div_region/chr3_highly_diff_positions.fst.out
cat highly_diff_variants_div_region/chr3_highly_diff_positions.fst.out | cut -f1-2 > highly_diff_variants_div_region/chr3_highly_diff_positions.list

bcftools view freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz  Scaffold61:56000000- | vcftools --vcf - --positions highly_diff_variants_div_region/chr3_highly_diff_positions.list --recode --recode-INFO-all --stdout > highly_diff_variants_div_region/chr3_highly_differentiated_variants.vcf


#Chr5
bcftools view freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz  Scaffold0:200000-4360000 | vcftools --vcf - --weir-fst-pop ~/chr5_southern.txt  --weir-fst-pop ~/chr5_northern.txt --stdout | perl -lane 'if($F[2]>=0.7){print $F[0],"\t",$F[1],"\t",$F[2]}' > highly_diff_variants_div_region/chr5_highly_diff_positions.fst.out
cat highly_diff_variants_div_region/chr5_highly_diff_positions.fst.out | cut -f1-2 > highly_diff_variants_div_region/chr5_highly_diff_positions.list

bcftools view freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.vcf.gz  Scaffold0:200000-4360000 | vcftools --vcf - --positions highly_diff_variants_div_region/chr5_highly_diff_positions.list --recode --recode-INFO-all --stdout > highly_diff_variants_div_region/chr5_highly_differentiated_variants.vcf
```

Use `snpeff` and `snpsift` to annotate the variants. As for annotation we will use manually curated genes from the divergent regions 

```
mkdir data
cd data
mkdir ww genomes

cd genomes
cp /../../ww_southern_hifi_bionano.filt.new.fasta . 
mv ww_southern_hifi_bionano.filt.new.fasta ww.fa

cd ../ww
cp /../../ww_hifi_southern_annotations_webapollo_20220208.gff3 genes.gff

java -jar snpEff.jar build -gff3 -v ww 


java -jar ../../annotation_data/snpEff/snpEff.jar ann ww chr1_highly_differentiated_variants.vcf  > chr1_highly_differentiated_variants.annotated.vcf
mv snpEff_genes.txt chr1_snpEff_genes.txt
mv snpEff_summary.html chr1_snpEff_summary.html

java -jar ../../annotation_data/snpEff/snpEff.jar ann ww chr3_highly_differentiated_variants.vcf  > chr3_highly_differentiated_variants.annotated.vcf
mv snpEff_genes.txt chr3_snpEff_genes.txt
mv snpEff_summary.html chr3_snpEff_summary.html

java -jar ../../annotation_data/snpEff/snpEff.jar ann ww chr5_highly_differentiated_variants.vcf  > chr5_highly_differentiated_variants.annotated.vcf
mv snpEff_genes.txt chr5_snpEff_genes.txt
mv snpEff_summary.html chr5_snpEff_summary.html

cat chr1_highly_differentiated_variants.annotated.vcf | ../../annotation_data/snpEff/scripts/vcfEffOnePerLine.pl | java -jar ../../annotation_data/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT ANN[*].GENE ANN[*].EFFECT ANN[*].IMPACT ANN[*].ERRORS | egrep "MODERATE|HIGH" | uniq > moderate_high_effect_variants.list 
cat chr3_highly_differentiated_variants.annotated.vcf | ../../annotation_data/snpEff/scripts/vcfEffOnePerLine.pl | java -jar ../../annotation_data/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT ANN[*].GENE ANN[*].EFFECT ANN[*].IMPACT ANN[*].ERRORS | egrep "MODERATE|HIGH" | uniq  >> moderate_high_effect_variants.list 
cat chr5_highly_differentiated_variants.annotated.vcf | ../../annotation_data/snpEff/scripts/vcfEffOnePerLine.pl | java -jar ../../annotation_data/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT ANN[*].GENE ANN[*].EFFECT ANN[*].IMPACT ANN[*].ERRORS | egrep "MODERATE|HIGH" | uniq  >> moderate_high_effect_variants.list 

```


## Call structural variants in the resequenced samples

First get a high-confidence set of variants from long-read alignments by using `delly`. For details on alignments see the breakpoint workflow document  

```
delly_v0.9.1_linux_x86_64bit lr -g ww_southern_hifi_bionano.filt.new.fasta --svtype ALL ww_southern_hifi_bionano.filt.fasta.ww_northern.pacbio_reads.sorted.bam -y pb --outfile delly.northern_pacbio.vs.southern_hifi_genome.bcf
```

Next, use graphtyper to genotype the resequenced willow warblers for the long read

```
bcftools view delly.northern_pacbio.vs.southern_hifi_genome.bcf | bgzip -c > $prefix.vcf.gz
tabix -p vcf $prefix.vcf.gz 

echo "Scaffold0" > regions.txt
echo "Scaffold19" >> regions.txt
echo "Scaffold61" >> regions.txt

graphtyper genotype_sv ww_southern_hifi_bionano.filt.new.fasta delly.northern_pacbio.vs.southern_hifi_genome.vcf.gz --sams=bam_files.list --region_file=regions.txt --threads=20 --output=graphtyper_northern_pb_southern_hifi_div_regions_scaffolds

find graphtyper_northern_pb_southern_hifi_div_regions_scaffolds -name "*.vcf.gz" > graphtyper_northern_pb_southern_hifi.list

bcftools concat --naive --file-list graphtyper_northern_pb_southern_hifi.list -Oz -o graphtyper_northern_pb_southern_hifi.vcf.gz

vcftools --gzvcf graphtyper_northern_pb_southern_hifi.vcf.gz --recode --recode-INFO-all --stdout --remove-filtered-all | gzip -c > graphtyper_northern_pb_southern_hifi.filt.vcf.gz

#Chr1

zcat graphtyper_northern_pb_southern_hifi.filt.vcf.gz | vcftools --vcf - --chr Scaffold19 --weir-fst-pop ~/chr1_southern.txt --weir-fst-pop ~/chr1_northern.txt --stdout | grep -v "nan" | awk '{if($3>=0.7){print $0}}' | sort -k2,2n > graphtyper_northern_pb_southern_hifi.filt.chr1.fst.0.7.out

cat graphtyper_northern_pb_southern_hifi.filt.chr1.fst.0.7.out | sed '1d' | cut -f1-2 | uniq > graphtyper_northern_pb_southern_hifi.filt.chr1.fst.0.7.list #N=43

#Find positions where at least 80 % of the genotypes among southern and nothern birds have a genotype 
vcftools --gzvcf graphtyper_northern_pb_southern_hifi.vcf.gz --positions graphtyper_northern_pb_southern_hifi.filt.chr1.fst.0.7.list --keep ~/chr1_southern.txt --max-missing 0.8 --recode --stdout | grep -v "#" | cut -f1-2 | uniq > graphtyper_northern_pb_southern_hifi.chr1.pos_ok_southern.list
vcftools --gzvcf graphtyper_northern_pb_southern_hifi.vcf.gz --positions graphtyper_northern_pb_southern_hifi.filt.chr1.fst.0.7.list --keep ~/chr1_northern.txt --max-missing 0.8 --recode --stdout | grep -v "#" | cut -f1-2 | uniq > graphtyper_northern_pb_southern_hifi.chr1.pos_ok_northern.list

cat graphtyper_northern_pb_southern_hifi.chr1.pos_ok_southern.list graphtyper_northern_pb_southern_hifi.chr1.pos_ok_northern.list | sort -k2,2n | uniq -d > graphtyper_northern_pb_southern_hifi.chr1.pos_ok_northern_southern.list #N=25

#Only keep variants with aggregated calls (N=20)
vcftools --gzvcf graphtyper_northern_pb_southern_hifi.vcf.gz --positions graphtyper_northern_pb_southern_hifi.chr1.pos_ok_northern_southern.list --recode --recode-INFO-all --stdout | grep -v "#" | grep "AGGREGATED" | perl -lane '$F[7]=~/END=(\d+).+SVSIZE=(\d+).+SVTYPE=([^;]+)/; print $F[0],"\t",$F[1]-1,"\t",$1,"\t",$3,"-",$2' | sort -k2,2n > graphtyper.chr1.high_diff_svs.bed 

#Chr3

zcat graphtyper_northern_pb_southern_hifi.filt.vcf.gz | vcftools --vcf - --chr Scaffold61 --from-bp 56069986 --weir-fst-pop ~/chr3_southern.txt --weir-fst-pop ~/chr3_northern.txt --stdout | grep -v "nan" | awk '{if($3>=0.7){print $0}}' | sort -k2,2n > graphtyper_northern_pb_southern_hifi.filt.chr3.fst.0.7.out

cat graphtyper_northern_pb_southern_hifi.filt.chr3.fst.0.7.out | sed '1d' | cut -f1-2 | uniq > graphtyper_northern_pb_southern_hifi.filt.chr3.fst.0.7.list

vcftools --gzvcf graphtyper_northern_pb_southern_hifi.vcf.gz --positions graphtyper_northern_pb_southern_hifi.filt.chr3.fst.0.7.list --keep ~/chr3_southern.txt --max-missing 0.8 --recode --stdout | grep -v "#" | cut -f1-2 | uniq > graphtyper_northern_pb_southern_hifi.chr3.pos_ok_southern.list
vcftools --gzvcf graphtyper_northern_pb_southern_hifi.vcf.gz --positions graphtyper_northern_pb_southern_hifi.filt.chr3.fst.0.7.list --keep ~/chr3_northern.txt --max-missing 0.8 --recode --stdout | grep -v "#" | cut -f1-2 | uniq > graphtyper_northern_pb_southern_hifi.chr3.pos_ok_northern.list

cat graphtyper_northern_pb_southern_hifi.chr3.pos_ok_southern.list graphtyper_northern_pb_southern_hifi.chr3.pos_ok_northern.list | sort -k2,2n | uniq -d > graphtyper_northern_pb_southern_hifi.chr3.pos_ok_northern_southern.list #N=36

vcftools --gzvcf graphtyper_northern_pb_southern_hifi.vcf.gz --positions graphtyper_northern_pb_southern_hifi.chr3.pos_ok_northern_southern.list --recode --recode-INFO-all --stdout | grep -v "#" | grep "AGGREGATED" | perl -lane '$F[7]=~/END=(\d+).+SVSIZE=(\d+).+SVTYPE=([^;]+)/; print $F[0],"\t",$F[1]-1,"\t",$1,"\t",$3,"-",$2' | sort -k2,2n > graphtyper.chr3.high_diff_svs.bed #N=25

#Chr5

zcat graphtyper_northern_pb_southern_hifi.filt.vcf.gz | vcftools --vcf - --chr Scaffold0 --from-bp 327097 --to-bp 4361809 --weir-fst-pop ~/chr5_southern.txt --weir-fst-pop ~/chr5_northern.txt --stdout | grep -v "nan" | awk '{if($3>=0.7){print $0}}' | sort -k2,2n > graphtyper_northern_pb_southern_hifi.filt.chr5.fst.0.7.out

cat graphtyper_northern_pb_southern_hifi.filt.chr5.fst.0.7.out | sed '1d' | cut -f1-2 | uniq > graphtyper_northern_pb_southern_hifi.filt.chr5.fst.0.7.list

vcftools --gzvcf graphtyper_northern_pb_southern_hifi.vcf.gz --positions graphtyper_northern_pb_southern_hifi.filt.chr5.fst.0.7.list --keep ~/chr5_southern.txt --max-missing 0.8 --recode --stdout | grep -v "#" | cut -f1-2 | uniq > graphtyper_northern_pb_southern_hifi.chr5.pos_ok_southern.list
vcftools --gzvcf graphtyper_northern_pb_southern_hifi.vcf.gz --positions graphtyper_northern_pb_southern_hifi.filt.chr5.fst.0.7.list --keep ~/chr5_northern.txt --max-missing 0.8 --recode --stdout | grep -v "#" | cut -f1-2 | uniq > graphtyper_northern_pb_southern_hifi.chr5.pos_ok_northern.list

cat graphtyper_northern_pb_southern_hifi.chr5.pos_ok_southern.list graphtyper_northern_pb_southern_hifi.chr5.pos_ok_northern.list | sort -k2,2n | uniq -d > graphtyper_northern_pb_southern_hifi.chr5.pos_ok_northern_southern.list

vcftools --gzvcf graphtyper_northern_pb_southern_hifi.vcf.gz --positions graphtyper_northern_pb_southern_hifi.chr5.pos_ok_northern_southern.list --recode --recode-INFO-all --stdout | grep -v "#" | grep "AGGREGATED" | perl -lane '$F[7]=~/END=(\d+).+SVSIZE=(\d+).+SVTYPE=([^;]+)/; print $F[0],"\t",$F[1]-1,"\t",$1,"\t",$3,"-",$2' | sort -k2,2n > graphtyper.chr5.high_diff_svs.bed #N=12

cat graphtyper.chr*.high_diff_svs.bed | sort -k1,1 -k2,2n > graphtyper.combined.high_diff_svs.bed
```

Check distance of the SVs to the closest annotated gene and determine if any SVs overlap with exons 

```
bedtools closest -a graphtyper.combined.high_diff_svs.bed -b ww_hifi_southern_annotations_webapollo_20220208.bed -d  > graphtyper.combined.high_diff_svs.gene_distance.bed

bedtools intersect -a graphtyper.combined.high_diff_svs.bed -b ww_hifi_southern_annotations_webapollo_20220208.exon_intervals.bed -wo > graphtyper.combined.high_diff_svs.exon_overlap.bed
```



## Get coverage for resequenced samples

Get coverage in 1 kb windows for the resequencing samples using `BEDTools`

```
cat ww_southern_hifi_bionano.filt.new.fasta.fai | cut -f1-2 > ww_southern_hifi_bionano.filt.new.fasta.scaffold_sizes.list
bedtools makewindows -g ../ww_southern_hifi_bionano.filt.new.fasta.scaffold_sizes.list -w 1000 > ww_southern_hifi_bionano.filt.new.fasta.1kb_windows.bed 
bedtools multicov -p -q 1 -bed $bed_file -bams 1A05.sorted.nodup.bam UK06.sorted.nodup.bam 0G03.sorted.dedup.bam 0G04.sorted.dedup.bam 0G10.sorted.dedup.bam 0J01.sorted.dedup.bam 1P02.sorted.dedup.bam 3K06.sorted.dedup.bam 7A12.sorted.dedup.bam 96A01.sorted.dedup.bam 96B07.sorted.dedup.bam 1M13.sorted.nodup.bam 1N12.sorted.nodup.bam 1K05.sorted.dedup.bam 1K10.sorted.dedup.bam 1L17.sorted.dedup.bam 1L19.sorted.dedup.bam 1L20.sorted.dedup.bam 1M08.sorted.dedup.bam 1O01.sorted.dedup.bam 1O04.sorted.dedup.bam 1O06.sorted.dedup.bam > ww_read_counts.southern_hifi.1kb.windows.out 
```


## Get divergent region haplotypes for resequenced samples

To determine the divergent region genotypes (southern/northern haplotypes) for each resequenced sample we extracted genotypes for SNPs present on the SNP array
in Lundberg et al. 2017.

First map the SNP array probes (final_probes.fasta) to the northern assembly using `gmap`

```
gmap_build -D ./ -d ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta

gmap -d ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt -D ./ -t5 -n1 -f coords final_probes.fasta > final_probes.new_genome.gmap.coords.out
```

Extract the position of the targeted SNPs from the alignments of the probes e

```
cat final_probes.new_genome.gmap.coords.out | extract_ref_pos_snps_array_gmap_coords.pl > final_probes.new_genome.gmap.coords.pos.info.out
cat final_probes.new_genome.gmap.coords.pos.info.out | grep -v "#" | cut -f3 > mapped_snps.txt 
```

Filter the original SNP array data file using `plink` as in Lundberg et al 2017 and only select SNPs with probes that map in the new genome

```
plink --file 2014-01_All_12_Plates_mod  --autosome-num 40 --maf 0.01 --geno 0.05 --mind 0.05 --exclude loci_to_remove.txt --extract mapped_snps.txt --remove ind_to_remove.txt --recode --out trimmed_new --allow-extra-chr --update-map /snps_plink_update_position.txt
```


Adjust the major allele to the reference allele, flip genotypes of SNPs found on the reverse strand and convert the data into a vcf file

```
plink --file trimmed_new --recode --allow-extra-chr --autosome-num 40 --a2-allele final_probes.new_genome.gmap.coords.pos.info.out 4 3 '#' --out trimmed_new_ref_flipped

cat trimmed_new_ref_flipped.log | perl -ne 'if($_=~/variant (\S+)\.$/){print "$1\n"}' > snps_to_flip.list

plink --file trimmed_new --recode --flip snps_to_flip.list --allow-extra-chr --autosome-num 40 --a2-allele final_probes.new_genome.gmap.coords.pos.info.out 4 3 '#' --out trimmed_new_ref_flipped_2

plink --file trimmed_new_ref_flipped_2 --allow-extra-chr --autosome-num 40 --out snparray_data_filt --update-chr final_probes.new_genome.gmap.coords.pos.info.out 1 3 '#' --update-map final_probes.new_genome.gmap.coords.pos.info.out 2 3 '#' --make-bed

plink --bfile snparray_data_filt --recode vcf --allow-extra-chr --out snparray_data_filt

cat snparray_data_filt.vcf | perl -ne 'if($_=~/CHROM/){$_=~s/\d+\_//g}; print $_' > snparray_data_filt.new.vcf
```

Genotype the resequenced samples and the two linked read samples for these variants using `freebayes` and use `vcflib` and `vcftools` to decompose and select the correct variants.
We will also remove eight SNPs that have more than one alternative allele and convert the data into plink format.

```
cat snparray_data_filt.new.vcf | grep -v "#" | perl -lane 'if($F[0] eq "Scaffold156" || $F[0] eq "Scaffold29" && $F[1]>55759500 || $F[0] eq "Scaffold68"){print $F[0],"\t",$F[1]-1,"\t",$F[1]}' > snparray_data_filt.new.pos.div.region.bed 

cat snparray_data_filt.new.pos.div.region.bed | cut -f1,3 > snparray_data_filt.new.pos.div.region.out

freebayes -f ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta --report-monomorphic --targets snparray_data_filt.new.pos.div.region.bed -L ~/bam_files_reseq_data.list > snp_array_snps_div_regions_reseq_10x_samples.vcf

vcfallelicprimitives -kg snp_array_snps_div_regions_reseq_10x_samples.vcf | vcftools --vcf - --positions snparray_data_filt.new.pos.div.region.out --recode --stdout | sed 's/_102/-102/;s/_101_uppmax_longranger_wgs/-101/' > snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.vcf

cat snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.vcf | cut -f1-5 | grep -v "#" | grep "," | cut -f1-2 > snp_positions_remove.list

vcftools --vcf snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.vcf  --recode --stdout --exclude-positions snp_positions_remove.list | add_snp_id_from_other_vcf.pl --ref_vcf_file snparray_data_filt.new.vcf > snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.new.vcf

plink --vcf snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.new.vcf --recode --allow-extra-chr --out snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.new
```

Combine the data from the SNParray with that of the resequenced samples

```
cat snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.new.vcf  | grep -v "#" | cut -f3 > snps_div_region_reseq_chromium.list

plink --bfile snparray_data_filt --extract snps_div_region_reseq_chromium.list --merge snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.new.ped  snp_array_snps_div_regions_reseq_10x_samples.decomposed.filt.new.map --recode --out merged_data --allow-extra-chr

plink --noweb --file merged_data --make-bed --out merged_bin --allow-extra-chr
```

Use `InvClust` in R to get inversion genotypes in the three divergent regions 


```
plinkdata=read.plink("merged_bin.bed","merged_bin.bim","merged_bin.fam")
plink_genos=plinkdata$genotypes
plink_annot=plinkdata$map

genos=as.data.frame(array(NA,c(nrow(plink_genos),4)))
genos[,1]=paste(plinkdata$fam[,2])

genos2=as.data.frame(array(NA,c(nrow(plink_genos),5)))
genos2[,1]=paste(plinkdata$fam[,2])


#chr1
roi=data.frame(chr="Scaffold156",LBP=1,RBP=12000000,reg= "inv1")
invcall=invClust(roi=roi,wh=1,geno=plink_genos,annot=plink_annot,dim=2)
values=invcall[[2]]$y[,1:2]
samples_plink=plinkdata$fam[,2]
write.table(cbind(samples_plink,values),"invclust_mds_coordinates_chrom1.txt",sep="\t",row.names=F,col.names=F,quote=F)

genos[,2]=paste(invGenotypes(invcall))


#chr3
roi=data.frame(chr="Scaffold29",LBP=55E6,RBP=80E6,reg= "inv3")
invcall=invClust(roi=roi,wh=1,geno=plink_genos,annot=plink_annot,dim=2)
values=invcall[[2]]$y[,1:2]
samples_plink=plinkdata$fam[,2]
write.table(cbind(samples_plink,values),"invclust_mds_coordinates_chrom3.txt",sep="\t",row.names=F,col.names=F,quote=F)

genos[,3]=paste(invGenotypes(invcall))

#chr5
roi=data.frame(chr="Scaffold68",LBP=1,RBP=6000000,reg= "inv5")
invcall=invClust(roi=roi,wh=1,geno=plink_genos,annot=plink_annot,dim=2)
values=invcall[[2]]$y[,1:2]
samples_plink=plinkdata$fam[,2]
write.table(cbind(samples_plink,values),"invclust_mds_coordinates_chrom5.txt",sep="\t",row.names=F,col.names=F,quote=F)


genos[,4]=paste(invGenotypes(invcall))

write.table(genos,file="invclust_genotypes_snparray_resequencing_new_samples.txt",sep="\t",row.names=F,col.names=F,quote=F)
```










