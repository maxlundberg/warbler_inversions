# Check reference bias

Here we will explore whether there is a reference bias when using the southern or the northern assembly as reference.
To this end, we will also map, call variants and calculate genetic diversity statistics using the northern reference


## Map reads to northern assembly

Map reads for each sample using `bwa mem` with a list of samples and read prefixes. Read duplicates will be removed using `picardtools`.

```
genome="ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta"
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
cat ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta.fai | grep -v MT | sort -k2,2nr | cut -f1 > ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta.scaffolds.list

ls *.bam > bamfiles.list

parallel --jobs 20 "freebayes --bam-list bamfiles.list -f ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta --region {} | gzip -c > ./vcf_files/{}.vcf.gz; echo {} processed" :::: ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta.scaffolds.list

vcf-concat vcf_files/*.vcf.gz | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.raw.vcf.gz
tabix -p vcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.raw.vcf.gz
```


## Filter variants

Use a combination of `vcflib` and `vcftools` to filter the raw set of variants


First, filter variants based on quality, strand support and read placement

```
vcffilter -f 'QUAL > 30 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0' freebayes_old_reseq_plus_new_samples_pacbio_bionano.raw.vcf.gz | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt1.vcf.gz&
tabix -p vcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt1.vcf.gz
```

Next, filter based on minimum coverage and remove sites that are missing a maximum of four genotypes in each population

```
vcffilter -g "DP > 4" freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt1.vcf.gz | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.vcf.gz

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.vcf.gz --keep ~/southern_samples_reseq_new.list --max-missing-count 4 --removed-sites --stdout > removed_sites_southern.txt #5454894 sites
vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.vcf.gz --keep ~/northern_samples_reseq_new.list --max-missing-count 4 --removed-sites --stdout > removed_sites_northern.txt #2328213 sites

cat removed_sites_southern.txt removed_sites_northern.txt | grep -v CHROM | sort -k1,1 -k2,2n | uniq > removed_sites_combined.txt  #5735154 sites

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.vcf.gz --exclude-positions removed_sites_combined.txt --recode --recode-INFO-all --stdout | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.vcf.gz&
```

Remove sites with a mean coverage twice that of the median mean coverage of all sites and monomorphic sites

```
vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.vcf.gz --site-mean-depth --stdout  --not-chr MT > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.mean_depths.out&

#Get median in R
depths=read.delim("freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.mean_depths.out")
median(depths[,3])
[1] 15.1818
#

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.vcf.gz --max-meanDP 30 --removed-sites --stdout > excessive_coverage_positions.out

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.vcf.gz --freq --stdout | sed '1d' | grep ":1" > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.monomorphic_pos.out

cat excessive_coverage_positions.out freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.monomorphic_pos.out | grep -v CHROM | cut -f1-2 | sort -k1,1 -k2,2n | uniq > excessive_coverage_monomorphic_positions.out 

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.vcf.gz --exclude-positions excessive_coverage_monomorphic_positions.out --recode --recode-INFO-all --stdout | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.vcf.gz
tabix -p vcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.vcf.gz
```

Decompose variants and remove variants overlapping annotated repeats

```
vcfallelicprimitives -kg freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.vcf.gz | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.vcf.gz
tabix -p vcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.vcf.gz

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.vcf.gz --recode --recode-INFO-all --exclude-bed ww_pacbio_bionano_final.repeats.corrected.new.bed --stdout | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.vcf.gz
tabix -p vcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.vcf.gz
```


For several analyses later on, such as fst calculations, we will need to fix how a small number of missing genotypes have been recoded in vcftools

```
zcat freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.vcf.gz | perl -ne 'if($_=~/\t\.:/){$_=~s/\t\.:/\t\.\/\.:/g} print $_' | bgzip -c > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.vcf.gz
tabix -p vcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.vcf.gz
```


## Phase variants and calculate genetic summary statistics 

Use `beagle` to phase variants 

```
#Exclude scaffolds with only one variants, otherwise the software will not work properly
echo "Scaffold158:57511" > exclude_markers.list
echo "Scaffold172:22204" >> exclude_markers.list
echo "Scaffold215:1410" >> exclude_markers.list
echo "Scaffold279:266" >> exclude_markers.list
echo "Scaffold284:3381" >> exclude_markers.list
echo "Scaffold331:72396" >> exclude_markers.list

java -Xmx240G -jar ~/beagle.22Jul22.46e.jar gt=freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.vcf.gz out=freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2 excludemarkers=exclude_markers.list
```


Calculate pi and Taj'D in 10 kb windows using `vcftools` 

```
vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2.vcf.gz  --keep ~/ww_samples_pure_northern.txt  --remove-indels --max-alleles 2 --TajimaD 10000 --stdout > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2.northern.tajd_10kb.out
vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2.vcf.gz  --keep ~/ww_samples_pure_northern.txt  --remove-indels --max-alleles 2 --window-pi 10000 --stdout > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2.northern.pi_10kb.out

vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2.vcf.gz  --keep ~/ww_samples_pure_southern.txt  --remove-indels --max-alleles 2 --TajimaD 10000 --stdout > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2.southern.tajd_10kb.out
vcftools --gzvcf freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2.vcf.gz  --keep ~/ww_samples_pure_southern.txt  --remove-indels --max-alleles 2 --window-pi 10000 --stdout > freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2.southern.pi_10kb.out
```


## Compare genetic summary statistics between the genomes

Summarize the diversity data in `R`

```
#In R

north_tajd.1=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.tajd_10kb.out",header=T)
south_tajd.1=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.tajd_10kb.out",header=T)

north_tajd.2=read.delim("freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2.northern.tajd_10kb.out",header=T)
south_tajd.2=read.delim("freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2.southern.tajd_10kb.out",header=T)

south_pi.1=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pi_10kb.out",header=T)
north_pi.1=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pi_10kb.out",header=T)

south_pi.2=read.delim("freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2.southern.pi_10kb.out",header=T)
north_pi.2=read.delim("freebayes_old_reseq_plus_new_samples_pacbio_bionano.filt2.min_depth.missingness.filt.decomposed.no_repeat.new.phased.2.northern.pi_10kb.out",header=T)
```


Southern assembly

```
chrom3_south_start=56E6
chrom5_south_start=3.27E5
chrom5_south_end=4.36E6

#North pi
north_pi.1_chr1=north_pi.1[which(north_pi.1$CHROM=="Scaffold19"),]
north_pi.1_rest=north_pi.1[-which(north_pi.1$CHROM=="Scaffold19"),]
north_pi.1_chr3=north_pi.1_rest[which(north_pi.1_rest$CHROM=="Scaffold61" & north_pi.1_rest$BIN_START>chrom3_south_start),]
north_pi.1_rest=north_pi.1_rest[-which(north_pi.1_rest$CHROM=="Scaffold61" & north_pi.1_rest$BIN_START>chrom3_south_start),]

north_pi.1_chr5=north_pi.1_rest[which(north_pi.1_rest$CHROM=="Scaffold0" & north_pi.1_rest$BIN_START<chrom5_south_end & north_pi.1_rest$BIN_START>chrom5_south_start),]
north_pi.1_rest=north_pi.1_rest[-which(north_pi.1_rest$CHROM=="Scaffold0" & north_pi.1_rest$BIN_START<chrom5_south_end & north_pi.1_rest$BIN_START>chrom5_south_start),]

north_pi.1_chr1=cbind(rep("chr1",nrow(north_pi.1_chr1)),north_pi.1_chr1)
names(north_pi.1_chr1)[1]="region"
north_pi.1_chr3=cbind(rep("chr3",nrow(north_pi.1_chr3)),north_pi.1_chr3)
names(north_pi.1_chr3)[1]="region"
north_pi.1_chr5=cbind(rep("chr5",nrow(north_pi.1_chr5)),north_pi.1_chr5)
names(north_pi.1_chr5)[1]="region"
north_pi.1_rest=cbind(rep("rest_genome",nrow(north_pi.1_rest)),north_pi.1_rest)
names(north_pi.1_rest)[1]="region"
north_pi.1_comb=rbind(north_pi.1_chr1,north_pi.1_chr3,north_pi.1_chr5,north_pi.1_rest)


#South pi
south_pi.1_chr1=south_pi.1[which(south_pi.1$CHROM=="Scaffold19"),]
south_pi.1_rest=south_pi.1[-which(south_pi.1$CHROM=="Scaffold19"),]
south_pi.1_chr3=south_pi.1_rest[which(south_pi.1_rest$CHROM=="Scaffold61" & south_pi.1_rest$BIN_START>chrom3_south_start),]
south_pi.1_rest=south_pi.1_rest[-which(south_pi.1_rest$CHROM=="Scaffold61" & south_pi.1_rest$BIN_START>chrom3_south_start),]

south_pi.1_chr5=south_pi.1_rest[which(south_pi.1_rest$CHROM=="Scaffold0" & south_pi.1_rest$BIN_START<chrom5_south_end & south_pi.1_rest$BIN_START>chrom5_south_start),]
south_pi.1_rest=south_pi.1_rest[-which(south_pi.1_rest$CHROM=="Scaffold0" & south_pi.1_rest$BIN_START<chrom5_south_end & south_pi.1_rest$BIN_START>chrom5_south_start),]

south_pi.1_chr1=cbind(rep("chr1",nrow(south_pi.1_chr1)),south_pi.1_chr1)
names(south_pi.1_chr1)[1]="region"
south_pi.1_chr3=cbind(rep("chr3",nrow(south_pi.1_chr3)),south_pi.1_chr3)
names(south_pi.1_chr3)[1]="region"
south_pi.1_chr5=cbind(rep("chr5",nrow(south_pi.1_chr5)),south_pi.1_chr5)
names(south_pi.1_chr5)[1]="region"
south_pi.1_rest=cbind(rep("rest_genome",nrow(south_pi.1_rest)),south_pi.1_rest)
names(south_pi.1_rest)[1]="region"
south_pi.1_comb=rbind(south_pi.1_chr1,south_pi.1_chr3,south_pi.1_chr5,south_pi.1_rest)

#North Tajima's D
north_tajd.1_chr1=north_tajd.1[which(north_tajd.1$CHROM=="Scaffold19"),]
north_tajd.1_rest=north_tajd.1[-which(north_tajd.1$CHROM=="Scaffold19"),]
north_tajd.1_chr3=north_tajd.1_rest[which(north_tajd.1_rest$CHROM=="Scaffold61" & north_tajd.1_rest$BIN_START>chrom3_south_start),]
north_tajd.1_rest=north_tajd.1_rest[-which(north_tajd.1_rest$CHROM=="Scaffold61" & north_tajd.1_rest$BIN_START>chrom3_south_start),]

north_tajd.1_chr5=north_tajd.1_rest[which(north_tajd.1_rest$CHROM=="Scaffold0" & north_tajd.1_rest$BIN_START<chrom5_south_end & north_tajd.1_rest$BIN_START>chrom5_south_start),]
north_tajd.1_rest=north_tajd.1_rest[-which(north_tajd.1_rest$CHROM=="Scaffold0" & north_tajd.1_rest$BIN_START<chrom5_south_end & north_tajd.1_rest$BIN_START>chrom5_south_start),]

north_tajd.1_chr1=cbind(rep("chr1",nrow(north_tajd.1_chr1)),north_tajd.1_chr1)
names(north_tajd.1_chr1)[1]="region"
north_tajd.1_chr3=cbind(rep("chr3",nrow(north_tajd.1_chr3)),north_tajd.1_chr3)
names(north_tajd.1_chr3)[1]="region"
north_tajd.1_chr5=cbind(rep("chr5",nrow(north_tajd.1_chr5)),north_tajd.1_chr5)
names(north_tajd.1_chr5)[1]="region"
north_tajd.1_rest=cbind(rep("rest_genome",nrow(north_tajd.1_rest)),north_tajd.1_rest)
names(north_tajd.1_rest)[1]="region"
north_tajd.1_comb=rbind(north_tajd.1_chr1,north_tajd.1_chr3,north_tajd.1_chr5,north_tajd.1_rest)



#South Tajima's D
south_tajd.1_chr1=south_tajd.1[which(south_tajd.1$CHROM=="Scaffold19"),]
south_tajd.1_rest=south_tajd.1[-which(south_tajd.1$CHROM=="Scaffold19"),]
south_tajd.1_chr3=south_tajd.1_rest[which(south_tajd.1_rest$CHROM=="Scaffold61" & south_tajd.1_rest$BIN_START>chrom3_south_start),]
south_tajd.1_rest=south_tajd.1_rest[-which(south_tajd.1_rest$CHROM=="Scaffold61" & south_tajd.1_rest$BIN_START>chrom3_south_start),]

south_tajd.1_chr5=south_tajd.1_rest[which(south_tajd.1_rest$CHROM=="Scaffold0" & south_tajd.1_rest$BIN_START<chrom5_south_end & south_tajd.1_rest$BIN_START>chrom5_south_start),]
south_tajd.1_rest=south_tajd.1_rest[-which(south_tajd.1_rest$CHROM=="Scaffold0" & south_tajd.1_rest$BIN_START<chrom5_south_end & south_tajd.1_rest$BIN_START>chrom5_south_start),]

south_tajd.1_chr1=cbind(rep("chr1",nrow(south_tajd.1_chr1)),south_tajd.1_chr1)
names(south_tajd.1_chr1)[1]="region"
south_tajd.1_chr3=cbind(rep("chr3",nrow(south_tajd.1_chr3)),south_tajd.1_chr3)
names(south_tajd.1_chr3)[1]="region"
south_tajd.1_chr5=cbind(rep("chr5",nrow(south_tajd.1_chr5)),south_tajd.1_chr5)
names(south_tajd.1_chr5)[1]="region"
south_tajd.1_rest=cbind(rep("rest_genome",nrow(south_tajd.1_rest)),south_tajd.1_rest)
names(south_tajd.1_rest)[1]="region"
south_tajd.1_comb=rbind(south_tajd.1_chr1,south_tajd.1_chr3,south_tajd.1_chr5,south_tajd.1_rest)
```



Northern assembly
```
chrom3_north_start=55.8E6

#North pi
north_pi.2_chr1=north_pi.2[which(north_pi.2$CHROM=="Scaffold156"),]
north_pi.2_rest=north_pi.2[-which(north_pi.2$CHROM=="Scaffold156"),]
north_pi.2_chr3=north_pi.2_rest[which(north_pi.2_rest$CHROM=="Scaffold29" & north_pi.2_rest$BIN_START>chrom3_north_start),]
north_pi.2_rest=north_pi.2_rest[-which(north_pi.2_rest$CHROM=="Scaffold29" & north_pi.2_rest$BIN_START>chrom3_north_start),]


north_pi.2_chr5=north_pi.2_rest[which(north_pi.2_rest$CHROM=="Scaffold68"),]
north_pi.2_rest=north_pi.2_rest[-which(north_pi.2_rest$CHROM=="Scaffold68"),]

north_pi.2_chr1=cbind(rep("chr1",nrow(north_pi.2_chr1)),north_pi.2_chr1)
names(north_pi.2_chr1)[1]="region"
north_pi.2_chr3=cbind(rep("chr3",nrow(north_pi.2_chr3)),north_pi.2_chr3)
names(north_pi.2_chr3)[1]="region"
north_pi.2_chr5=cbind(rep("chr5",nrow(north_pi.2_chr5)),north_pi.2_chr5)
names(north_pi.2_chr5)[1]="region"
north_pi.2_rest=cbind(rep("rest_genome",nrow(north_pi.2_rest)),north_pi.2_rest)
names(north_pi.2_rest)[1]="region"
north_pi.2_comb=rbind(north_pi.2_chr1,north_pi.2_chr3,north_pi.2_chr5,north_pi.2_rest)



#South pi
south_pi.2_chr1=south_pi.2[which(south_pi.2$CHROM=="Scaffold156"),]
south_pi.2_rest=south_pi.2[-which(south_pi.2$CHROM=="Scaffold156"),]
south_pi.2_chr3=south_pi.2_rest[which(south_pi.2_rest$CHROM=="Scaffold29" & south_pi.2_rest$BIN_START>chrom3_north_start),]
south_pi.2_rest=south_pi.2_rest[-which(south_pi.2_rest$CHROM=="Scaffold29" & south_pi.2_rest$BIN_START>chrom3_north_start),]

south_pi.2_chr5=south_pi.2_rest[which(south_pi.2_rest$CHROM=="Scaffold68"),]
south_pi.2_rest=south_pi.2_rest[-which(south_pi.2_rest$CHROM=="Scaffold68"),]

south_pi.2_chr1=cbind(rep("chr1",nrow(south_pi.2_chr1)),south_pi.2_chr1)
names(south_pi.2_chr1)[1]="region"
south_pi.2_chr3=cbind(rep("chr3",nrow(south_pi.2_chr3)),south_pi.2_chr3)
names(south_pi.2_chr3)[1]="region"
south_pi.2_chr5=cbind(rep("chr5",nrow(south_pi.2_chr5)),south_pi.2_chr5)
names(south_pi.2_chr5)[1]="region"
south_pi.2_rest=cbind(rep("rest_genome",nrow(south_pi.2_rest)),south_pi.2_rest)
names(south_pi.2_rest)[1]="region"
south_pi.2_comb=rbind(south_pi.2_chr1,south_pi.2_chr3,south_pi.2_chr5,south_pi.2_rest)


#North Tajima's D
north_tajd.2_chr1=north_tajd.2[which(north_tajd.2$CHROM=="Scaffold156"),]
north_tajd.2_rest=north_tajd.2[-which(north_tajd.2$CHROM=="Scaffold156"),]
north_tajd.2_chr3=north_tajd.2_rest[which(north_tajd.2_rest$CHROM=="Scaffold29" & north_tajd.2_rest$BIN_START>chrom3_north_start),]
north_tajd.2_rest=north_tajd.2_rest[-which(north_tajd.2_rest$CHROM=="Scaffold29" & north_tajd.2_rest$BIN_START>chrom3_north_start),]

north_tajd.2_chr5=north_tajd.2_rest[which(north_tajd.2_rest$CHROM=="Scaffold68"),]
north_tajd.2_rest=north_tajd.2_rest[-which(north_tajd.2_rest$CHROM=="Scaffold68"),]

north_tajd.2_chr1=cbind(rep("chr1",nrow(north_tajd.2_chr1)),north_tajd.2_chr1)
names(north_tajd.2_chr1)[1]="region"
north_tajd.2_chr3=cbind(rep("chr3",nrow(north_tajd.2_chr3)),north_tajd.2_chr3)
names(north_tajd.2_chr3)[1]="region"
north_tajd.2_chr5=cbind(rep("chr5",nrow(north_tajd.2_chr5)),north_tajd.2_chr5)
names(north_tajd.2_chr5)[1]="region"
north_tajd.2_rest=cbind(rep("rest_genome",nrow(north_tajd.2_rest)),north_tajd.2_rest)
names(north_tajd.2_rest)[1]="region"
north_tajd.2_comb=rbind(north_tajd.2_chr1,north_tajd.2_chr3,north_tajd.2_chr5,north_tajd.2_rest)

#South Tajima's D
south_tajd.2_chr1=south_tajd.2[which(south_tajd.2$CHROM=="Scaffold156"),]
south_tajd.2_rest=south_tajd.2[-which(south_tajd.2$CHROM=="Scaffold156"),]
south_tajd.2_chr3=south_tajd.2_rest[which(south_tajd.2_rest$CHROM=="Scaffold29" & south_tajd.2_rest$BIN_START>chrom3_north_start),]
south_tajd.2_rest=south_tajd.2_rest[-which(south_tajd.2_rest$CHROM=="Scaffold29" & south_tajd.2_rest$BIN_START>chrom3_north_start),]

south_tajd.2_chr5=south_tajd.2_rest[which(south_tajd.2_rest$CHROM=="Scaffold68"),]
south_tajd.2_rest=south_tajd.2_rest[-which(south_tajd.2_rest$CHROM=="Scaffold68"),]

south_tajd.2_chr1=cbind(rep("chr1",nrow(south_tajd.2_chr1)),south_tajd.2_chr1)
names(south_tajd.2_chr1)[1]="region"
south_tajd.2_chr3=cbind(rep("chr3",nrow(south_tajd.2_chr3)),south_tajd.2_chr3)
names(south_tajd.2_chr3)[1]="region"
south_tajd.2_chr5=cbind(rep("chr5",nrow(south_tajd.2_chr5)),south_tajd.2_chr5)
names(south_tajd.2_chr5)[1]="region"
south_tajd.2_rest=cbind(rep("rest_genome",nrow(south_tajd.2_rest)),south_tajd.2_rest)
names(south_tajd.2_rest)[1]="region"
south_tajd.2_comb=rbind(south_tajd.2_chr1,south_tajd.2_chr3,south_tajd.2_chr5,south_tajd.2_rest)
```

Calculate means for the different chromosome regions, subspecies and assemblies

```
stats=as.data.frame(array(0,c(16,3)))

stats[,1]=c("pi_chr1_south","pi_chr1_north","pi_chr3_south","pi_chr3_north","pi_chr5_south","pi_chr5_north","pi_rest_south","pi_rest_north","tajd_chr1_south","tajd_chr1_north","tajd_chr3_south","tajd_chr3_north","tajd_chr5_south","tajd_chr5_north","tajd_rest_south","tajd_rest_north")

stats[1,2]=mean(south_pi.1_chr1$PI)
stats[2,2]=mean(north_pi.1_chr1$PI)
stats[3,2]=mean(south_pi.1_chr3$PI)
stats[4,2]=mean(north_pi.1_chr3$PI)
stats[5,2]=mean(south_pi.1_chr5$PI)
stats[6,2]=mean(north_pi.1_chr5$PI)
stats[7,2]=mean(south_pi.1_rest$PI)
stats[8,2]=mean(north_pi.1_rest$PI)

stats[1,3]=mean(south_pi.2_chr1$PI)
stats[2,3]=mean(north_pi.2_chr1$PI)
stats[3,3]=mean(south_pi.2_chr3$PI)
stats[4,3]=mean(north_pi.2_chr3$PI)
stats[5,3]=mean(south_pi.2_chr5$PI)
stats[6,3]=mean(north_pi.2_chr5$PI)
stats[7,3]=mean(south_pi.2_rest$PI)
stats[8,3]=mean(north_pi.2_rest$PI)


stats[9,2]=mean(south_tajd.1_chr1$TajimaD,na.rm=T)
stats[10,2]=mean(north_tajd.1_chr1$TajimaD,na.rm=T)
stats[11,2]=mean(south_tajd.1_chr3$TajimaD,na.rm=T)
stats[12,2]=mean(north_tajd.1_chr3$TajimaD,na.rm=T)
stats[13,2]=mean(south_tajd.1_chr5$TajimaD,na.rm=T)
stats[14,2]=mean(north_tajd.1_chr5$TajimaD,na.rm=T)
stats[15,2]=mean(south_tajd.1_rest$TajimaD,na.rm=T)
stats[16,2]=mean(north_tajd.1_rest$TajimaD,na.rm=T)


stats[9,3]=mean(south_tajd.2_chr1$TajimaD,na.rm=T)
stats[10,3]=mean(north_tajd.2_chr1$TajimaD,na.rm=T)
stats[11,3]=mean(south_tajd.2_chr3$TajimaD,na.rm=T)
stats[12,3]=mean(north_tajd.2_chr3$TajimaD,na.rm=T)
stats[13,3]=mean(south_tajd.2_chr5$TajimaD,na.rm=T)
stats[14,3]=mean(north_tajd.2_chr5$TajimaD,na.rm=T)
stats[15,3]=mean(south_tajd.2_rest$TajimaD,na.rm=T)
stats[16,3]=mean(north_tajd.2_rest$TajimaD,na.rm=T)

names(stats)=c("region","southern_genome","northern_genome")

write.table(stats,file="comparison_div_stats_southern_northern_genome.txt",sep="\t",quote=F,col.names=T,row.names=F)
```
