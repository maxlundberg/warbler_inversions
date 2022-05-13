# Breakpoint analyses


## Whole-genome alignments

Use `mummer4` to align the southern willow warbler scaffolds in or adjacent to the divergent regions to chromosomes 1,3 and 5 
in the flycatcher genome

``` 
chroms="1 3 5"

for chrom in $chroms
do
samtools faidx FicAlb.1.5.genome.ncbi.fasta $chrom > FicAlb.chr${chrom}.fasta
done

scaffolds="Scaffold19 Scaffold11 Scaffold12"

for scaff in $scaffolds
do
samtools faidx ww_southern_hifi_bionano.filt.new.fasta $scaff 
done > chr1.scaffolds.ww_south.fasta

scaffolds="Scaffold61 Scaffold38"

for scaff in $scaffolds
do
samtools faidx ww_southern_hifi_bionano.filt.new.fasta $scaff 
done > chr3.scaffolds.ww_south.fasta

samtools faidx ww_southern_hifi_bionano.filt.new.fasta Scaffold0 > chr5.scaffolds.ww_south.fasta


for chrom in $chroms
do
nucmer --threads=1 FicAlb.chr$chrom.fasta chr$chrom.scaffolds.ww_south.fasta --prefix=chr$chrom.scaffolds.ww_south.vs.flycatcher
delta-filter -1 chr$chrom.scaffolds.ww_south.vs.flycatcher.delta > chr$chrom.scaffolds.ww_south.vs.flycatcher.filt
show-coords chr$chrom.scaffolds.ww_south.vs.flycatcher.vs.flycatcher.delta.filt -dT | sed '1,4d' > chr$chrom.scaffolds.ww_south.vs.flycatcher.delta.filt.coords.out
done
```

Also align the northern willow warbler and chiffchaff scaffolds that contain the divergent region and predicted adjacent scaffolds to the homologous
scaffolds in the southern willow warbler assembly 

```
scaffolds="Scaffold156 Scaffold374 Scaffold143"

for scaff in $scaffolds
do
samtools faidx ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.2.fasta $scaff 
done > chr1.scaffolds.ww_north.fasta

scaffolds="Scaffold29 Scaffold29b Scaffold139"

for scaff in $scaffolds
do
samtools faidx ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.2.fasta $scaff 
done > chr3.scaffolds.ww_north.fasta

samtools faidx ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.2.fasta Scaffold68 > chr5.scaffolds.ww_north.fasta  


scaffolds="ptg000006l ptg000007l ptg000024l"

for scaff in $Scaffolds
do
samtools faidx hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta $scaff
done > chr1.scaffolds.cc.fasta

scaffolds="ptg000026l ptg000040l"

for scaff in $Scaffolds
do
samtools faidx hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta $scaff
done > chr3.scaffolds.cc.fasta

samtools faidx hifiasm_chiffchaff.asm.bp.p_ctg.fasta ptg000051l > mummer_div_regions/chr5.scaffolds.cc.fasta  

chroms="1 3 5"
sample_list="ww_north cc"

for chrom in $chroms
do
for sample in $sample_list
do
nucmer --threads=2 chr$chrom.scaffolds.ww_south.fasta chr$chrom.scaffolds.$sample.fasta --prefix=chr$chrom.scaffolds.$sample.vs.southern_ww
delta-filter -1 chr$chrom.scaffolds.$sample.vs.southern_ww.delta > chr$chrom.scaffolds.$sample.vs.southern_ww.filt.2.delta
show-coords -rdHT chr$chrom.scaffolds.$sample.vs.southern_ww.filt.2.delta > chr$chrom.scaffolds.$sample.vs.southern_ww.filt.coords.out
done
done
```


## Optical maps


Map the southern willow warbler optical map to the northern willow warbler assembly

```
perl Solve3.2.2_08222018/Pipeline/08222018/fa2cmap_multi_color.pl -i ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -e BspQI 1 -o ww_pacbio_bionano

python Solve3.2.2_08222018/Pipeline/08222018/runCharacterize.py -t Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -r ww_pacbio_bionano/ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt_BSPQI_0kb_0labels.cmap -q exp_refineFinal1_contigs.cmap -n19 -a optArguments_human.xml

python Solve3.2.2_08222018/Pipeline/08222018/runSV.py -t Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -r ww_pacbio_bionano/ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt_BSPQI_0kb_0labels.cmap -q exp_refineFinal1_contigs.cmap -T 19 -j 19 -a optArguments_human.xml -e alignref/exp_refineFinal1_contigs.err -o southern_bionano_vs_ref_genome_sv
 
python smap_to_vcf_v2.py -s exp_refineFinal1_contigs.smap -o southern_bionano_vs_northern_ref -b False
cat southern_vs_reference_bionano_sv.vcf | change_name_sequences_bionano_sv_vcf.pl --key_file ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt_BSPQI_0kb_0labels_key.txt > southern_vs_reference_bionano_sv.new.vcf 
```

Map the northern willow warbler optical map to the southern willow warbler assembly

```
perl Solve3.2.2_08222018/Pipeline/08222018/fa2cmap_multi_color.pl -i ww_southern_hifi_bionano.filt.new.fasta -e BspQI 1 -o ww_southern_hifi_bionano.filt.new

python Solve3.2.2_08222018/Pipeline/08222018/runCharacterize.py -t Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -r ww_southern_hifi_bionano.filt.new/ww_southern_hifi_bionano.filt.new_BSPQI_0kb_0labels.cmap -q exp_refineFinal1_contigs.cmap -n19 -a optArguments_human.xml

python Solve3.2.2_08222018/Pipeline/08222018/runSV.py -t Solve3.2.2_08222018/RefAligner/7782.7865rel/RefAligner -r ww_southern_hifi_bionano.filt.new/ww_southern_hifi_bionano.filt.new_BSPQI_0kb_0labels.cmap -q exp_refineFinal1_contigs.cmap -T 19 -j 19 -a optArguments_human.xml -e alignref/exp_refineFinal1_contigs.err -o northern_bionano_vs_southern_hifi_genome_sv

python ~/smap2vcf/smap_to_vcf_v2.py -s exp_refineFinal1_contigs.smap -o northern_bionano_vs_southern_hifi_sv -b False
cat northern_bionano_vs_southern_hifi_sv.vcf | ~/change_name_sequences_bionano_sv_vcf.pl --key_file ../ww_southern_hifi_bionano.filt.new/ww_southern_hifi_bionano.filt.new_BSPQI_0kb_0labels_key.txt > northern_bionano_vs_southern_hifi_sv.new.vcf                
```


Also hybridize the northern willow warbler optical map to the southern willow warbler assembly in order to visualize the differences close to breakpoints. This time don't perform any cuts in the reference
assembly

```
genome="ww_southern_hifi_bionano.filt.new.fasta"

perl Solve3.6.1_11162020/HybridScaffold/11162020/hybridScaffold.pl -f -B 2 -N 1 -r Solve3.6.1_11162020/RefAligner/11442.11643rel/RefAligner -n $genome -b exp_refineFinal1_contigs.cmap -c hybridScaffold_config_solve_aggressive.xml -o $genome.bionano.northern.solve_config.aggressive.no_cut.output

```


## Map 10x chromium data to detect structural variants

Here we will use the `longranger wgs` pipeline to check for structural variants in the divergent regions, in particular inversion breakpoints

Map the linked reads from the three willow warbler samples to the northern willow warbler assembly

```
genome="ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta"
genome_prefix=$(echo $genome | sed 's/.fasta//')

longranger mkref $genome

#Set sex to female as it has two copies of the sex chromosome (Z), as would be the case in the female of mammals
longranger wgs --reference=refdata-$genome_prefix --id=P6352_101_uppmax_longranger_wgs --fastqs=./ --sex=female
longranger wgs --reference=refdata-$genome_prefix --id=P6352_101_uppmax_longranger_wgs_southern_genome --fastqs=./ --sex=female
longranger wgs --reference=refdata-$genome_prefix --id=P13103_302_longranger_wgs --fastqs=./ --sex=female
```

Map the linked reads from the three willow warbler samples to the southern willow warbler assembly. In this case we will extract the 500 largest scaffolds from the assembly

```
samtools faidx ww_southern_hifi_bionano.filt.new.fasta

cat ww_southern_hifi_bionano.filt.new.fasta.fai | sort -k2,2nr | head -n498 | cut -f1 > ww_southern_hifi_bionano.filt.new.fasta.largest_499.scaffold.ids.list
echo "mito" >>  ww_southern_hifi_bionano.filt.new.fasta.largest_499.scaffold.ids.list

cat ww_southern_hifi_bionano.filt.new.fasta.fai | sort -k2,2nr | tail -n+499 | grep -v "mito" | cut -f1 > ww_southern_hifi_bionano.filt.new.fasta.rest.scaffolds.list

echo ">Rem_scaffolds" > ww_southern_hifi_bionano.filt.new.rest.scaffolds.fasta

cat ww_southern_hifi_bionano.filt.new.fasta.rest.scaffolds.list | while read scaffold
do
samtools faidx ww_southern_hifi_bionano.filt.new.fasta $scaffold | perl -ne 'if($_!~/>/){chomp($_);print $_}'    
perl -e 'print "N"x500'      
done >> ww_southern_hifi_bionano.filt.new.rest.scaffolds.fasta

faSomeRecords ww_southern_hifi_bionano.filt.new.fasta ww_southern_hifi_bionano.filt.new.fasta.largest_499.scaffold.ids.list ww_southern_hifi_bionano.filt.new.fasta.largest_499.scaffold.fasta

genome="ww_southern_hifi_bionano.filt.new.500.scaffolds.fasta"
genome_prefix=$(echo $genome | sed 's/.fasta//')

longranger mkref $genome

longranger wgs --reference=refdata-$genome_prefix --id=P6352_101_longranger_wgs_southern_hifi_genome --fastqs=./ --sex=female
longranger wgs --reference=refdata-$genome_prefix --id=P6352_102_longranger_wgs_southern_hifi_genome --fastqs=./ --sex=female
longranger wgs --reference=refdata-$genome_prefix --id=P13103_302_longranger_wgs_southern_hifi_genome --fastqs=./ --sex=female
```


## Calculate coverage of linked reads 

Map linked reads of the three warbler samples to the southern genome

```
bwa mem ww_southern_hifi_bionano.filt.fasta -p northern_barcoded.fastq.gz -t 19 -C | samtools view -bS - > P6352_101.southern_genome.bam
samtools sort -@ 20 P6352_101.southern_genome.bam > P6352_101.southern_genome.sorted.bam
samtools index P6352_101.southern_genome.sorted.bam

bwa mem ww_southern_hifi_bionano.filt.fasta -p southern_barcoded.fastq.gz -t 19 -C | samtools view -bS - > P6352_102.southern_genome.bam
samtools sort -@ 20 P6352_102.southern_genome.bam > P6352_102.southern_genome.sorted.bam
samtools index P6352_102.southern_genome.sorted.bam

bwa mem ww_southern_hifi_bionano.filt.fasta -p southern.302_barcoded.fastq.gz -t 19 -C | samtools view -bS - > P13103_302.southern_genome.bam
samtools sort -@ 20 P13103_302.southern_genome.bam > P13103_302.southern_genome.sorted.bam
samtools index P13103_302.southern_genome.sorted.bam
```


Extract alignments from the scaffolds in the divergent region on chromosome 5 (Scaffold68 in the northern genome and Scaffold36 in the southern genome)
and in the divergent region on chromosome 1 in the southern genome (Scaffold11).

```
#Northern genome
samtools view -h P6352_101.sorted.bam Scaffold68 | samtools sort -tBX - > P6352_101.Scaffold68.bc.bam
samtools view -h P6352_102.sorted.bam Scaffold68 | samtools sort -tBX - > P6352_102.Scaffold68.bc.bam

#Southern genome 
samtools view -h P6352_101.southern_genome.sorted.bam Scaffold11 | samtools sort -tBX - > P6352_101.southern_genome.Scaffold11.bc.bam
samtools view -h P6352_102.southern_genome.sorted.bam Scaffold11 | samtools sort -tBX - > P6352_102.southern_genome.Scaffold11.bc.bam

samtools view -h P6352_101.southern_genome.sorted.bam Scaffold36 | samtools sort -tBX - > P6352_101.southern_genome.Scaffold36.bc.bam
samtools view -h P6352_102.southern_genome.sorted.bam Scaffold36 | samtools sort -tBX - > P6352_102.southern_genome.Scaffold36.bc.bam
```

Extract molecule positions (bed file) for each sample and region using `tigmint`

```
bam_files=$(ls *.bam)

molecule_size=10000

for file in $bam_files
do
samtools sort -tBX $file > $file.bc.bam
tigmint-molecule $file.bc.bam -s $molecule_size -q 1  | sort -k1,1 -k2,2n -k3,3n > $file.molecules.sorted.bed
done
``` 


Quantify the number of molecules overlapping non-overlapping 1 kb windows using `BEDTools`

```
prefixes=$(ls *.molecules.bed | sed 's/.bed//g')

for prefix in $prefixes
do
bedtools coverage -a ../ww_southern_hifi_bionano.filt.fasta.1kb_windows.sorted.bed -b $prefix.sorted.bed > $prefix.1kb_coverage.bed
done

```


## Map long reads 

Map long reads to the southern assembly and the chiffchaff assembly using `minimap2`. These alignments will be used to visualize breakpoints at the 
start of the chromosome 3 divergent region. The alignment of the reads of the northern willow warbler to the southern assembly will also be used
to identify high-confidence SVs to be genotyped with delly. 

```
#Southern willow warbler assembly 
minimap2 -a -t 19 --MD -x map-hifi ww_southern_hifi_bionano.filt.new.fasta ww_southern_hifi_reads.fastq.gz  | samtools view -bS - > ww_southern_hifi_bionano.filt.fasta.ww_southern.hifi_reads.bam
samtools sort -@ 20 ww_southern_hifi_bionano.filt.fasta.ww_southern.hifi_reads.bam > ww_southern_hifi_bionano.filt.fasta.ww_southern.hifi_reads.sorted.bam
samtools index ww_southern_hifi_bionano.filt.fasta.ww_southern.hifi_reads.sorted.bam

minimap2 -a -t 19 --MD -x map-pb ww_southern_hifi_bionano.filt.new.fasta ps_036_subreads.fastq.gz | samtools view -bS - > ww_southern_hifi_bionano.filt.fasta.ww_northern.pacbio_reads.bam
samtools sort -@ 20 ww_southern_hifi_bionano.filt.fasta.ww_northern.pacbio_reads.bam > ww_southern_hifi_bionano.filt.fasta.ww_northern.pacbio_reads.sorted.bam
samtools index ww_southern_hifi_bionano.filt.fasta.ww_northern.pacbio_reads.sorted.bam

minimap2 -a -t 19 --MD -x map-hifi ww_southern_hifi_bionano.filt.new.fasta chiffchaff_hifi_reads.fastq.gz  | samtools view -bS - > ww_southern_hifi_bionano.filt.fasta.chiffchaff.hifi_reads.bam
samtools sort -@ 20 ww_southern_hifi_bionano.filt.fasta.chiffchaff.hifi_reads.bam > ww_southern_hifi_bionano.filt.fasta.chiffchaff.hifi_reads.sorted.bam
samtools index ww_southern_hifi_bionano.filt.fasta.chiffchaff.hifi_reads.sorted.bam

#Chiffchaff assembly
minimap2 -a -t 19 --MD -x map-hifi hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta ww_southern_hifi_reads.fastq.gz  | samtools view -bS - > hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta.ww_southern.hifi_reads.bam
samtools sort -@ 20 hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta.ww_southern.hifi_reads.bam > hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta.ww_southern.hifi_reads.sorted.bam
samtools index hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta.ww_southern.hifi_reads.sorted.bam

minimap2 -a -t 19 --MD -x map-pb hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta ps_036_subreads.fastq.gz | samtools view -bS - > hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta.ww_northern.pacbio_reads.bam
samtools sort -@ 20 hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta.ww_northern.pacbio_reads.bam > hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta.ww_northern.pacbio_reads.sorted.bam
samtools index hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta.ww_northern.pacbio_reads.sorted.bam

minimap2 -a -t 19 --MD -x map-hifi hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta chiffchaff_hifi_reads.fastq.gz  | samtools view -bS - > hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta.chiffchaff.hifi_reads.bam
samtools sort -@ 20 hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta.chiffchaff.hifi_reads.bam > hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta.chiffchaff.hifi_reads.sorted.bam
samtools index hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta.chiffchaff.hifi_reads.sorted.bam
```


## Repeatmask the northern willow warbler and chiffchaff assembly  

The northern willow warbler assembly is annotated based on the same repeats as identified in the southern assembly

```
RepeatMasker -pa 20 -s -lib combined_aves_rmodeler.lib.fa -dir ./ ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.2.fasta
```

For the chiffchaff, we will first create a species-specific *de novo* library with `repeatmodeler` and then combine this with the focal
chr1-3 tandem repeat-associated repeat (rnd-5_family-604) from the southern willow warbler and repeats from other birds

```
BuildDatabase -name chiffchaff_hifiasm ../hifiasm_chiffchaff.asm.bp.p_ctg.fasta
RepeatModeler -engine ncbi -pa 20 -database chiffchaff_hifiasm

cat RM_20643.WedSep221803572021/consensi.fa.classified aves_repeats.fasta rnd-5_family-604.fasta > combined_aves_rmodeler.with_rnd-5-604.lib.fa 

RepeatMasker -pa 20 -s -lib combined_aves_rmodeler.with_rnd-5-604.lib.fa -dir ./ hifiasm_chiffchaff.asm.bp.p_ctg.new.fasta

```




