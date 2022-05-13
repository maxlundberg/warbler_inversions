## Call variants

Use `freebayes` to call variants in the four high coverage resequenced samples

```
freebayes-v1.3.1 -f ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.fasta -b 1A05.sorted.nodup.bam -b UK06.sorted.nodup.bam -b 1L19.sorted.dedup.bam -b 1M13.sorted.nodup.bam -b DW83.sorted.nodup.bam -v pacbio.vcf
gzip -c pacbio.vcf > pacbio.vcf.gz
```


## Run gIMble

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

