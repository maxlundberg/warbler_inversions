# Annotation of the southern assembly


## Repeatmask genome

First use `repeatmodeler` for *de novo* identification of repeats

```
BuildDatabase -name repeatmodeler_ww_southern_hifi_bionano.filt.new ww_southern_hifi_bionano.filt.new.fasta
RepeatModeler -engine ncbi -pa 20 -database ww_southern_hifi_bionano.filt
```

Next use `Repeatmasker` to mask the genome. For this part we will download the repbase repeat release from 20181026 and the Dfam consensus 20171107 library.
Combine this with the output from repeatmodeler (consensi.fa.classified).

```
RepeatMasker/util/queryRepeatDatabase.pl -species aves > aves_repeats.fasta 

cat aves_repeats.fasta consensi.fa.classified > combined_aves_rmodeler.lib.fa 

RepeatMasker -pa 20 -s -lib combined_aves_rmodeler.lib.fa -dir ./ ww_southern_hifi_bionano.filt.new.fasta

cat ww_southern_hifi_bionano.filt.new.fasta.out | perl -lane 'next if($.<4);$,="\t";print $F[4],$F[5]-1,$F[6],$F[9]' | sort -k1,1 -k2,2n > ww_southern_hifi_bionano.filt.new.fasta.bed
```

Softmask the genome with `BEDTools`

```
bedtools maskfasta -fi ww_southern_hifi_bionano.filt.new.fasta -bed ww_southern_hifi_bionano.filt.new.fasta.bed -soft -fo ww_southern_hifi_bionano.filt.new.sm.fasta
```


## Map protein data

Download protein data from NCBI for chicken, zebra finch and great tit, as well as bird proteins from uniprot with transcript evidence, protein evidence or that are manually curated. Divide
this dataset into batches of 100 sequences using `fasta-splitter` (http://kirill-kryukov.com/study/tools/fasta-splitter/) to facilitate parallelization


```
prot_files="GGa6a.proteins.fasta Parus_major1.1.proteins.fasta TaeGut1.4.protein.fasta swissprot_aves_20211025.fasta uniprot_aves_protein_evidence_20211025.fasta uniprot_aves_transcript_evidence_20211025.fasta"

mkdir prot_parts
cd prot_parts

for file in $prot_files
do
fasta-splitter.pl --part-size 100 --measure count ../$file
done
```

Map batches of proteins to the softmasked assembly using `exonerate` and speed up the process using `gnu-parallel`

```
ls prot_parts/ |  sed 's/.fasta//' > protein_parts.list
parallel -j 20 'exonerate -t ww_southern_hifi_bionano.filt.new.sm.fasta  -q ./protein_parts/{}.fasta --model protein2genome --bestn 1 --showtargetgff --fsmmemory 5000 --showvulgar no --showalignment no --showquerygff no --ryo "AveragePercentIdentity: %pi\n" --softmasktarget > ./exonerate_gff/{}.gff' :::: protein_parts.list 

cat exonerate_gff/*.gff | egrep -v "Hostname|Average" > exonerate_for_augustus.gff
```

Extract hints for augustus using scripts from `BRAKER`

```
align2hints.pl  --in=exonerate_for_augustus.gff  --out=exonerate_for_augustus.hints.new.gff --prg=exonerate --genome ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm.fasta

cat exonerate_for_augustus.hints.new.gff | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1 | join_mult_hints.pl > exonerate_hints.nr.new.gff
```


## RNAseq data

Use `STAR` to map the RNAseq data to the softmasked assembly. This data has previously been trimmed using `trimgalore` with default settings

```
mkdir ww_southern_hifi_bionano.filt.new.sm.fasta.star
star --runThreadN 8 --runMode genomeGenerate --genomeDir ./ww_southern_hifi_bionano.filt.new.sm.fasta.star --genomeFastaFiles ww_southern_hifi_bionano.filt.new.sm.fasta
star --runMode alignReads --genomeDir ./ww_southern_hifi_bionano.filt.new.sm.fasta.star --runThreadN 18 --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outSAMattributes Standard --outFileNamePrefix merged_ww. --readFilesIn 8_160311_BC8B3GANXX_P3956_112_1_val_1.fq.gz,8_160311_BC8B3GANXX_P3956_128_1_val_1.fq.gz,8_160311_BC8B3GANXX_P3956_144_1_val_1.fq.gz,8_160311_BC8B3GANXX_P3956_160_1_val_1.fq.gz,8_160311_BC8B3GANXX_P3956_176_1_val_1.fq.gz,8_160311_BC8B3GANXX_P3956_192_1_val_1.fq.gz 8_160311_BC8B3GANXX_P3956_112_2_val_2.fq.gz,8_160311_BC8B3GANXX_P3956_128_2_val_2.fq.gz,8_160311_BC8B3GANXX_P3956_144_2_val_2.fq.gz,8_160311_BC8B3GANXX_P3956_160_2_val_2.fq.gz,8_160311_BC8B3GANXX_P3956_176_2_val_2.fq.gz,8_160311_BC8B3GANXX_P3956_192_2_val_2.fq.gz --readFilesCommand zcat
samtools index merged_ww.Aligned.sortedByCoord.out.bam
```


Create hints for augustus using auxiliary augustus scripts, `samtools` and `strand_cov` (https://github.com/pmenzel/stranded-coverage)

```
samtools sort -@20 -n merged_ww.Aligned.sortedByCoord.out.bam > merged_ww_rnaseq.ns.bam
filterBam --uniq --paired --pairwiseAlignment --in merged_ww_rnaseq.ns.bam --out rnaseq.ns.filt.bam
samtools sort rnaseq.ns.filt.bam  > rnaseq.filt.sorted.bam
bam2hints --intronsonly --in=rnaseq.filt.sorted.bam --out=intron_hints_rnaseq.gff

strand_cov -o rnaseq.filt.sorted rnaseq.filt.sorted.bam

cat rnaseq.filt.sorted.plus.wig | wig2hints.pl width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep  --UCSC=stranded.track --radius=4.5 --pri=4 --strand="+" > exonpart_hints_rnaseq_plus.gff
cat rnaseq.filt.sorted.minus.wig | wig2hints.pl width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep  --UCSC=stranded.track --radius=4.5 --pri=4 --strand="-" > exonpart_hints_rnaseq_minus.gff
cat exonpart_hints_rnaseq_plus.gff exonpart_hints_rnaseq_minus.gff > exonpart_hints_rnaseq_both_strands.gff
```

Combine RNAseq hints and protein hints and divide them into separate files for each scaffold to make gene prediction more efficient

```
cat intron_hints_rnaseq.gff exonpart_hints_rnaseq_both_strands.gff exonerate_hints.nr.new.gff | sort -k1,1 -k4,4n> rnaseq_exonerate_protein.hints.new.gff


mkdir augustus_rnaseq_protein_hints_new

cat scaffolds.list | while read scaffold
do
touch augustus_rnaseq_protein_hints_new/$scaffold.hints.gff
done

cat rnaseq_exonerate_protein.hints.new.gff | perl -lane '$,="\t"; open(OUT,">>augustus_rnaseq_protein_hints_new/$F[0].hints.gff");print OUT @F; close(OUT)'
```


## Generate gene models 

Use `augustus` to predict genes. For the prediction we will use specific-specific parameters that originate from training based on transcript data mapped
to a short-read genome. For details about the training see the augustus training workflow in the scripts folder. 

```
mkdir scaffold_fasta

cat scaffolds.list | while read scaffold
do
samtools faidx ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm.fasta $scaffold > scaffold_fasta/$scaffold.fasta
done

mkdir augustus_rnaseq_intron_exonpart_exonerate_new_gffs
parallel -j 20 "augustus --alternatives-from-evidence=true --softmasking=true --UTR=on --gff3=on --hintsfile=augustus_rnaseq_protein_hints_new/{}.hints.gff --species=ww_old_genome_webapollo_1200 --allow_hinted_splicesites=atac --extrinsicCfgFile=augustus_config/extrinsic/extrinsic.M.RM.E.W.P.own.12.cfg scaffold_fasta/{}.fasta > augustus_rnaseq_intron_exonpart_exonerate_new_gffs/{}.gff3" :::: scaffolds.list

cat augustus_rnaseq_intron_exonpart_exonerate_new_gffs/*gff3 | join_aug_pred.pl > augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.new.gff3
```


## Assemble transcripts

Also assemble transcripts from the RNAseq data using `Hisat2` and `Stringtie`

```
hisat2-build ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm.fasta ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.sm.hisat2

prefixes=$(ls *.gz | cut -b1-29 | sort | uniq)

for prefix in $prefixes
do
hisat2 -x $genome_database -1 ${prefix}_1_val_1.fq.gz -2 ${prefix}_2_val_2.fq.gz -p 19 --rna-strandness RF --dta | samtools view -bS - > $prefix.bam
samtools sort -@ 19 $prefix.bam > $prefix.sorted.bam
samtools index $prefix.sorted.bam
done

ls *.sorted.bam > bamfiles.list
samtools merge -@ 19 -b bamfiles.list merged_ww_rnaseq.hisat2.new.bam
samtools index merged_ww_rnaseq.hisat2.new.bam

stringtie --rf -p 4 merged_ww_rnaseq.hisat2.new.bam > stringtie.new.gtf
gffread -E stringtie.new.gtf -o- > stringtie.new.gff3
gffread -w transcripts.stringtie.new.fasta -g ww_pacbio_bionano_arcs_pbjelly_pilon.filt_with_mt.filt.sm.fasta stringtie.new.gtf
```


## Synteny-transferred models

Transfer chicken annotations to the willow warbler genome by first aligning the warbler and zebra finch (bTaeGut1.4.pri) with `satsuma` and then using `kraken`

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/565/GCF_003957565.2_bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/565/GCF_003957565.2_bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_genomic.gtf.gz

zcat GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz | perl -lane 'if($F[0]=~/(>\S+)/){print $1}else{print @F}' > bTaeGut1.4.pri_genomic.2.fasta

SatsumaSynteny -q ww_southern_hifi_bionano.filt.new.fasta -t bTaeGut1.4.pri_genomic.2.fasta -n 19 -o ww_southern_hifi_bionano.filt.new.fasta.zf_new.2.satsuma > ww_southern_hifi_bionano.filt.new.fasta.zf_new.2.satsuma.stdout

RunKraken -c kraken_zf_ncbi_ww_south_hifi.config -s GCF_003957565.2_bTaeGut1.4.pri_genomic.gtf  -S zf -T warbler -o kraken_zf_ww_hifi.gtf
CleanKrakenFiles -i kraken_zf_ww_hifi.gtf > kraken_zf_ww_hifi.filt.gtf
```

Intersect CDS annotation of the synteny-transferred genes with the CDS annotation of the augustus predictions using `bedtools`

```
cat kraken_zf_ww_hifi.filt.gtf | awk '{if($3=="CDS"){print $0}}' > kraken_zf_ww_hifi.filt.cds.gtf

cat augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.gff3 | awk '{if($3=="CDS"){print $0}}' > augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.cds.gff3

bedtools intersect -wao -s -a augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.cds.gff3 -b kraken_zf_ww_hifi.filt.cds.gtf > augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.cds.zf_kraken_overlap.out
```


Summarize the intersection

```
~/summarize_gene_models_synteny_intersect.pl --intersection_file augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.cds.zf_kraken_overlap.out > summary_overlap_augustus_kraken.out
```



## Blast augustus gene models to bird proteins

Extract amino acid sequence for the longest isoform using scripts from `augustus` and customized scripts

```
getAnnoFasta.pl augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.new.gff3
mv augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.new3.aa augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.fasta
cat augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.fasta | extract_longest_aug_transcript.pl > augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.proteins.longest_transcript.fasta 
```

Extract the longest protein for each gene from chicken, zebra finch and great tit downloaded from NCBI

```
#Extract ID of longest transcripts from the feature tables
cat  GCF_003957565.2_bTaeGut1.4.pri_feature_table.txt | awk '{if($1=="mRNA"){print $0}}' | cut -f13,15,19 | sort -k2,2 -k3,3nr | sort -u -k2,2 > zf_proteins.longest_transcript_info.txt
zcat GCF_001522545.3_Parus_major1.1_feature_table.txt.gz | awk '{if($1=="mRNA"){print $0}}' | cut -f13,15,19 | sort -k2,2 -k3,3nr | sort -u -k2,2 > gt_proteins.longest_transcript_info.txt
zcat GCF_000002315.5_GRCg6a_feature_table.txt.gz |awk '{if($1=="mRNA"){print $0}}' | cut -f13,15,19 | sort -k2,2 -k3,3nr | sort -u -k2,2  > chicken_proteins.longest_transcript_info.txt

#Modify the headers of the amino acid fasta files and index them
files="GCF_001522545.3_Parus_major1.1_protein GCF_003957565.2_bTaeGut1.4.pri_protein GCF_000002315.6.chicken.protein.faa"
for file in $files
do
cat $file.faa | perl -lane 'if($F[0]=~/(>\S+)/){print $1}else{print @F}' > $file.filt.fasta
samtools faidx $file.filt.fasta
done

cat gt_proteins.longest_transcript_info.txt | cut -f1-2 | while read transcript gene
do
export gene transcript
samtools faidx GCF_001522545.3_Parus_major1.1_protein.filt.fasta $transcript | perl -lane 'if($F[0]=~/(>\S+)/){print ">",$ENV{gene},":",$ENV{transcript}}else{print @F}'
done > gt_proteins.longest_transcript.fasta

cat zf_proteins.longest_transcript_info.txt | cut -f1-2 | while read transcript gene
do
export gene transcript
samtools faidx GCF_003957565.2_bTaeGut1.4.pri_protein.filt.fasta $transcript | perl -lane 'if($F[0]=~/(>\S+)/){print ">",$ENV{gene},":",$ENV{transcript}}else{print @F}'
done > zf_proteins.longest_transcript.fasta

cat chicken_proteins.longest_transcript_info.txt | egrep "NP|XP" | cut -f1-2 | while read transcript gene
do
export gene transcript
samtools faidx GCF_000002315.6.chicken.protein.filt.fasta $transcript | perl -lane 'if($F[0]=~/(>\S+)/){print ">",$ENV{gene},":",$ENV{transcript}}else{print @F}'
done > chicken_proteins.longest_transcript.fasta

cat chicken_proteins.longest_transcript.fasta gt_proteins.longest_transcript.fasta zf_proteins.longest_transcript.fasta > bird_proteins.combined.fasta
```


Perform a reciprocal blast search between the willow warbler proteins and the bird proteins

```
makeblastdb -in bird_proteins.combined.fasta -out bird_proteins.combined.fasta -title bird_proteins.combined.fasta -dbtype prot
blastp -num_threads 10 -query augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.fasta -db bird_proteins.combined.fasta -outfmt 6 -evalue 1e-5 > augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.bird_proteins.blast.out 

makeblastdb -in augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.fasta -out augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.fasta -title augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.fasta -dbtype prot
blastp -num_threads 20 -query bird_proteins.combined.fasta -db augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.fasta -outfmt 6 -evalue 1e-5 > bird_proteins.vs.ww_south_proteins.blast.out

cat augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.bird_proteins.blast.out |  sort -u -k1,1 > augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.bird_proteins.blast.filt.out     
cat bird_proteins.vs.ww_south_proteins.blast.out | sort -u -k1,1 > bird_proteins.vs.ww_south_proteins.blast.filt.out 

#Summarize the output
~/parse_reciprocal_blast_data.new.pl --query_vs_target augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.bird_proteins.blast.filt.out --target_vs_query bird_proteins.vs.ww_south_proteins.blast.filt.out > reciprocal_blast_summary.out
```


Also blast against vertebrate swissprot proteins (N=86,131)

```
makeblastdb -in swissprot_vertebrata_20211118.prot.fasta -out swissprot_vertebrata_20211118.prot.fasta -title swissprot_vertebrata_20211118.prot.fasta -dbtype prot
blastp -num_threads 20 -query augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.fasta -db swissprot_vertebrata_20211118.prot.fasta -evalue 1e-5 -outfmt 6 > augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.fasta.blast.swissprot.vertebrate.out
cat augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.fasta.blast.swissprot.vertebrate.out  | sort -u -k1,1 > augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.fasta.blast.swissprot.vertebrate.filt.out   
```


## Identify protein domains


Identify protein domains using `Interproscan`. First split the data into 20 batches using `fasta-splitter`

```
interproscan.sh -i augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.fasta  -dp -pa -appl Pfam,PANTHER,SMART,CDD --goterms --iprlookup

summarize_interpro_output.pl --interpro_output_file augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.fasta.tsv > interpro_summary_augustus.ww_southern.hifi.txt                          
```


## Integrate annotation data and filter gene models


Collect general info from the augustus predictions

```
cat augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.corrected.new.gff3 | collect_info_from_augustus_gff.pl > summary_augustus.out
```

Extract information about the bird proteins used in the blast search

```
cat GCF_003957565.2_bTaeGut1.4.pri_feature_table.txt | awk '{if($1=="mRNA"){print $0}}' | cut -f13,14,15,19 | sort -k3,3 -k4,4nr  -t $'\t' | sort -u -k3,3 -t $'\t' | while read line
do
echo -e "Zebra finch\t$line"
done > longest_proteins.ncbi.birds.info.txt

zcat GCF_001522545.3_Parus_major1.1_feature_table.txt.gz | awk '{if($1=="mRNA"){print $0}}' | cut -f13,14,15,19 | sort -k3,3 -k4,4nr  -t $'\t' | sort -u -k3,3 -t $'\t' | while read line
do
echo -e "Great tit\t$line"
done >> longest_proteins.ncbi.birds.info.txt

zcat GCF_000002315.5_GRCg6a_feature_table.txt.gz | awk '{if($1=="mRNA"){print $0}}' | cut -f13,14,15,19 | sort -k3,3 -k4,4nr  -t $'\t' | sort -u -k3,3 -t $'\t' | while read line
do
echo -e "Chicken\t$line"
done >> longest_proteins.ncbi.birds.info.txt
```

Summarize all the information 

```
summarize_annotation_data.new.pl --gene_info_file summary_augustus.out --blast_summary_file reciprocal_blast_summary.out --blast_target_info_file longest_proteins.ncbi.birds.info.txt --kraken_intersect_file summary_overlap_augustus_kraken.out --swissprot_blast augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.protein.longest.fasta.blast.swissprot.vertebrate.filt.out  --swissprot_info_file swissprot_vertebrata_20211118.prot.info.txt --interpro_summary interpro_summary_augustus.ww_southern.hifi.txt > augustus_annotation_combined.info.ww_southern.hifi.txt
```

Extract a list of old augustus names and new names after filtering genes containing no protein domain and use this list to
create a new filtered gff3 file with renamed gene models

```
change_gene_names_aug_models_gff3.pl --aug_gff augustus_pacbio_bionano_rnaseq_intron_exonpart_exonerate.combined.gff3 --name_list new_names_filt_augustus.txt > ww_southern_hifi_augustus_new_names.filt.gff3
```

