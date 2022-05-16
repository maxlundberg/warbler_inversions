#This workflow describes the steps that were used to create transcripts used for training augustus

#For this purpose we used the trimmed RNAseq data published in this study and a previous short-read assembly from a northern willow warbler


##########################################

#de novo assembly of RNAseq data

##########################################


#First concatenate the RNAseq reads from each sample

zcat *1.fq.gz | gzip -c > ww_rnaseq_trimmed_all_1.fq.gz
zcat *2.fq.gz | gzip -c > ww_rnaseq_trimmed_all_2.fq.gz


#Use Trinity (2.0.2) to create a de novo transcriptome assembly

Trinity --seqType fq --max-memory 250G --left combined_rnaseq_1.fq.gz --right combined_rnaseq_1.fq.gz --CPU 20 --SS_lib_type RF --output ww_trinity




##########################################

#genome-guided assembly of RNAseq data

##########################################

#In addition to the de novo assembly, we will perform a genome-guided assembly with trinity

#First map the reads to the genome using gsnap using a previous database (gmap_build -d $database $genome.fasta)

database="ww_genome_final_masked"
threads=15
npaths=1

samples=$(ls *fq.gz | sed 's/_[12]_val_[12].fq.gz//g' | uniq)
for sample in $samples 
do
gsnap -d $database -D . ${sample}_1_val_1.fq.gz ${sample}_2_val_2.fq.gz -N1 --gunzip -A sam -t $threads --read-group-name $sample --npaths 1 | samtools view -bS - > ${sample}.bam 
samtools sort ${sample}.bam ${sample}_sorted
samtools index ${sample}_sorted.bam
done

#merge the bam files into a single one
ls *_sorted.bam > bamfiles.list

samtools merge -@6 -b bamfiles.list merged_ww_rnaseq.bam


#Use Trinity to perform a genome-guided assembly

Trinity --max_memory 200G --CPU 30 --genome_guided_bam merged_ww_rnaseq.bam --genome_guided_max_intron 120000 --SS_lib_type RF --verbose --output ww_warbler_new_genome_guided_trinity



##########################################

#get transcripts using PASA 

##########################################

pasadir="/home/max/Downloads/PASApipeline-2.0.2/scripts"
transcriptdir="/home/max/Downloads/PASApipeline-2.0.2/scripts"

#Combine de novo and genome-guided transcripts

cat ww_trinity.fasta ww_warbler_new_genome_guided_trinity.fasta > Trinity_combined.fasta

#Use the seqclean script to remove low-quality transcripts

PASApipeline-2.0.2/seqclean/seqclean/seqclean Trinity_combined.fasta -c6


#extract IDs from trinity transcripts
cat Trinity.fasta | /home/max/Downloads/PASApipeline-2.0.2/misc_utilities/accession_extractor.pl > tdn.accs

cd PASApipeline-2.0.2/scripts/

./Launch_PASA_pipeline.pl --ALIGNERS gmap --MAX_INTRON_LENGTH 120000 -c ../pasa_conf/alignAssembly.config -g ww_genome_gapclosed_cleaned.fasta.masked.fasta -t Trinity_combined.fasta.clean -u Trinity_combined.fasta -T -R --TDN tdn.accs --CPU 6 --transcribed_is_aligned_orient

#Extract ORFs using transdecoder

./pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta ww_6.assemblies.fasta --pasa_transcripts_gff3 ww_6.pasa_assemblies.gff3


The pasa transdecoder transcripts were imported into webapollo and compared to synteny-transferred 
chicken genes (using satsuma and kraken) to evaluate their completeness.
From this data set we selected transcripts from ~1200 genes for training. 


##########################################

#get willow warbler parameters for augustus

##########################################


## Get parameters for augustus

Use scripts provided with `augustus` to filter the preliminary gene set and to perform the training


new_species.pl --species ww_old_genome_webapollo_1200

genome="ww_genome_gapclosed_cleaned.fasta"
distance=1000

#We must filter the webapollo names 
cat Annotations.ww_training_old_genome.prot.fasta | perl -lane 'if($F[0]=~/^>(\S+)/){print ">$1"} else{print @F}'  > Annotations.ww_training_old_genome.prot.filt_names.fasta

#Remove transcripts that have >80 % sequence similarity
aa2nonred.pl Annotations.ww_training_old_genome.prot.filt_names.fasta Annotations.ww_training_old_genome.prot.filt.id.fasta
gff2gbSmallDNA.pl Annotations.ww_training_old_genome.gff3 $genome $distance Annotations.ww.training_old_genome.gb

#Remove genes that are found to be problematic in the training - 1228 genes survive
etraining --species=generic Annotations.ww.training_old_genome.gb 2>&1 | grep "n sequence" | perl -pe 's/.*n sequence (\S+):.*/$1/' | sort -u > bad.ww.genes.lst
filterGenes.pl bad.ww.genes.lst Annotations.ww.training_old_genome.gb > Annotations.ww.training_old_genome.filt.gb

#Randomly split the gene set into a training and a test set
randomSplit.pl  Annotations.ww.training_old_genome.filt.gb 200

etraining --species=ww_old_genome_webapollo_1200 Annotations.ww.training_old_genome.filt.gb.train

optimize_augustus.pl --UTR=on --cpus=8 --species=ww_old_genome_webapollo_1200 Annotations.ww.training_old_genome.filt.gb.train









