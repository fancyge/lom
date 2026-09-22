# Longitudinal twins oral microbiome
**This file serves to display the general commands used to process childhood oral microbiome in a longitudinal twin cohort; statistical codes are in the scripts folder**


# Pre-processing of the sequencing reads
Summary of the datasets:

Totally 964 samples from 380 individuals and 189 families

Sample numbers for each time point : T1 - 271; T2 - 313; T3 - 380

Each sample has paired end reads

## 1. Quality Control
To check the quality of raw sequencing reads using FastQC for each sample:

#### Tools
fastqc/0.11.7s
#### Commands 
```
cd /project
mkdir FastQC
mkdir Logs/FastQC
fastqc -t 1 --extract -o ./FastQC/ ./Fastq/${sample}_R1.fastq.gz >> ./Logs/FastQC/${sample}_R1.fastq.gz.log  2>&1
fastqc -t 1 --extract -o ./FastQC/ ./Fastq/${sample}_R2.fastq.gz >> ./Logs/FastQC/${sample}_R2.fastq.gz.log  2>&1
```

## 2. Host reads removal
Remove human reads from each sample using reference sequence Homo_sapiens_assembly38.fasta. 
#### Tools
bbmap/38.93
#### Commands 
Quality-based masking:
```
time bbmask.sh in=./Reference/Homo_sapiens_assembly38.fasta out=./Reference/Homo_sapiens_assembly38.bbMasked.fasta overwrite=t
```		
Index the masked reference:
```
time bbmap.sh -Xmx161g ref=./Reference/Homo_sapiens_assembly38.bbMasked.fasta
```
```
bbmap.sh -Xmx42g \
    in=./Fastq/${sample}_R1.fastq.gz \
    in2=./Fastq/${sample}_R2.fastq.gz \
    outu=./Target_reads/Interleaved/${sample}.targetReads.interleaved.fastq.gz \
    ref=./Reference/Homo_sapiens_assembly38.bbMasked.fasta \
    threads=12 \
    overwrite=t \
    unpigz=t \
    usejni=t > ./Logs/bbmap/${sample}.bbmap.log
```

# Analysis on metagenomic reads 

## 1. Species taxonomic profiling using MetaPhlAn
To profile relative abundance of microbial communities
#### Tools
MetaPhlAn/3.0.13
samtools/1.12 
bowtie2/2.3.5.1

#### Inputs
paired end clean reads from 920 samples:

{sample}.targetReads.R1.fastq.gz; {sample}.targetReads.R2.fastq.gz

#### Commands
```
git clone https://github.com/biobakery/MetaPhlAn.git
cd MetaPhlAn
python -m pip install .
```
Install bowtie2db
```
metaphlan --install --bowtie2db ./Bowtie_DB
bowtie2-build --quiet --threads 4 -f ./Bowtie_DB/mpa_v30_CHOCOPhlAn_201901.fna ./Bowtie_DB/mpa_v30_CHOCOPhlAn_201901
```
Run MetaPhlAn
```
metaphlan ./Targeted_reads/${sample}.targetReads.R1.fastq.gz,./Targeted_reads/${sample}.targetReads.R2.fastq.gz --input_type fastq -s ./metaphlan/sams/${sample}.sam.bz2 --bowtie2out ./bowtie2/${sample}.bowtie2.bz2 -o ./profiles/${sample}_profile.tsv --bowtie2db ./Bowtie_DB --nproc 10 > ./Logs/mpa/${sample}.metaphlan.log
```

merge profile for all the samples
```
merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt
```
Generate the species level and genus level abundance table
```
grep -E "s__|clade" merged_abundance_table.txt | sed 's/^.*s__//g'\| cut -f1,3-933 > merged_abundance_table_species.txt
grep -v "|s__" merged_abundance_table_genus_sp.txt > merged_abundance_table_genus.txt
```

ordination plot using weighted unifrac distance generated from metaphlan abundance table
```
Rscript calculate_unifrac.R merged_abundance_table.txt mpa_v30_CHOCOPhlAn_201901_species_tree.nwk unifrac_merged_mpa3_profiles_weighted.tsv log10
```


## 2. Functional pathways profiling
To profile the abundance of functional pathways from microbial communities for each sample. 

#### Tools
HUMAnN3

#### Inputs
Interleaved reads

#### Commands 
```
humann --threads ${NCPUS} \
	--input ${sample}.targetReads.interleaved.fastq \
	--taxonomic-profile ./Metaphlan/profiles/${sample}_profile.tsv\
	--output ./Function_humann/${sample} \
	--metaphlan-options "--bowtie2db ./Bowtie_DB" >>${logfile} 2>&1
```
For each sample, the output includes 3 files: ${sample}_genefamilies.tsv; ${sample}_pathabundance.tsv; ${sample}_pathcoverage.tsv

Join the output files (gene families, coverage, and abundance) from all samples into three files
```
humann_join_tables -i ./Function_humann/genefamilies_pathabundance_pathcoverage -o joined_genefamilies.tsv --file_name genefamilies
humann_join_tables -i ./Function_humann/genefamilies_pathabundance_pathcoverage -o joined_pathabundance.tsv --file_name pathabundance
humann_join_tables -i ./Function_humann/genefamilies_pathabundance_pathcoverage -o joined_pathcoverage.tsv --file_name pathcoverage
```

Normalise joined pathway abundance files to copies per million (cpm); Split merged file into unstratified and stratified files
```
humann_renorm_table --input joined_pathabundance.tsv --output merged_CPM_pathabundance.tsv --units cpm --update-snames 
humann_split_stratified_table --input merged_CPM_pathabundance.tsv --output merged_CPM_pathabundance
```

## 3. Strain-level analysis for selected species

#### Tools
StrainPhlAn3

#### Inputs
Paired end reads

#### Commands

Obtain sam output files from MetaPhlAn 
```
metaphlan ./Targeted_reads/${sample}.targetReads.R1.fastq.gz,./Targeted_reads/${sample}.targetReads.R2.fastq.gz --input_type fastq -s ./metaphlan/sams/${sample}.sam.bz2 --bowtie2out ./bowtie2/${sample}.bowtie2.bz2 -o ./profiles/${sample}_profile.tsv --bowtie2db ./Bowtie_DB --nproc 10 > ./Logs/mpa/${sample}.metaphlan.log
```

For each sample, generate a marker file using sample_to_markers script. The marker file (.pkl) will contain the consensus of unique marker genes for each species found in that sample

```
mkdir consensus_markers
sample2markers.py -i ./metaphlan/sams/${sample}.sam.bz2 -o consensus_markers --nproc 10 ## T1551A_7032008_D1_D2.pkl 
```

extract markers for selected species; such as for __Streptococcus mitis__
```
mkdir clade_markers
extract_markers.py -c s__Streptococcus_mitis -o clade_markers -d ./Bowtie_DB/mpa_v30_CHOCOPhlAn_201901.pkl
```
Run strainphlan to generate concatenated alignments and polymorphic rates information for each species
```
strainphlan -s consensus_markers/*.pkl\
	-m clade_markers/s__${species}.fna\
	-r reference_genomes/${species}/*.fna\
	-o output/output_${species}\
	-c s__${species}\
	--phylophlan_mode fast
	--nprocs 14 
	-d ./Bowtie_DB/mpa_v30_CHOCOPhlAn_201901.pkl 
	--marker_in_n_samples 20 
	--sample_with_n_markers 20 
	--debug > ./Logs/strainphlan/spl_${species}_20.log 2>&1
```

## 4. Gene compositional analysis of Streptococcus mitis

#### Tools
PanPhlAn3

#### Inputs
paired end clean reads

#### Commands

install software; download reference pangenome for Streptococcus mitis
```
git clone https://github.com/SegataLab/panphlan.git
panphlan_download_pangenome.py -i Streptococcus_mitis -o ./Panphlan
```
map samples against pangenome
```
panphlan_map.py -p ./panphlan/Streptococcus_mitis/Streptococcus_mitis_pangenome.tsv\
	--indexes ./panphlan/Streptococcus_mitis/Streptococcus_mitis\
	-i ./panphlan/Interleaved_reads/${sample}.targetReads.interleaved.fastq.gz\
	-o ./panphlan/map_results/result_map_${sample}_smitis.csv\
	--nproc 10 > ./Logs/panphlan_map_${sample}_smitis.log 2>&1
```

profile mapped results to generate presence/absence matrix of gene families
```
panphlan_profiling.py -i ./panphlan/map_results --o_matrix ./panphlan/result_profile_smitis.tsv -p ./panphlan/Streptococcus_mitis/Streptococcus_mitis_pangenome.tsv --add_ref
```

# Analysis on metagenome assembled genomes (MAGs)


## 1. Assemble contigs using clean reads
To assemble contigs from microbial communities using paired end reads for each sample

#### Tools
MEGAHIT v1.2.9

#### Inputs
Interleaved target reads

#### Commands 
```
megahit --12 ./Target_reads/{sample}.targetReads.interleaved.fastq.gz -o ./Assembly_Megahit/{sample} -t 6 --continue --out-prefix {sample} > ./Logs/Assembly_Megahit/{sample}.log 2>&1
```

## 2. Metagenome assemble 
To reconstruct single genomes from microbial communities using metagenome binning tool

#### Tools
MetaBAT2

#### Inputs
Assembled contigs from Megahit: {sample}.contigs.fa

Sorted bam file from Align_to_assembly: {sample}.sort.bam (refer to **align_reads_to_contigs.sh**)

#### Commands 
```
runMetaBat.sh ./Assembly_Megahit/${sample}/${sample}.contigs.fa ./Align_to_assembly/${sample}/${sample}.sort.bam > ../Logs/metabat/${sample}.log
```

## 3. Quality check of MAGs
To check the quality of all MAGs and select qualified MAGs based on criteria

#### Tools 
Checkm v1.1.3

#### Inputs
output folder from MetaBAT2

#### Commands:
```
checkm lineage_wf -x fa ./Metabat2/${sample}.contigs.fa.metabat-bins-* ./Checkm/checkm_results/${sample} --reduced_tree -t 8 --tab_table -f qa_bins_${sample}.txt > ../Logs/checkm/${sample}.log 2>&1
checkm qa --tab_table ./Checkm/Checkm_results/${sample}/lineage.ms ./Checkm/Checkm_results/${sample} -o 2 -f ./Checkm/Checkm_qa/qa_bins_${sample}.txt
```


## 4. Dereplication of qualified MAGs

#### Tools
dRep v3.2.2

#### Inputs
all qualified MAGs

#### Commands:
minimum completeness is 75%; use default ANI threshold for primary clutering (0.9) and secondary clustering (0.95)
```
dRep dereplicate ./output_drep_15393MAGs_comp75/ -g ./dereplicate_MAGs/15393_MAGs_path.txt --S_algorithm fastANI -p 40 -comp 75 --genomeInfo 15393_MAGs_Info.csv > dRep_15393MAGs_comp75.log
```

## 5. taxonomy assignment of dereplicated MAGs

#### Inputs
species-level genome bins (SGBs)

#### Tools
MASH v2.3

#### Commands:

download complete sequences from reference sequences from NCBI and genbank
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt 
awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary_refseq.txt > assembly_summary_refseq_complete_genomes.txt
sed -i 's/https/ftp/' assembly_summary_refseq_all_genomes_cmd.txt 
for next in $(cat assembly_summary_refseq_complete_genomes_cmd.txt); do wget -P refseq "$next"; done # 38917 seqs
```
create sketch

```
cd ./refseq/
mash sketch -o Refseq *.fna.gz
cd ./gb/
mash sketch -o genbank *.fna.gz
```
generate mash distance

```
mash dist .msh ./output_drep_15393MAGs_comp75/dereplicated_genomes/{magID} ./mash_result/mash_result_{magID}.txt
```
#### Outputs
mash result txt file for each MAG, with best-matched reference sequence from database, mash distances and pvalue.

take the best-matched accession number from each mash result file and obtain the taxonomy ID using esearch command
```
esearch -db assembly -query '{accessionID}' | esummary | xtract -pattern DocumentSummary -element Taxid
```


## 6. Phylogeny construction using species-level genome bins
  
#### Inputs
species-level genome bins (SGBs)

#### Tools
PhyloPhlAn v3

#### Commands:
```
phylophlan -i ./dereplicate_MAGs/output_drep_15393MAGs_comp75/dereplicated_genomes\
	-o output_965MAGs\
	--nproc 40\
	-d phylophlan\
	--diversity high --fast\
	--genome_extension .fa\
	-f ./Phylophlan_MAGs/references_config.cfg --verbose 2>&1 |tee phylophlan_968.log 
```


## 7. Replication rates assessment 

map the MAGs to coverage files; fit a linear regression to the log-transformed coverage values along the genome.

#### Inputs
selected MAGs of particular species

#### Tools
iRep v1.1.14

#### Commands:
```
samtools view -h -o ./irep/input/${sample}.sam ./irep/input/${sample}.sort.bam > ./logs/irep_${sample}.log 2>&1
iRep -f ./irep/input/${mag}.fa -s ./irep/input/${sample}.sam -o ./irep/output/${sample}.iRep.txt > ./logs/${sample}.iRep.log
```

## 8. phylogeny construction using SGBs
#### Inputs
SGBs

#### Tools
PhyloPhlAn v3

#### Commands
```
phylophlan -i ./dereplicate_MAGs/output_drep_15463MAGs_comp75/dereplicated_genomes -o output --nproc 40 -d phylophlan --min_num_markers 50 --diversity high --fast  --genome_extension .fa -f ./references_config.cfg --verbose 2>&1 |tee phylophlan.log 
raxmlHPC-PTHREADS-SSE3 -p 1989 -m PROTCATLG -T 10 -t output/dereplicated_genomes_resolved.tre -w output -s output/dereplicated_genomes_concatenated.aln -n dereplicated_genomes_refined.tre
```


