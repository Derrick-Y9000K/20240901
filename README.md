# G4TotiSC

Except for specific introduction, all reference genomde are downloaded from Ensembl database (http://asia.ensembl.org/index.html), while reference gtf are downloaded from Gencode database (https://www.gencodegenes.org/). <br>
hg38 genome: https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz <br>
hg38 gtf: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz <br>
mm39 genome: https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz <br>
mm39 gtf: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.chr_patch_hapl_scaff.annotation.gtf.gz <br>


## 1. RNAseq
- environment construction <br>
Use anaconda to create a new environment with python=3.7. <br>
Install fastqc (v0.11.9) and multiqc (v1.13) via anaconda if necessary. <br>
Install trim_galore (v0.6.6) and stringtie (v2.1.7) via apt. <br>
Install hisat2 (v2.1.0) from https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download. <br>

- index build <br>
Use hisat2-build to build reference index under /path/to/reference. <br>

- script run <br>
Place all raw sequencing files under /run_directory/fastq. <br>
Run script: `bash RNAseq.sh -p ${num_threads} -d /run_directory -r mm39` <br>


## 2. RNA-TE
- environment construction <br>
Can use the same environment as RNAseq after TEtranscripts installed. <br>
Install TEtranscripts (v2.2.3) from github (https://github.com/mhammell-laboratory/TEtranscripts/). <br>
TE gtf is downloaded from https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/. <br>

- script run <br>
Place all raw sequencing files under /run_directory/fastq. <br>
Run script: `bash TEcount.sh -p ${num_threads} -d /run_directory -r mm39` <br>


## 3. PCA
- environment construction <br>
packages: python (v3.10), numpy (v1.26.4), pandas (v2.2.2), scikit-learn (v1.5.1), matplotlib (v3.8.4), adjustText (v0.7.3). <br>

- script run <br>
Provide infile and outfile filename and run the cell. <br>


## 4. PCA-3D
- to do <br>


## 5. ATAC-seq/CUTtag
- environment construction <br>
Install trim_galore (v0.6.6) and bowtie2 (v2.4.4) via apt. <br>
Install samtools (v1.13), bamCoverage (v3.5.4) and macs2 (v2.2.9.1) via anaconda. <br>
Install picard (v2.25.5) from github (https://github.com/broadinstitute/picard/releases/latest/download/picard.jar). <br>

- index build <br>
Use bowtie2-build to build reference index under /path/to/bowtie2/reference. <br>

- script run <br>
Place raw sequencing files of each sample under separate directory. <br>
Run script: `bash cutandtag_mm39.sh` <br>


## 6. DNA methylation
- environment construction <br>
Use anaconda to create a new environment with python=2.7. <br>
Install fastqc (v0.11.9) and multiqc (v1.13) via anaconda if necessary. <br>
Install samtools (v1.13) via anaconda. <br>
Install BSseeker2 (v2.1.8) and CGmaptools (v0.1.2) from github (git://github.com/BSSeeker/BSseeker2.git, git://github.com/guoweilong/cgmaptools.git). <br>

- index build <br>
Use bs_seeker2-build.py to build reference index under /path/to/reference. <br>

- script run <br>
Place all raw sequencing files under /run_directory/fastq. <br>
Run script: `bash call_methylation.sh -p ${num_threads} -d /run_directory -r mm39 -b 500` <br>


## 7. scRNA-seq
- environment construction <br>
Install cellranger (v7.2.0) from official website (https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome). <br>
packages: ggplot2 (), dplyr (), Seurat (), tidyr (), paletteer (). <br>

- script run <br>
Place all raw sequencing files under /run_directory/data and cd /run_directory. <br>
Rename all the files in [sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz format. <br>
Run upstream analysis: `cellranger count --id=$sample --localcores=16 --transcriptome=$db --fastqs=./data --sample=$sample --nosecondary --expect-cells=5000` <br>
Run downstream analysis: use ${sample}/outs/filtered_feature_bc_matrix/ to conduct downstream analysis with singlecell_downstream.R. <br>
