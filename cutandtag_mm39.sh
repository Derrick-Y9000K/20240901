ls *_1.fq.gz  > 0
ls *_2.fq.gz  > 1
ls *_1.fq.gz|cut -d"_" -f 1  > 2
ls *_1.fq.gz|cut -d"." -f 1  > 3
ls *_2.fq.gz|cut -d"." -f 1  > 4

read fq1raw < 0 
read fq2raw < 1
read sample < 2
read sample1 < 3
read sample2 < 4


trim_galore --trim1 --paired --dont_gzip $fq1raw $fq2raw 

fq1=${sample1}_val_1.fq
fq2=${sample2}_val_2.fq


bowtie2 -N 1 -L 25 -p 12  --un-conc ${sample}_filtered.fq -x /path/to/repeat/index -1 $fq1  -2 $fq2  -S mapped_and_unmapped.sam 2>> ${sample}_filter_log.txt
rm mapped_and_unmapped.sam
fq1_filtered=${sample}_filtered.1.fq
fq2_filtered=${sample}_filtered.2.fq
python2 fastqCombinePairedEnd.py $fq1_filtered $fq2_filtered
bowtie2 -t -q -p 12 -I 25 -X 900 -N 1 -L 25 --no-unal -x /path/to/bowtie2/reference/index -1 $fq1 -2 $fq2 -S ${sample}.sam 2>> ${sample}_mapping_log.txt


sam=${sample}.sam
samtools view -b -S $sam > ${sample}.tmp.bam
samtools sort ${sample}.tmp.bam -o ${sample}.sort.bam
picard MarkDuplicates \
    I=${sample}.sort.bam \
    O=${sample}_marked_duplicates.bam \
    REMOVE_DUPLICATES=true M=${sample}_metrics.txt AS=true
samtools index ${sample}_marked_duplicates.bam

bamCoverage --bam ${sample}_marked_duplicates.bam -o ${sample}.RPKM.bw \
    --binSize 10  \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2654621783

macs2 callpeak -t ${sample}_marked_duplicates.bam -f BAMPE --g mm --nomodel --nolambda -q 0.05 -n ${sample}


rm $fq1_filtered
rm $fq2_filtered

rm  ${sample}.sam  
rm ${sample}.tmp.bam
rm ${sample}.sort.bam
rm  ${sample}_filtered.1.fq_pairs_R1.fastq
rm ${sample}_filtered.2.fq_pairs_R2.fastq

rm ${sample1}_val_1.fq
rm ${sample2}_val_2.fq


rm 0
rm 1
rm 2
rm 3
rm 4
