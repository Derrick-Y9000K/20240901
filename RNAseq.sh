#! /bin/bash
set -e

dir="."
declare -a input
data=""
stage="qc"
include_qc=0
pca=0
transcripts_reconstruct=1
gene_length=""
gzip=0
thread=12
fifo="multithread_processing"
retain=0
reference_dir="/path/to/reference"
reference="hg38"
output=""
    
while getopts "d:i:f:r:o:s:l:qatngp:" params; do
    case $params in
        d) dir="$OPTARG" ;;
        i) input[${#input[@]}]="$OPTARG" ;;
        f) data="$OPTARG" ;;
        r) reference="$OPTARG" ;;
        o) output="$OPTARG" ;;
        s) stage="$OPTARG" ;;
        q) include_qc=1 ;;
        a) pca=1 ;;
        t) retain=1 ;;
        n) transcripts_reconstruct=0 ;;
        l) gene_length="$OPTARG" ;;
        g) gzip=1 ;;
        p) thread="$OPTARG" ;;
        *) exit 1 ;;
    esac
done

if [[ $thread -lt 2 ]]; then
    thread=2
fi

dir=`realpath $dir`
if [[ ${!input[@]} -gt 0 ]]; then
    stage="multi-assemble"
    for i in ${!input[@]}; do
        input[i]=`realpath ${input[i]}`
    done
fi

case $stage in
    merge) start=0 ;;
    merge-only) start=0 ;;
    qc) start=1 ;;
    qc-only) start=1 ;;
    trim) start=2 ;;
    trim-only) start=2 ;;
    alignment) start=3 ;;
    alignment-only) start=3 ;;
    assemble) start=4 ;;
    multi-assemble) start=4 ;;
    assemble-only) start=4 ;;
    quantity) start=5 ;;
    quantity-only) start=5 ;;
    postprocess) start=6 ;;
    postprocess-only) start=6 ;;
    analysis) start=7 ;;
    clean) start=8 ;;
    *) echo "unrecognized stage $stage"; exit 1 ;;
esac

# procedure
if [ -d $dir ]; then
    # information output
    start_time=$(date "+%s")
    echo $(date "+%Y-%m-%d %H:%M:%S") "pepline start"
    echo "Working directory: $dir"

    case $reference in
    hg38)
        echo "Reference: hg38"
        echo "directory: ${reference_dir}/hg38/"
        reference_gtf="${reference_dir}/hg38/Homo_sapiens.GRCh38.107.chr.gtf"
        reference_index="${reference_dir}/hg38/Homo_sapiens.GRCh38_tran"
    ;;
    mm39)
        echo "Reference: mm39"
        echo "directory: ${reference_dir}/mm39/"
        reference_gtf="${reference_dir}/mm39/Mus_musculus.GRCm39.107.chr.gtf"
        reference_index="${reference_dir}/mm39/Mus_musculus.GRCm39_tran"
    ;;
    mm39-MuERV)
        echo "Reference: mm39-MuERV"
        echo "directory: ${reference_dir}/mm39-MuERV/"
        reference_gtf="${reference_dir}/mm39-MuERV/Mus_musculus.GRCm39.MuERV.chrN.gtf"
        reference_index="${reference_dir}/mm39-MuERV/Mus_musculus.GRCm39.MuERV_tran"
    ;;
    tdTamato)
        echo "Reference: tdTamato"
        echo "directory: ${reference_dir}/tdTamato/"
        reference_gtf="${reference_dir}/tdTamato/tdTamato.gtf"
        reference_index="${reference_dir}/tdTamato/tdTamato_tran"
    ;;
    *)
        echo "unrecognized reference genome $reference; must be one of hg38, mm39, mm39-MuERV, tdTamato!"
        exit 1
    ;;
    esac

    # pre-process
    if [[ $gzip -eq 1 ]] && [[ $start -le 2 ]]; then
        if [ $start -gt 0 ]; then
            find $dir/fastq -type f -name "*.f*q" | xargs -P $thread gzip
        else
            find $dir/$data -type f -name "*.f*q" | xargs -P $thread gzip
        fi
    fi

    if [[ $start -gt 0 ]]; then
        rename 's/-L.*_1\.f/_R1\.f/' $dir/fastq/*
        rename 's/-L.*_2\.f/_R2\.f/' $dir/fastq/*

        rename 's/_1\.f/_R1\.f/' $dir/fastq/*
        rename 's/_2\.f/_R2\.f/' $dir/fastq/*
    else
        rename 's/-L.*_1\.f/_R1\.f/' $dir/$data/*
        rename 's/-L.*_2\.f/_R2\.f/' $dir/$data/*

        rename 's/_1\.f/_R1\.f/' $dir/$data/*
        rename 's/_2\.f/_R2\.f/' $dir/$data/*
    fi

    # merge
    if [[ $start -gt 0 ]]; then
        echo "fastq file merge process skipped"
    elif [ -d $dir/$data ]; then
        echo $(date "+%Y-%m-%d %H:%M:%S") "fastq file merge started"

        if [ -d $dir/fastq ]; then
            rm -r $dir/fastq
        fi
        mkdir $dir/fastq

        declare -a single_end
        declare -a paired_end

        for file in `ls $dir/$data | sort`; do
            if [[ $file == *_R1*.gz ]]; then
                if [[ ! " ${paired_end[*]} " =~ " ${file%-*} " ]]; then
                    paired_end[${#paired_end[@]}]=${file%-*}
                fi
            elif [[ $file == *.gz ]] && [[ $file != *R*.gz ]]; then
                if [[ ! " ${single_end[*]} " =~ " ${file%-*} " ]]; then
                    single_end[${#single_end[@]}]=${file%-*}
                fi
            fi
        done

        if [ -e $dir/$fifo ]; then
            rm $dir/$fifo
        fi

        mkfifo -m 600 $dir/$fifo
        exec 8<>$dir/$fifo
        for ((i=1;i<=${thread};i++)); do
            echo -n "p" >&8
        done 

        for file in ${single_end[*]}; do
            read -u 8 -n 1; {
                echo $(date "+%Y-%m-%d %H:%M:%S") "Running: cat \$wd/$data/${file}*.f*q.gz >\$wd/fastq/$file.fastq.gz &"
                cat $dir/$data/${file}*.f*q.gz >$dir/fastq/$file.fastq.gz
                echo -n "p" >&8
            } &
        done
        for file in ${paired_end[*]}; do
            read -u 8 -n 1; {
                echo $(date "+%Y-%m-%d %H:%M:%S") "Running: cat \$wd/$data/${file}*_R1.f*q.gz >\$wd/fastq/${file}_R1.fastq.gz &"
                cat $dir/$data/${file}*_R1.f*q.gz >$dir/fastq/${file}_R1.fastq.gz
                echo -n "p" >&8
            } &
            read -u 8 -n 1; {

                echo $(date "+%Y-%m-%d %H:%M:%S") "Running: cat \$wd/$data/${file}*_R2.f*q.gz >\$wd/fastq/${file}_R2.fastq.gz &"
                cat $dir/$data/${file}*_R2.f*q.gz >$dir/fastq/${file}_R2.fastq.gz
                echo -n "p" >&8
            } &
        done
        wait

        exec 8>&-
        rm $dir/$fifo

        echo $(date "+%Y-%m-%d %H:%M:%S") "fastq file merge finished"
    else
        echo "the data directory $dir/$data don't exist"
    fi

    if [[ $stage == "merge-only" ]]; then
        end_time=$(date "+%s")
        elapsed_time=$((end_time - start_time))
        echo "fastq file merge process completed"
        echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
        exit 0
    fi

    # qc
    if [[ $start -gt 1 ]] || [[ $include_qc -eq 0 ]]; then
        echo "qc process skipped"
    elif [ -d $dir/fastq ]; then 
        echo $(date "+%Y-%m-%d %H:%M:%S") "qc started"

        if [ -d $dir/qc ]; then
            rm -r $dir/qc
        fi
        mkdir $dir/qc

        echo $(date "+%Y-%m-%d %H:%M:%S") "Running: fastqc -t $thread -o \$wd/qc/ \$wd/fastq/*.gz >/dev/null 2>&1"
        fastqc -t $thread -o $dir/qc/ $dir/fastq/*.gz >/dev/null 2>&1
        echo $(date "+%Y-%m-%d %H:%M:%S") "Running: multiqc -o $dir/qc/ \$wd/qc/* >/dev/null 2>&1"
        multiqc -o $dir/qc/ $dir/qc/* >/dev/null 2>&1

        echo $(date "+%Y-%m-%d %H:%M:%S") "qc finished"
    else
        echo "the data directory $dir/fastq don't exist"
    fi

    if [[ $staege == "qc-only" ]]; then
        end_time=$(date "+%s")
        elapsed_time=$((end_time - start_time))
        echo "qc process completed"
        echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
        exit 0
    fi

    # trim adapter
    if [[ $start -gt 2 ]]; then
        echo "trim adapter process skipped"
    elif [ -d $dir/fastq ]; then 
        echo $(date "+%Y-%m-%d %H:%M:%S") "trim adapter started"

        if [ -d $dir/trim ]; then
            rm -r $dir/trim
        fi
        mkdir $dir/trim

        if [ -e $dir/$fifo ]; then
            rm $dir/$fifo
        fi

        mkfifo -m 600 $dir/$fifo
        exec 8<>$dir/$fifo
        for ((i=1;i<=${thread};i++)); do
            echo -n "p" >&8
        done

        for file in `ls $dir/fastq | sort`; do
            if [[ $file == *R1*.gz ]]; then
                read -u 8 -n 2; {
                    echo $(date "+%Y-%m-%d %H:%M:%S") "Running: trim_galore --paired --trim-n --gzip \$wd/fastq/$file \$wd/fastq/${file%R1*}R2${file##*R1} -o \$wd/trim >/dev/null 2>&1 &"
                    trim_galore --paired --trim-n --gzip $dir/fastq/$file $dir/fastq/${file%R1*}R2${file##*R1} -o $dir/trim >/dev/null 2>&1
                    echo -n "pp" >&8
                } &
            elif [[ $file == *.gz ]] && [[ $file != *R2*.gz ]]; then
                read -u 8 -n 2; {
                    echo $(date "+%Y-%m-%d %H:%M:%S") "Running: trim_galore --trim-n --gzip \$wd/fastq/$file -o \$wd/trim >/dev/null 2>&1 &"
                    trim_galore --trim-n --gzip $dir/fastq/$file -o $dir/trim >/dev/null 2>&1
                    echo -n "pp" >&8
                } &
            fi
        done
        wait

        exec 8>&-
        rm $dir/$fifo

        echo $(date "+%Y-%m-%d %H:%M:%S") "trim adapter finished"
    else
        echo "the data directory $dir/fastq don't exist"
    fi

    if [[ $stage == "trim-only" ]]; then
        end_time=$(date "+%s")
        elapsed_time=$((end_time - start_time))
        echo "trim adapter process completed"
        echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
        exit 0
    fi

    # alignment
    if [[ $start -gt 3 ]]; then
        echo "alignment process skipped"
    elif [ -d $dir/trim ]; then
        echo $(date "+%Y-%m-%d %H:%M:%S") "alignment started"

        if [ -d $dir/alignment ]; then
            rm -r $dir/alignment
        fi
        mkdir $dir/alignment

        for file in `ls $dir/trim | sort`; do
            if [[ $file == *R1_val_1*.gz ]]; then
                echo $(date "+%Y-%m-%d %H:%M:%S") "Running: hisat2 -p $thread --dta -x $reference_index -1 \$wd/trim/$file -2 \$wd/trim/${file%R1_val_1*}R2_val_2${file##*R1_val_1} --summary-file \$wd/alignment/${file%%_R1*}_alignment.txt 2>/dev/null | samtools sort -@ $thread -o \$wd/alignment/${file%%_R1*}_sorted.bam >/dev/null 2>&1"
                hisat2 -p $thread --dta -x $reference_index -1 $dir/trim/$file -2 $dir/trim/${file%R1_val_1*}R2_val_2${file##*R1_val_1} --summary-file $dir/alignment/${file%%_R1*}_alignment.txt 2>/dev/null | samtools sort -@ $thread -o $dir/alignment/${file%%_R1*}_sorted.bam >/dev/null 2>&1
            elif [[ $file == *.gz ]] && [[ $file != *R2*.gz ]]; then
                echo $(date "+%Y-%m-%d %H:%M:%S") "Running: hisat2 -p $thread --dta -x $reference_index -U \$wd/trim/$file --summary-file \$wd/alignment/${file%_*}_alignment.txt 2>/dev/null | samtools sort -@ $thread -o \$wd/alignment/${file%%_*}_sorted.bam >/dev/null 2>&1"
                hisat2 -p $thread --dta -x $reference_index -U $dir/trim/$file --summary-file $dir/alignment/${file%_*}_alignment.txt 2>/dev/null | samtools sort -@ $thread -o $dir/alignment/${file%%_*}_sorted.bam >/dev/null 2>&1
            fi
        done

        echo $(date "+%Y-%m-%d %H:%M:%S") "alignment finished"
    else
        echo "the data directory $dir/fastq don't exist"
    fi

    if [[ $stage == "alignment-only" ]]; then
        end_time=$(date "+%s")
        elapsed_time=$((end_time - start_time))
        echo "alignment process completed"
        echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
        exit 0
    fi

    # assemble
    if [[ $start -lt 6 ]] && [[ $transcripts_reconstruct -eq 0 ]]; then
        start=5
        quantity_gtf=$reference_gtf
    else
        quantity_gtf=$dir/assemble/merged.gtf
    fi

    if [[ $start -gt 4 ]]; then
        echo "assemble process skipped"
    elif [ -d $dir/alignment ] || [[ $stage == "multi-assemble" ]]; then 
        if [[ $stage == "multi-assemble" ]]; then
            if [[ ${#input[@]} -eq 0 ]]; then
                echo "missing parameter -i !"
                exit 1
            fi

            if [[ ! " ${input[*]} " =~ " $dir " ]] && [ -d $dir/assemble ]; then
                rm -r $dir/assemble
            fi
        fi

        echo $(date "+%Y-%m-%d %H:%M:%S") "assembly started"
        if [ ! -d $dir/assemble ]; then
            mkdir $dir/assemble
        fi

        if [ -e $dir/assemble/sample.txt ]; then
            rm $dir/assemble/sample.txt
        fi
        if [ -e $dir/assemble/merged.gtf ]; then
            rm $dir/assemble/merged.gtf
        fi

        if [[ $stage == "multi-assemble" ]]; then
            if [ ! -d $dir/alignment ]; then
                mkdir $dir/alignment
            fi

            for input_dir in ${input[*]}; do
                for file in `ls $input_dir/alignment | sort`; do
                    if [[ $file == *.bam ]] && [[ ! -e $dir/alignment/$file ]]; then
                        ln -s $input_dir/alignment/$file $dir/alignment/$file
                    fi
                done
                for file in `ls $input_dir/assemble | sort`; do
                    if [[ $file == *.gtf ]] && [[ $file != merged.gtf ]]; then
                        echo $input_dir/assemble/$file >>$dir/assemble/sample.txt
                    fi
                done
            done
        else
            if [ -e $dir/$fifo ]; then
                rm $dir/$fifo
            fi

            mkfifo -m 600 $dir/$fifo
            exec 8<>$dir/$fifo
            for ((i=1;i<=${thread};i++)); do
                echo -n "p" >&8
            done

            for file in `ls $dir/alignment | sort`; do
                if [[ $file == *.bam ]]; then
                    read -u 8 -n 2; {
                        echo $(date "+%Y-%m-%d %H:%M:%S") "Running: stringtie -p 2 -G $reference_gtf -l ${file%_*}_assemble -o \$wd/assemble/${file%_*}_assemble.gtf \$wd/alignment/$file >/dev/null 2>&1 &"
                        stringtie -p 2 -G $reference_gtf -l ${file%_*}_assemble -o $dir/assemble/${file%_*}_assemble.gtf $dir/alignment/$file >/dev/null 2>&1
                        echo -n "pp" >&8
                    } &
                fi
            done
            wait

            exec 8>&-
            rm $dir/$fifo

            for file in `ls $dir/assemble | sort`; do
                if [[ $file == *.gtf ]]; then
                    echo $dir/assemble/$file >>$dir/assemble/sample.txt
                fi
            done
        fi

        if [ -e $dir/assemble/sample.txt ]; then
            echo $(date "+%Y-%m-%d %H:%M:%S") "Running: stringtie --merge -p $thread -G $reference_gtf -o \$wd/assemble/merged.gtf \$wd/assemble/sample.txt"
            stringtie --merge -p $thread -G $reference_gtf -o $dir/assemble/merged.gtf $dir/assemble/sample.txt
        fi

        echo $(date "+%Y-%m-%d %H:%M:%S") "assembly finished"
    else
        echo "the data directory $dir/alignment don't exist"
    fi

    if [[ $stage == "assemble-only" ]]; then
        end_time=$(date "+%s")
        elapsed_time=$((end_time - start_time))
        echo "assembly process completed"
        echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
        exit 0
    fi

    # quantity
    if [[ $start -gt 5 ]]; then
        echo "quantity process skipped"
    elif ([ $transcripts_reconstruct -eq 0 ] || [[ -d $dir/assemble ]]) && [[ -d $dir/alignment ]]; then 
        echo $(date "+%Y-%m-%d %H:%M:%S") "quantity started"

        if [ -d $dir/quantity ]; then
            rm -r $dir/quantity
        fi
        mkdir $dir/quantity

        if [ -e $dir/$fifo ]; then
            rm $dir/$fifo
        fi

        mkfifo -m 600 $dir/$fifo
        exec 8<>$dir/$fifo
        for ((i=1;i<=${thread};i++)); do
            echo -n "p" >&8
        done

        for file in `ls $dir/alignment | sort`; do
            if [[ $file == *.bam ]]; then
                read -u 8 -n 2; {
                    echo $(date "+%Y-%m-%d %H:%M:%S") "Running: stringtie -p 2 -e -G $quantity_gtf -o \$wd/quantity/${file%_*}_quantity.gtf -A \$wd/quantity/${file%_*}_gene_abundances.tsv \$wd/alignment/$file >/dev/null 2>&1 &"
                    stringtie -p 2 -e -G $quantity_gtf -o $dir/quantity/${file%_*}_quantity.gtf -A $dir/quantity/${file%_*}_gene_abundances.tsv $dir/alignment/$file >/dev/null 2>&1
                    echo -n "pp" >&8
                } &
            fi
        done
        wait

        exec 8>&-
        rm $dir/$fifo

        echo $(date "+%Y-%m-%d %H:%M:%S") "quantity finished"
    else
        if [ $transcripts_reconstruct -eq 0 ]; then
            echo "the data directory $dir/alignment don't exist"
        else
            echo "the data directory $dir/alignment or $dir/assemble don't exist"
        fi
    fi

    if [[ $stage == "quantity-only" ]]; then
        end_time=$(date "+%s")
        elapsed_time=$((end_time - start_time))
        echo "quantity process completed"
        echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
        exit 0
    fi

    # postprocess
    if [[ $start -gt 6 ]]; then
        echo "postprocess process skipped"
    elif [[ -d $dir/quantity ]]; then 
        echo $(date "+%Y-%m-%d %H:%M:%S") "postprocess started"

        if [ -e $dir/quantity/sample.txt ]; then
            rm $dir/quantity/sample.txt
        fi
        for file in `ls $dir/quantity | sort`; do
            if [[ $file == *.gtf ]]; then
                echo ${file%_*} $dir/quantity/$file >>$dir/quantity/sample.txt
            fi
        done

        if [ -d $dir/result ]; then
            rm -r $dir/result
        fi
        mkdir $dir/result

        if [ -e $dir/quantity/sample.txt ]; then
            echo $(date "+%Y-%m-%d %H:%M:%S") "Running: prepDE.py -i \$wd/quantity/sample.txt -g \$wd/result/${output}_gene_count_matrix.csv -t \$wd/result/${output}_transcript_count_matrix.csv"
            prepDE.py -i $dir/quantity/sample.txt -g $dir/result/${output}_gene_count_matrix.csv -t $dir/result/${output}_transcript_count_matrix.csv
            echo $(date "+%Y-%m-%d %H:%M:%S") "Running: python /utility/anaconda/envs/RNAseq/bin/gene_count_filter.py --input $dir/result/${output}_gene_count_matrix.csv --output $dir/result/${output}_gene_count_matrix_filtered.txt"
            python gene_count_filter.py --input $dir/result/${output}_gene_count_matrix.csv --output $dir/result/${output}_gene_count_matrix_filtered.txt
        fi

        if [ -z "$gene_length" ]; then
            echo $(date "+%Y-%m-%d %H:%M:%S") "Running: GetGeneLength --database ensembl --gtffile $reference_gtf --lengthfile \$wd/result/gene_length.txt >/dev/null 2>&1"
            GetGeneLength --database ensembl --gtffile $reference_gtf --lengthfile $dir/result/gene_length.txt >/dev/null 2>&1
            gene_length=$dir/result/gene_length.txt
        fi

        echo $(date "+%Y-%m-%d %H:%M:%S") "Running: python normalization.py --directory \$wd/result --input ${output}_gene_count_matrix.csv --length $gene_length --rpkm-result ${output}_rpkm_gene_count.csv --tpm-result ${output}_tpm_gene_count.csv"
        python normalization.py --directory $dir/result --input ${output}_gene_count_matrix.csv --length $gene_length --rpkm-result ${output}_rpkm_gene_count.csv --tpm-result ${output}_tpm_gene_count.csv

        echo $(date "+%Y-%m-%d %H:%M:%S") "postprocess finished"
    else
        echo "the data directory $dir/quantity don't exist"
    fi

    if [[ $stage == "postprocess-only" ]]; then
        end_time=$(date "+%s")
        elapsed_time=$((end_time - start_time))
        echo "postprocess process completed"
        echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
        exit 0
    fi

    # pca analysis
    if [[ $start -gt 7 ]] || [[ $pca -eq 0 ]]; then
        echo "principal component analysis skipped"
    elif [[ -d $dir/result ]]; then 
        echo $(date "+%Y-%m-%d %H:%M:%S") "principal component analysis started"

        echo $(date "+%Y-%m-%d %H:%M:%S") "Running: python pca_analysis.py --directory \$wd/result --input ${output}_rpkm_gene_count.csv --separator , --output ${output}_rpkm_result.png"
        python pca_analysis.py --directory $dir/result --input ${output}_rpkm_gene_count.csv --separator , --output ${output}_rpkm_result.png

        echo $(date "+%Y-%m-%d %H:%M:%S") "Running: python pca_analysis.py --directory \$wd/result --input ${output}_tpm_gene_count.csv --separator , --output ${output}_tpm_result.png"
        python pca_analysis.py --directory $dir/result --input ${output}_tpm_gene_count.csv --separator , --output ${output}_tpm_result.png

        echo $(date "+%Y-%m-%d %H:%M:%S") "principal component analysis finished"
    else
        echo "the data directory $dir/result don't exist"
    fi

    echo "all process completed"

    # clean
    if [[ $retain -eq 0 ]]; then
        if [ -d $dir/trim ]; then
            rm -r $dir/trim
        fi
        if [ -d $dir/quantity ]; then
            rm -r $dir/quantity
        fi
        echo "redundant files deleted"
    fi

    # done
    end_time=$(date "+%s")
    echo $(date "+%Y-%m-%d %H:%M:%S") "pipeline done"
    elapsed_time=$((end_time - start_time))
    echo "pipeline elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"

else
    echo "the directory $dir don't exist!"
fi
