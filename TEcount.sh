#! /bin/bash
set -e

dir="."
declare -a input
data=""
stage="qc"
include_qc=0
gene_length=""
gzip=0
thread=12
fifo="multithread_processing"
retain=0
reference="hg38"
output=""
    
while getopts "d:i:f:r:o:s:qtgp:" params; do
    case $params in
        d) dir="$OPTARG" ;;
        i) input[${#input[@]}]="$OPTARG" ;;
        f) data="$OPTARG" ;;
        r) reference="$OPTARG" ;;
        o) output="$OPTARG" ;;
        s) stage="$OPTARG" ;;
        q) include_qc=1 ;;
        t) retain=1 ;;
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
    quantity) start=4 ;;
    quantity-only) start=4 ;;
    postprocess) start=5 ;;
    postprocess-only) start=5 ;;
    clean) start=6 ;;
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
        reference_index="/utility/reference/RNAseq/hg38-TE/Homo_sapiens.GRCh38.TE_tran"
        reference_gtf="/utility/reference/gtf/Homo_sapiens.GRCh38.chrN.gtf"
        TE_gtf="/utility/reference/gtf/Homo_sapiens.GRCh38.TE.chrN.gtf"
    ;;
    mm39)
        echo "Reference: mm39"
        reference_index="/utility/reference/RNAseq/mm39-TE/Mus_musculus.GRCm39.TE_tran"
        reference_gtf="/utility/reference/gtf/Mus_musculus.GRCm39.chrN.gtf"
        TE_gtf="/utility/reference/gtf/Mus_musculus.GRCm39.TE.chrN.gtf"
    ;;
    *)
        echo "unrecognized reference genome $reference; must be one of hg38, mm39!"
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
                echo $(date "+%Y-%m-%d %H:%M:%S") "Running: hisat2 -p $thread --dta -k 100 -x $reference_index -1 $dir/trim/$file -2 $dir/trim/${file%R1_val_1*}R2_val_2${file##*R1_val_1} --summary-file $dir/alignment/${file%%_R1*}_alignment.txt 2>/dev/null | samtools sort -@ $thread -o $dir/alignment/${file%%_R1*}_sorted.bam >/dev/null 2>&1"
                hisat2 -p $thread --dta -k 100 -x $reference_index -1 $dir/trim/$file -2 $dir/trim/${file%R1_val_1*}R2_val_2${file##*R1_val_1} --summary-file $dir/alignment/${file%%_R1*}_alignment.txt 2>/dev/null | samtools sort -@ $thread -o $dir/alignment/${file%%_R1*}_sorted.bam >/dev/null 2>&1
            elif [[ $file == *.gz ]] && [[ $file != *R2*.gz ]]; then
                echo $(date "+%Y-%m-%d %H:%M:%S") "Running: hisat2 -p $thread --dta -k 100 -x $reference_index -U $dir/trim/$file --summary-file $dir/alignment/${file%_*}_alignment.txt 2>/dev/null | samtools sort -@ $thread -o $dir/alignment/${file%%_*}_sorted.bam >/dev/null 2>&1"
                hisat2 -p $thread --dta -k 100 -x $reference_index -U $dir/trim/$file --summary-file $dir/alignment/${file%_*}_alignment.txt 2>/dev/null | samtools sort -@ $thread -o $dir/alignment/${file%%_*}_sorted.bam >/dev/null 2>&1
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

    # quantity
    if [[ $start -gt 4 ]]; then
        echo "quantity process skipped"
    elif [[ -d $dir/alignment ]]; then 
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
            if [[ $file == *_sorted.bam ]]; then
                read -u 8 -n 1; {
                    echo $(date "+%Y-%m-%d %H:%M:%S") "Running: TEcount -b $dir/alignment/$file --GTF $reference_gtf --TE $TE_gtf --mode multi --project ${file%_sorted*}_count --outdir $dir/quantity --sortByPos >/dev/null 2>&1"
                    TEcount -b $dir/alignment/$file --GTF $reference_gtf --TE $TE_gtf --mode multi --project ${file%_sorted*}_count --outdir $dir/quantity --sortByPos >/dev/null 2>&1
                    echo -n "p" >&8
                } &
            fi
        done
        wait

        exec 8>&-
        rm $dir/$fifo

        echo $(date "+%Y-%m-%d %H:%M:%S") "quantity finished"
    else
        echo "the data directory $dir/alignment don't exist"
    fi

    if [[ $stage == "quantity-only" ]]; then
        end_time=$(date "+%s")
        elapsed_time=$((end_time - start_time))
        echo "quantity process completed"
        echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
        exit 0
    fi

    # postprocess
    if [[ $start -gt 5 ]]; then
        echo "postprocess process skipped"
    elif [[ -d $dir/quantity ]]; then 
        echo $(date "+%Y-%m-%d %H:%M:%S") "postprocess started"

        echo $(date "+%Y-%m-%d %H:%M:%S") "Running: python /utility/anaconda/envs/RNAseq/bin/TE_count_merge.py --input-directory \$wd/quantity/ --output-directory \$wd/result/"
        python /utility/anaconda/envs/RNAseq/bin/TE_count_merge.py --input-directory $dir/quantity/ --output-directory $dir/result/

        mv $dir/result/TE_count.csv $dir/result/${output}_TE_count.csv
        mv $dir/result/TE_count.txt $dir/result/${output}_TE_count.txt

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

    echo "all process completed"

    # clean
    if [[ $retain -eq 0 ]]; then
        if [ -d $dir/trim ]; then
            rm -r $dir/trim
        fi
        if [ -d $dir/alignment ]; then
            rm -r $dir/alignment
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
