#! /bin/bash

dir="."
data=""
reference="hg38"
stage="qc"
include_qc=0
retain=0
gzip=0
thread=32
block_size=16000000
bin_size=100000
refdir="/path/to/reference"
tmpdir="/tmp"
fifo="multithread_processing"

while getopts "d:f:r:s:qtgp:b:z:x:" params; do
	case $params in
		d) dir="$OPTARG" ;;
		f) data="$OPTARG" ;;
		r) reference="$OPTARG" ;;
		s) stage="$OPTARG" ;;
		q) include_qc=1 ;;
		t) retain=1 ;;
		g) gzip=1 ;;
		p) thread="$OPTARG" ;;
		b) bin_size="$OPTARG" ;;
		z) block_size="$OPTARG" ;;
		x) tmpdir="$OPTARG" ;;
		*) echo "exit" ;;
	esac
done

if [[ $thread -lt 2 ]]; then
	thread=2
fi

case $stage in
	qc) start=1 ;;
	qc-only) start=1 ;;
	trim) start=2 ;;
	trim-only) start=2 ;;
	transform) start=3 ;;
	transform-only) start=3 ;;
	alignment) start=4 ;;
	alignment-only) start=4 ;;
	merge) start=5 ;;
	merge-only) start=5 ;;
	call-methylation) start=6 ;;
	call-methylation-only) start=6 ;;
	statistics) start=7 ;;
	clean) start=8 ;;
esac

# procedure
if [ -d $dir ]; then
	# information output
	start_time=$(date "+%s")
	echo $(date "+%Y-%m-%d %H:%M:%S") "pepline start"
	echo "Working directory: $dir"

	case $reference in
	hg38)
		index="${refdir}/hg38.fa_bowtie2/"
		reference="${refdir}/hg38.fa"
	;;
	mm39)
		index="${refdir}/mm39.fa_bowtie2/"
		reference="${refdir}/mm39.fa"
	;;
	*)
		echo "unrecognized reference genome $reference; must be one of hg38, mm39!"
		exit 1
	;;
	esac

	echo "Reference used: $reference"
	echo "index directory: $index"
	echo "reference path: $reference"

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
		echo "qc process completed"
		exit 0
	fi

	# trim adaptor
	if [[ $start -gt 2 ]]; then
		echo "trim adaptor process skipped"
	elif [ -d $dir/fastq ]; then 
		echo $(date "+%Y-%m-%d %H:%M:%S") "trim adaptor started"

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
					echo $(date "+%Y-%m-%d %H:%M:%S") "Running: trim_galore --paired --trim-n --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --gzip \$wd/fastq/$file \$wd/fastq/${file%R1*}R2${file##*R1} -o \$wd/trim >/dev/null 2>&1 &"
					trim_galore --paired --trim-n --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --gzip $dir/fastq/$file $dir/fastq/${file%R1*}R2${file##*R1} -o $dir/trim >/dev/null 2>&1
					echo -n "pp" >&8
				} &
			elif [[ $file == *.gz ]] && [[ $file != *R2*.gz ]]; then
				read -u 8 -n 2; {
					echo $(date "+%Y-%m-%d %H:%M:%S") "Running: trim_galore --trim-n --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --gzip \$wd/fastq/$file -o \$wd/trim >/dev/null 2>&1 &"
					trim_galore --trim-n --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --gzip $dir/fastq/$file -o $dir/trim >/dev/null 2>&1
					echo -n "pp" >&8
				} &
			fi
		done
		wait

		exec 8>&-
		rm $dir/$fifo

		echo $(date "+%Y-%m-%d %H:%M:%S") "trim adaptor finished"
	else
		echo "the data directory $dir/fastq don't exist"
	fi

	if [[ $stage == "trim-only" ]]; then
		end_time=$(date "+%s")
		elapsed_time=$((end_time - start_time))
		echo "trim adaptor process completed"
		echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
		exit 0
	fi

	# transform
	if [[ $start -gt 3 ]]; then
		echo "transform process skipped"
	elif [ -d $dir/trim ]; then 
		echo $(date "+%Y-%m-%d %H:%M:%S") "trim adaptor started"

		if [ -d $dir/transform ]; then
			rm -r $dir/transform
		fi
		mkdir $dir/transform

		if [ -e $dir/$fifo ]; then
			rm $dir/$fifo
		fi

		mkfifo -m 600 $dir/$fifo
		exec 8<>$dir/$fifo
		for ((i=1;i<=${thread};i++)); do
			echo -n "p" >&8
		done

		for file in `ls $dir/trim | sort`; do
			if [[ $file == *R1*.gz ]]; then
				read -u 8 -n 1; {
					echo $(date "+%Y-%m-%d %H:%M:%S") "Running: zcat \$wd/trim/$file | split -l 16000000 --additional-suffix=\"${file##*R1_val_1}\" --filter='gzip >\${FILE%.gz*}.gz' - \$wd/transform/${file%R1_val_1*}R1_val_1.part >/dev/null 2>&1 &"
					zcat $dir/trim/$file | split -l 16000000 --additional-suffix="${file##*R1_val_1}" --filter='gzip >${FILE%.gz*}.gz' - $dir/transform/${file%R1_val_1*}R1_val_1.part # >/dev/null 2>&1
					echo -n "p" >&8
				} &
				read -u 8 -n 1; {
					echo $(date "+%Y-%m-%d %H:%M:%S") "Running: zcat \$wd/trim/${file%R1_val_1*}R2_val_2${file##*R1_val_1} | split -l 16000000 --additional-suffix=\"${file##*R1_val_1}\" --filter='gzip >\${FILE%.gz*}.gz' - \$wd/transform/${file%R1_val_1*}R2_val_2.part >/dev/null 2>&1 &"
					zcat $dir/trim/${file%R1_val_1*}R2_val_2${file##*R1_val_1} | split -l 16000000 --additional-suffix="${file##*R1_val_1}" --filter='gzip >${FILE%.gz*}.gz' - $dir/transform/${file%R1_val_1*}R2_val_2.part >/dev/null 2>&1
					echo -n "p" >&8
				} &

			elif [[ $file == *.gz ]] && [[ $file != *R2*.gz ]]; then
				read -u 8 -n 1; {
					echo $(date "+%Y-%m-%d %H:%M:%S") "Running: zcat \$wd/trim/$file | split -l 16000000 --additional-suffix=\"${file##*trimmed}\" --filter='gzip >\${FILE%.gz*}.gz' - \$wd/transform/${file%trimmed*}trimmed.part >/dev/null 2>&1 &"
					zcat $dir/trim/$file | split -l 16000000 --additional-suffix="${file##*trimmed}" --filter='gzip >${FILE%.gz*}.gz' - $dir/transform/${file%trimmed*}trimmed.part >/dev/null 2>&1
					echo -n "p" >&8
				} &
			fi
		done
		wait

		for file in `ls $dir/transform | sort`; do
			if [[ $file == *R2_val_2.part* ]]; then
				read -u 8 -n 1; {
					echo $(date "+%Y-%m-%d %H:%M:%S") "Running: Antisense.py -i \$wd/transform/$file -o \$wd/transform/${file%R2_val_2.part*}R2_val_2.part.antisense${file##*R2_val_2.part}  >/dev/null 2>&1 &"
					Antisense.py -i $dir/transform/$file -o $dir/transform/${file%R2_val_2.part*}R2_val_2.antisense.part${file##*R2_val_2.part} >/dev/null 2>&1
					rm $dir/transform/$file
					echo -n "p" >&8
				} &
			fi
		done
		wait

		exec 8>&-
		rm $dir/$fifo

		echo $(date "+%Y-%m-%d %H:%M:%S") "transform finished"
	else
		echo "the data directory $dir/fastq don't exist"
	fi

	if [[ $stage == "transform-only" ]]; then
		end_time=$(date "+%s")
		elapsed_time=$((end_time - start_time))
		echo "transform process completed"
		echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
		exit 0
	fi


	# alignment
	if [[ $start -gt 4 ]]; then
		echo "alignment process skipped"
	elif [ -d $dir/transform ]; then
		echo $(date "+%Y-%m-%d %H:%M:%S") "alignment started"

		if [ -d $dir/alignment ]; then
			rm -r $dir/alignment
		fi
		mkdir $dir/alignment

		for file in `ls $dir/transform | sort`; do
			echo $(date "+%Y-%m-%d %H:%M:%S") "Running: bs_seeker2-align.py -i \$wd/transform/$file -g $reference --aligner=bowtie2 --bt2-p $((thread / 2)) --bt2--end-to-end -m 0.1 --XSteve --temp_dir=$tmpdir -o \$wd/alignment/${file%.f*q*}_alignment.bam >/dev/null 2>&1"
			bs_seeker2-align.py -i $dir/transform/$file -g $reference --aligner=bowtie2 --bt2-p $((thread / 2)) --bt2--end-to-end -m 0.1 --XSteve --temp_dir=$tmpdir -o $dir/alignment/${file%.f*q*}_alignment.bam >/dev/null 2>&1
		done

		echo $(date "+%Y-%m-%d %H:%M:%S") "alignment finished"
	else
		echo "the data directory $dir/transform don't exist"
	fi

	if [[ $stage == "alignment-only" ]]; then
		end_time=$(date "+%s")
		elapsed_time=$((end_time - start_time))
		echo "alignment process completed"
		echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
		exit 0
	fi

	# merge and sort
	if [[ $start -gt 5 ]]; then
		echo "merge process skipped"
	elif [ -d $dir/alignment ]; then 
		echo $(date "+%Y-%m-%d %H:%M:%S") "trim adaptor started"

		if [ -d $dir/merge ]; then
			rm -r $dir/merge
		fi
		mkdir $dir/merge

		declare -a sample_list
		for file in `ls $dir/alignment | sort`; do
			if [[ $file == *R1_val_1.part*.bam ]]; then
				sample=${file%R1_val_1.part*}
			elif [[ $file == *R2_val_2.antisense.part*.bam ]]; then
				sample=${file%R2_val_2.antisense.part*}
			elif [[ $file == *trimmed.part*.bam ]]; then
				sample=${file%trimmed.part*}
			fi
			if [[ -n $sample ]] && [[ ! " ${sample_list[*]} " =~ " $sample " ]]; then
				sample_list[${#sample_list[@]}]=$sample
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

		for sample in ${sample_list[*]}; do
			read -u 8 -n 4; {
				echo $(date "+%Y-%m-%d %H:%M:%S") "Running: samtools merge -f -@ 4 \$wd/alignment/${sample}merge.bam \$wd/alignment/${sample}*.bam >/dev/null 2>&1 &"
				samtools merge -f -@ 4 $dir/merge/${sample}merge.bam $dir/alignment/${sample}*.bam >/dev/null 2>&1
				echo $(date "+%Y-%m-%d %H:%M:%S") "Running: samtools sort -@ 4 --reference $reference -o \$wd/merge/${sample}sort.bam \$wd/merge/${sample}merge.bam >/dev/null 2>&1 &"
				samtools sort -@ 4 --reference $reference -o $dir/merge/${sample}sort.bam $dir/merge/${sample}merge.bam >/dev/null 2>&1
				echo -n "pppp" >&8
			} &
		done
		wait

		exec 8>&-
		rm $dir/$fifo

		echo $(date "+%Y-%m-%d %H:%M:%S") "merge finished"
	else
		echo "the data directory $dir/fastq don't exist"
	fi

	if [[ $stage == "merge-only" ]]; then
		end_time=$(date "+%s")
		elapsed_time=$((end_time - start_time))
		echo "merge process completed"
		echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
		exit 0
	fi

	# call methylation
	if [[ $start -gt 6 ]]; then
		echo "call methylation process skipped"
	elif [ -d $dir/merge ]; then
		echo $(date "+%Y-%m-%d %H:%M:%S") "call methylation started"

		if [ -d $dir/methylation ]; then
			rm -r $dir/methylation
		fi
		mkdir $dir/methylation

		if [ -e $dir/$fifo ]; then
			rm $dir/$fifo
		fi

		mkfifo -m 600 $dir/$fifo
		exec 8<>$dir/$fifo
		for ((i=1;i<=${thread};i++)); do
			echo -n "p" >&8
		done

		for file in `ls $dir/merge | sort`; do
			if [[ $file == *sort.bam ]]; then
				read -u 8 -n 1; {
					echo $(date "+%Y-%m-%d %H:%M:%S") "Running: bs_seeker2-call_methylation.py -i \$wd/merge/$file -d $index -o \$wd/methylation/${file%_*} --rm-overlap --sorted >/dev/null 2>&1 &"
					bs_seeker2-call_methylation.py -i $dir/merge/$file -d $index -o $dir/methylation/${file%_*} --rm-overlap --sorted >/dev/null 2>&1
					echo -n "p" >&8
				} &
			fi
		done
		wait

		exec 8>&-
		rm $dir/$fifo

		echo $(date "+%Y-%m-%d %H:%M:%S") "call methylation finished"
	else
		echo "the data directory $dir/merge don't exist"
	fi

	if [[ $stage == "call-methylation-only" ]]; then
		end_time=$(date "+%s")
		elapsed_time=$((end_time - start_time))
		echo "call methylation process completed"
		echo "process elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"
		exit 0
	fi

	# statistics
	if [[ $start -gt 7 ]]; then
		echo "methylation statistics process skipped"
	elif [ -d $dir/methylation ]; then
		echo $(date "+%Y-%m-%d %H:%M:%S") "methylation statistics started"

		if [ -d $dir/average ]; then
			rm -r $dir/average
		fi
		mkdir $dir/average

		if [ -e $dir/$fifo ]; then
			rm $dir/$fifo
		fi

		mkfifo -m 600 $dir/$fifo
		exec 8<>$dir/$fifo
		for ((i=1;i<=${thread};i++)); do
			echo -n "p" >&8
		done

		for file in `ls $dir/methylation | sort`; do
			if [[ $file == *.CGmap.gz ]]; then
				read -u 8 -n 1; {
					echo $(date "+%Y-%m-%d %H:%M:%S") "Running: CGmapMethInBins -i $dir/methylation/$file -B $bin_size -c 1  -p $dir/average/${file%.CGmap.gz*}. >$dir/average/${file%.CGmap.gz*}.average 2>&1 &"
					CGmapMethInBins -i $dir/methylation/$file -B $bin_size -c 1  -p $dir/average/${file%.CGmap.gz*}. >$dir/average/${file%.CGmap.gz*}.average 2>&1
					echo -n "p" >&8
				} &
			fi
		done
		wait

		exec 8>&-
		rm $dir/$fifo

		echo $(date "+%Y-%m-%d %H:%M:%S") "methylation statistics finished"
	else
		echo "the data directory $dir/methylation don't exist"
	fi

	# clean

	# done
	end_time=$(date "+%s")
	echo $(date "+%Y-%m-%d %H:%M:%S") "pepline done"
	elapsed_time=$((end_time - start_time))
	echo "pepline elapsed" $((elapsed_time / 3600)) "hours" $((elapsed_time / 60 % 60)) "minutes" $((elapsed_time % 60)) "seconds"

else
	echo "the directory $dir don't exist!"
fi

