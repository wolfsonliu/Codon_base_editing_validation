#! /usr/bin/bash

function printusage {
    echo "Usage: $0 [-h] -r ref.fa -p input_1.fq -q input_2.fq -o result" 1>&2
}

function printhelp {
    printusage
    echo "" 1>&2
    echo "Parameters:" 1>&2
    echo "    h: print this help" 1>&2
    echo "    r: the indexed reference fasta file." 1>&2
    echo "    p: the input forward fastq file." 1>&2
    echo "    q: the input reversed fastq file." 1>&2
    echo "    o: the output prefix." 1>&2
    echo "" 1>&2
    echo "Input:" 1>&2
    echo "    The reference file should be index by bwa index." 1>&2
    echo "    The input file should be fastq files." 1>&2
    echo "    The output prefix should be path and basic file name label" 1>&2
    echo "    For example -o ./out results in ./out.bam ./out.vcf and so on" 1>&2
    echo "" 1>&2
}

while getopts "hr:p:q:o:" opt; do
  case ${opt} in
    h)
      printhelp; exit 0
      ;;
    r)
      REF_FA=${OPTARG}
      ;;
    p)
      FQ1=${OPTARG}
      ;;
    q)
      FQ2=${OPTARG}
      ;;
    o)
      OUT_PREFIX=${OPTARG}
      ;;
    \?)
      echo "Invalid option: -"${OPTARG}
      printusage; exit 1
      ;;
  esac
done


function splitfq_insertion {
    awk '$6 ~ /I/ {
        print $0;
    }' $1
}

function splitfq_deletion {
    awk '$6 ~ /D/ {
        print $0;
    }' $1
}

function splitfq_other {
    awk '$6 !~ /(I|D)/ {
        print $0;
    }' $1
}


bwa mem -t 4 ${REF_FA} ${FQ1} ${FQ2} \
    | samtools view -b - | samtools sort -O bam > ${OUT_PREFIX}.sorted.bam
samtools index ${OUT_PREFIX}.sorted.bam
# samtools view ${OUT_PREFIX}.sorted.bam | cut -f1-11 > ${OUT_PREFIX}.sorted.sam

# splitfq_insertion ${OUT_PREFIX}.sorted.sam > ${OUT_PREFIX}.insertion.sam
# splitfq_deletion ${OUT_PREFIX}.sorted.sam > ${OUT_PREFIX}.deletion.sam
# splitfq_other ${OUT_PREFIX}.sorted.sam > ${OUT_PREFIX}.other.sam

####################
