#! /usr/bin/bash
#SBATCH --job-name=t1_samstack
#SBATCH --output=job.%x.%j_%A_%a.%N.out
#SBATCH --error=job.%x.%j_%A_%a.%N.err
#SBATCH --partition=C032M0128G
#SBATCH --qos=low
#SBATCH --get-user-env
#SBATCH -n 4
#SBATCH --cpu-freq=high
#SBATCH --time=120:00:00

set -e -u -o pipefail


####################
# info.txt
# number type label gene lysine sgRNA a.first.ref R1 R2
# number: id of data
# type: experiment or control
# gene: gene
# lysine: number of lysine in the protein to be edited
# sgRNA: sgRNA sequence
# a.first.ref: first A position in the reference file
# R1: file path of FASTQ R1 file
# R2: file path of FASTQ R2 file


DIR_WK=$(pwd)
DIR_CODE=${DIR_WK}/code
DIR_REF=${DIR_WK}/ref
DIR_RFQ=${DIR_WK}/data/rawdata/fq
DIR_CFQ=${DIR_WK}/data/cleanfq && mkdir -p ${DIR_CFQ}
DIR_SAM=${DIR_WK}/data/sam && mkdir -p ${DIR_SAM}
DIR_TXT=${DIR_WK}/data/txt && mkdir -p ${DIR_TXT}
DIR_SAMSTACK=${DIR_WK}/data/samstack && mkdir -p ${DIR_SAMSTACK}
DIR_TMP=${DIR_WK}/tmp && mkdir -p ${DIR_TMP}

lab=$(awk -F "\t" -v line=${SLURM_ARRAY_TASK_ID} 'FNR == (line + 1) {print $1;}' ${DIR_WK}/info.txt)

FQ1=$(awk -F "\t" -v line=${SLURM_ARRAY_TASK_ID} 'FNR == (line + 1) {print $8;}' ${DIR_WK}/info.txt)
FQ2=$(awk -F "\t" -v line=${SLURM_ARRAY_TASK_ID} 'FNR == (line + 1) {print $9;}' ${DIR_WK}/info.txt)

# bwa index of FASTA reference
FA=$(awk -F "\t" -v line=${SLURM_ARRAY_TASK_ID} 'FNR == (line + 1) {print $1;}' ${DIR_WK}/info.txt)

# unzip
echo "[[" `date` "]]----Unzipping"

zcat ${FQ1} > ${DIR_TMP}/${lab}.1.fq
zcat ${FQ2} > ${DIR_TMP}/${lab}.2.fq

echo "[[" `date` "]]----Unzip finished"

# cutadapt
echo "[[" `date` "]]----Cutadapt"

forward=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
backward=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA

cutadapt -a ${forward} -A ${backward} \
    -u 10 \
    -q 10,20 \
    -o ${DIR_CFQ}/${lab}.1.clean.fq \
    -p ${DIR_CFQ}/${lab}.2.clean.fq \
    ${DIR_TMP}/${lab}.1.fq ${DIR_TMP}/${lab}.2.fq

echo "[[" `date` "]]----Cutadapt finished"

# BWA

echo "[[" `date` "]]----BWA"

# gunzip ${DIR_CFQ}/${lab}.1.clean.fq.gz
# gunzip ${DIR_CFQ}/${lab}.2.clean.fq.gz

bwa mem -t 4  ${DIR_REF}/${FA}.fa  \
    ${DIR_CFQ}/${lab}.1.clean.fq \
    ${DIR_CFQ}/${lab}.2.clean.fq | \
    samtools view -b - | \
    samtools sort -O bam > ${DIR_SAM}/${lab}.sorted.bam

samtools index ${DIR_SAM}/${lab}.sorted.bam

echo "[[" `date` "]]----BWA finished"

# zip

echo "[[" `date` "]]----Zipping"

gzip ${DIR_CFQ}/${lab}.1.clean.fq
gzip ${DIR_CFQ}/${lab}.2.clean.fq

echo "[[" `date` "]]----Zip finished"


# # Split bam and stack

echo "[[" `date` "]]----Split and stack"

samtools view -F 12 -q 30 ${DIR_SAM}/${lab}.sorted.bam ${FA} | cut -f1-11 > ${DIR_TXT}/${lab}.txt

python3 ${DIR_CODE}/samstack.py -r ${DIR_REF}/${FA}.fa \
    -s ${DIR_TXT}/${lab}.txt -o ${DIR_SAMSTACK}/${lab}

echo "[[" `date` "]]----Split and stack finished"

################################################################################
