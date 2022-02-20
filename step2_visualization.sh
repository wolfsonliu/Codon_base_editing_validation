#! /usr/bin/bash
#SBATCH --job-name=t2_vis
#SBATCH --output=%x.%j_%A_%a.%N.out
#SBATCH --error=%x.%j_%A_%a.%N.err
#SBATCH --partition=C032M0128G
#SBATCH --qos=low
#SBATCH --get-user-env
#SBATCH -n 2
#SBATCH --time=120:00:00

set -e -u -o pipefail

DIR_WK=$(pwd)

mkdir -p ${DIR_WK}/AA_seq
mkdir -p ${DIR_WK}/AA_count
mkdir -p ${DIR_WK}/AA_weblogo
mkdir -p ${DIR_WK}/DNA_seq
mkdir -p ${DIR_WK}/DNA_count
mkdir -p ${DIR_WK}/DNA_weblogo

lab=$(awk -F "\t" -v line=${SLURM_ARRAY_TASK_ID} 'FNR == (line + 1) {print $1;}' ${DIR_WK}/info.txt)

Rscript ${DIR_WK}/code/compute.r \
    ${DIR_WK}/data/samstack \
    ${DIR_WK}/info.txt \
    $((SLURM_ARRAY_TASK_ID)) \
    ${DIR_WK}

####################
echo -e "DNA_weblogo\n"

function plotweblogo1 {
    echo ${1}
    cat ${DIR_WK}/DNA_seq/${1} | awk 'length($1) == 21' | \
        weblogo -A dna -U "bits" -c classic \
        > ${DIR_WK}/DNA_weblogo/${1}.bits.eps
}

export -f plotweblogo1

plotweblogo1 ${lab}.dna

function plotweblogo2 {
    echo ${1}
    cat ${DIR_WK}/DNA_seq/${1} | awk 'length($1) == 21' | \
        weblogo -A dna -U "probability" -c classic \
        > ${DIR_WK}/DNA_weblogo/${1}.probability.eps
}

export -f plotweblogo2

plotweblogo2 ${lab}.dna

####################
echo "AA_weblogo"

function AAweblogo1 {
    echo ${1}
    cat ${DIR_WK}/AA_seq/${1} | awk 'length($1) == 7' | \
        weblogo -A protein -U "bits" -c chemistry \
        > ${DIR_WK}/AA_weblogo/${1}.bits.eps

}

export -f AAweblogo1

AAweblogo1 ${lab}.aa

function AAweblogo2 {
    echo ${1}
    cat ${DIR_WK}/AA_seq/${1} | awk 'length($1) == 7' | \
        weblogo -A protein -U "probability" -c chemistry \
        > ${DIR_WK}/AA_weblogo/${1}.probability.eps

}

export -f AAweblogo2

AAweblogo2 ${lab}.aa

################################################################################
