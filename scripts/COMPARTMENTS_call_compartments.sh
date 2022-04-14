#!/bin/bash

WORKDIR=".."
DATADIR="${WORKDIR}/output/CONVERT/converted_hic/hic_all_intra"
OUTDIR="${WORKDIR}/output/COMPARTMENTS"

res=100000
condition=("as" "gfp")

chr_arr=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X)

for cond in "${condition[@]}"
do  
    hicfile="${DATADIR}/${cond}/${cond}_${res}_merged_iced_all_intra.hic"
    for chr in "${chr_arr[@]}"
    do
        echo "${chr}"
        eigen_file="${OUTDIR}/${cond}_all_intra_eigen_chr${chr}.txt"
        java -jar "${WORKDIR}/tools/juicer_tools_1.19.02.jar" eigenvector -p NONE "${hicfile}" "$chr" BP "${res}" "${eigen_file}" 
    done 

done
