#!/bin/bash

cd ./tools

chr_arr=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chrX")
size_arr=(282763074 266435125 177699992 184226339 173707219 147991367 145729302 133307652 122095297 112626471 90463843 52716770 114033958 115493446 111246239 90668790 90843779 88201929 62275575 56205956 159970021)

for i in ${!chr_arr[*]}; do
    echo "element $i is ${chr_arr[$i]} and ${size_arr[$i]}"

    ./OnTAD ../output/CONVERT/converted_hic/hic_all_intra/as_10000_merged_iced_all_intra.hic -bedout "${chr_arr[$i]}" "${size_arr[$i]}" 10000 -o "../output/TADS/as_10000_merged_iced_all_intra_${chr_arr[$i]}.ontad"

    ./OnTAD ../output/CONVERT/converted_hic/hic_all_intra/gfp_10000_merged_iced_all_intra.hic -bedout "${chr_arr[$i]}" "${size_arr[$i]}" 10000 -o "../output/TADS/gfp_10000_merged_iced_all_intra_${chr_arr[$i]}.ontad"
done

cd ../
