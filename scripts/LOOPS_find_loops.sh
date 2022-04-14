#!/bin/bash 

#SBATCH --job-name=job_06Apr2022
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hyeyeon_hwang@brown.edu
#SBATCH --mem=200G
#SBATCH --time=20:00:00
#SBATCH --output=/users/hhwang16/tapinos_lab/egr2as-hic/slurm/v1_oct20_%j.out

#WORKDIR="/Users/hyeyeonhwang/tapinos_lab/research_projects/HiC-Analysis-of-Egr2-AS-RNA-v1"
WORKDIR=".."
DATADIR="${WORKDIR}/data/hicPro/10000"

SCRIPTDIR="${WORKDIR}/scripts"
script_features="${SCRIPTDIR}/LOOPS_generate_features.R"
script_sigloops="${SCRIPTDIR}/LOOPS_call_significant_loops.R"
script_diffloops="${SCRIPTDIR}/LOOPS_call_differential_loops.R"
script_loopsets="${SCRIPTDIR}/LOOPS_loops_sets_criteria.R"

OUTDIR="${WORKDIR}/output/LOOPS"
outdir_features="${OUTDIR}/features"
outdir_sig_analysis="${OUTDIR}/significant_loops"
outdir_diff_analysis="${OUTDIR}/differential_loops"
outdir_loopsets="${OUTDIR}/loops_sets"


module load R/4.1.0
module load gcc/10.2 pcre2/10.35 intel/2020.2 texlive/2018

res=10000

# specified: 10G,  1:00:00
# used: 3.05 GB, 00:19:12
# do: 10G, 1:00:00
# Already ran
start_features=$(date +'%d%b%Y %H:%M:%S')
echo "********** START | Generate features | ${start_features}"
Rscript --vanilla "${script_features}" --outdir "${outdir_features}" --res "${res}"
end_features=$(date +'%d%b%Y %H:%M:%S')
echo "********** END | Generate features | ${end_features}"

# rn6_10000_GATC_GANTC_features_bintolen.txt.gz
# size: 3540116 = 3.4M
# features file has ending *features_bintolen.txt.gz
# get features file extension in ${outdir_features}
features_path=$(find "${outdir_features}" -type f -name "*_features_bintolen.txt.gz")
echo "$features_path"

# specified 100G, 10:00:00
# used: 14.07 GB, 00:39:12 
# do: 30GB, 2:00:00
start_sigloops=$(date +'%d%b%Y %H:%M:%S')
echo "********** START | Call significant loops | ${start_sigloops}"
Rscript --vanilla "${script_sigloops}" --datadir "${DATADIR}" --features "${features_path}" --outdir "${outdir_sig_analysis}" --res "${res}" --distrange "50Kb-2Mb" --dmin 50000 --sampsize 0.01 --qval 0.05 
end_sigloops=$(date +'%d%b%Y %H:%M:%S')
echo "********** END | Call significant loops | ${end_sigloops}"

# total sigindices has ending *sigindices.txt.gz
# get sigindices file extension in ${outdir_sig_analysis}
sigindices=$(find "${outdir_sig_analysis}" -type f -name "*_sigindices.txt.gz")
echo "$sigindices"

# specified 100G, 10:00:00
start_diffloops=$(date +'%d%b%Y %H:%M:%S')
echo "********** START | Call differential loops | ${start_diffloops}"
Rscript --vanilla "${script_diffloops}" --sigindices "${sigindices}" --sigdir "${outdir_sig_analysis}" --outdir "${outdir_diff_analysis}" --res "${res}" --distrange "50Kb-2Mb" --dmin 50000 
end_diffloops=$(date +'%d%b%Y %H:%M:%S')
echo "********** END | Call differential loops | ${end_diffloops}"

start_loopsets=$(date +'%d%b%Y %H:%M:%S')
echo "********** START | Make loops sets | ${start_loopsets}"
Rscript --vanilla "${script_loopsets}" --diffdir "${outdir_diff_analysis}" --sigloops "${sigindices}" --outdir "${outdir_loopsets}" --res "${res}" --distrange "50Kb-2Mb" --pval 0.05 --log2fc 1
end_loopsets=$(date +'%d%b%Y %H:%M:%S')
echo "********** END | Make loops sets | ${end_loopsets}"

