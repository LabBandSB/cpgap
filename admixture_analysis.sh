#!/bin/bash
#
# author: -*- askhat -*-
# date  : -*- 201707171211 -*-
# env   : -*- qsymphony -*-
# 
# __VERSION__ = "0.0.1_raw"


function run_admixture_from_vcf_to_admixture {
    ulimit -n 300000
    #
    plink="/home/Pipeline/plink/plink-1.07-x86_64/plink"
    admixture="/home/Pipeline/admixture/admixture_linux-1.3.0/admixture"
    #
    vcf_file=${1}
    #
    plink_file=${vcf_file}.plink
    prune_file=${vcf_file}.pruning
    admixture_file=${vcf_file}.admixture
    #
    echo ""
    echo "      vcf_file=${vcf_file}"
    echo "    plink_file=${plink_file}"
    echo "    prune_file=${prune_file}"
    echo "admixture_file=${admixture_file}"
    echo ""
    #
    echo "${plink} --noweb --indep-pairwise 50 10 0.1 --file ${plink_file} --out ${prune_file}" &&
          ${plink} --noweb --indep-pairwise 50 10 0.1 --file ${plink_file} --out ${prune_file} &&
    #
    echo "${plink} --noweb --extract ${prune_file}.prune.in --make-bed --file ${plink_file} --out ${admixture_file}" &&
          ${plink} --noweb --extract ${prune_file}.prune.in --make-bed --file ${plink_file} --out ${admixture_file} &&
    #
    echo "for K in 2 3 4 5 6 7 8 9 10 11 12; do ${admixture} --cv ${admixture_file}.bed \$K > log_\$K.prune.out ; done" &&
          for K in 2 3 4 5 6 7 8 9 10 11 12; do ${admixture} --cv ${admixture_file}.bed $K > log_$K.prune.out & done &&
    #
    echo -e "\n\nADMIXTUERE: -*- DONE -*-" ||
    echo -e "\n\nADMIXTUERE: -*- FAILED -*-"
}

# run_admixture_from_vcf_to_admixture
