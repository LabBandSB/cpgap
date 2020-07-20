__VERSION__ = "0.0.1_raw"

import argparse


def run_admixture_from_plink_to_admixture(vcf_file):
    plink_file = f'{vcf_file}.plink'
    prune_file = f'{vcf_file}.pruning'
    admixture_file = f'{vcf_file}.admixture'

    print('\n')
    print('*' * 40)
    print(f'input file: {vcf_file}')
    print('*' * 40)
    print(f'plink_file: {plink_file}')
    print(f'bin_file  : {prune_file}')
    print(f'gcta_file : {admixture_file}')
    print('*' * 40)
    print('\n')

    print('ulimit -n 300000')
    print(f"plink --noweb --indep-pairwise 50 10 0.1 --file {plink_file} --out {prune_file}")
    print(f"plink --noweb --extract {prune_file}.prune.in --make-bed --file {plink_file} --out {admixture_file}")
    print(f"for K in 2 3 4 5 6 7 8 9 10 11 12; do admixture --cv {admixture_file}.bed \$K > log_\$K.prune.out ; done")
    print('echo -e "\\n\\nPCA_calculations: -*- DONE -*-" ||')
    print('echo -e "\\n\\nPCA_calculations: -*- FAILED -*-"')
    print('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='run admixture pipeline')
    parser.add_argument('-i', '--input_vcf', required=True)
    args = parser.parse_args()
    run_admixture_from_plink_to_admixture(args.input_vcf)
