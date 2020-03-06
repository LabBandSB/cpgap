__VERSION__ = "0.0.1_raw"

import argparse


def run_gcta_from_vcf_to_plink_to_pca(vcf_file):
	plink_file = f'{vcf_file}.plink'
	bin_file   = f'{vcf_file}.bin'
	gcta_file  = f'{vcf_file}.gcta'
	pca_file   = f'{vcf_file}.pca'

	print('\n')
	print('*' * 40)
	print(f'input file: {vcf_file}')
	print('*' * 40)
	print(f'plink_file: {plink_file}')
	print(f'bin_file  : {bin_file}')
	print(f'gcta_file : {gcta_file}')
	print(f'pca_file  : {pca_file}')
	print('*' * 40)
	print('\n')

	print('ulimit -n 300000')
	print(f'vcftools --vcf {vcf_file} --plink --out {plink_file} &&')
	print(f'plink --file {plink_file} --out {bin_file} --allow-no-sex --make-bed --noweb &&')
	print(f'gcta64 --bfile {bin_file} --autosome --maf 0.05 --make-grm --out {gcta_file} --thread-num 8 &&')
	print(f'gcta64 --grm {gcta_file} --pca 20 --out {gcta_file} &&')
	print(f'$plink2evec --eigenval {gcta_file}.eigenval --eigenvec {gcta_file}.eigenvec --out {pca_file} &&')
	print('echo -e "\\n\\nPCA_calculations: -*- DONE -*-" ||')
	print('echo -e "\\n\\nPCA_calculations: -*- FAILED -*-"')
	print('\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='run gcta pca pipeline')
	parser.add_argument('-i', '--input_vcf', required = True)
	args = parser.parse_args()
	run_gcta_from_vcf_to_plink_to_pca(args.input_vcf)
