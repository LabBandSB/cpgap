from collections import defaultdict

"""
pop_file = "pop.kaz.1kgp.tab.txt"
vcf_file = "merged.1000g.kaz24567.attempt3.filtered_dup_rs.filtered_dup_chrom_pos.modified_gt_alt.autosomes.vcf"
split -l 40 run_vcftools_fst.sh
"""

__fst_pop_dir__ = "fst_pop_dir"
__fst_index_dir__ = "fst_index_dir"
__fst_run_cmd__ = "fst_index_run_vcftools.sh"
__VCFTOOLS_PATH__ = "# export PATH=$PATH:/home/Pipeline/vcftools/vcftools_0.1.13/bin/"


def is_file_exists(file_name) -> bool:
    import os
    return os.path.isfile(file_name)


def mkdir(dir_name) -> None:
    import os
    try:
        os.makedirs(dir_name)
    except:
        pass


def parse_arguments() -> dict:  # dict of string_key to string_value
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_pop_file", "-p", default="")
    parser.add_argument("--input_vcf_file", "-v", default="")
    parser.add_argument("--output_pop_dir", default=__fst_pop_dir__)
    parser.add_argument("--output_fst_index_dir", default=__fst_index_dir__)
    parser.add_argument("--output_bash_script", default=__fst_run_cmd__)
    args = parser.parse_args()
    d = {
        "pop_file": args.input_pop_file,
        "vcf_file": args.input_vcf_file,
        "pop_dir": args.output_pop_dir,
        "fst_dir": args.output_fst_index_dir,
        "out_bash": args.output_bash_script,
    }
    return d


def get_header_list(vcf_file) -> list:  # list of string_elemens
    for line in open(vcf_file, "r"):
        if "CHROM" in line:
            header_list = line.strip().split("\t")
            return header_list


def main() -> None:
    d = parse_arguments()
    assert is_file_exists(d["pop_file"])
    assert is_file_exists(d["vcf_file"])
    mkdir(d["pop_dir"])
    mkdir(d["fst_dir"])
    d["pop_dict"] = get_pop_dict_with_samples(d)
    extract_pop_to_single_files(d)
    prepare_vcftools_bash_script(d)


def get_pop_dict_with_samples(d) -> dict:  # defaultdict of string_key to list_of_strings_value
    pop_file = d["pop_file"]
    vcf_file = d["vcf_file"]
    header_list = get_header_list(vcf_file)
    res = defaultdict(list)
    for line in open(pop_file, "r"):
        if line:
            arr = line.strip().split("\t")
            sample, pop = arr[0], arr[2]
            if sample in header_list:
                res[pop].append(sample)
    return res


def extract_pop_to_single_files(d) -> None:  #
    pop_dict = d["pop_dict"]
    pop_dir = d["pop_dir"]
    for pop, samples_list in pop_dict.items():
        out_file = "{:s}/{:s}.txt".format(pop_dir, pop)
        with open(out_file, "w") as f:
            for sample in samples_list:
                new_line = sample + "\n"
                f.write(new_line)


def prepare_vcftools_bash_script(d) -> None:  #
    from itertools import combinations
    pop_dir = d["pop_dir"]
    out_dir = d["fst_dir"]
    pop_list = sorted(d["pop_dict"])
    vcf_file = d["vcf_file"]
    out_file = d["out_bash"]
    c = 0
    with open(out_file, "w") as f:
        for p, q in combinations(pop_list, 2):
            r = {
                "v": vcf_file,
                "p": "{:s}/{:s}.txt".format(pop_dir, p),
                "q": "{:s}/{:s}.txt".format(pop_dir, q),
                "o": "{:s}/fst_{:s}_vs_{:s}_.txt".format(out_dir, p, q)
            }
            new_line = "vcftools --vcf {v} --weir-fst-pop {p} --weir-fst-pop {q} --out {o} > {o}.log\n".format(**r)
            f.write(new_line)
            c += 1
    print("# {:s} has {:d} lines".format(out_file, c))
    print("# " + __VCFTOOLS_PATH__)
    print("# # split -l ({:d}/cores)+2 {:s}".format(out_file, c))
    print("# # cat {:s} | parallel --pipe bash".format(out_file))


if __name__ == "__main__":
    main()
