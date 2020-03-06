"""
a script to vcf-merge files with common chr pos ref alt

it is adviced to place small ones in the begining and bigger to the end of the input list

python vcf-merge-for-pca-and-admixture.py \
    -i \
    /home/PublicData/kazwg5_hgdp/kazwg5.hgdp.common_merged.vcf \
    /home/PublicData/jordelab/jordelab_2014_11_25/2010-Xing.TowardAmoreUniformSamplingOfHumanGeneticDiversity/vcf/jordelab.chr_gt.gold.vcf \
    /home/PublicData/hgdp/hgdp940_2014_10_15/plink2vcf/hgdp.ALL.gold.vcf \
    /home/PublicData/1000g_phase3/1000G.ALL.chr_gt.gold.vcf \
    -o /home/PublicData/merged/

"""

__VERSION__ = "0.0.2"

import os
import argparse


bgzip = "/root/anaconda3/envs/merge/bin/bgzip"
tabix = "/root/anaconda3/envs/merge/bin/tabix"
vcf_sort = "/usr/local/bin/vcf-sort"
vcf_merge = "/usr/local/bin/vcf-merge"


def timeit(method):
    import time
    def timed(*args, **kw):
        print("#     %r " % (method.__name__))
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print("# +++ %r %2.2f sec" % (method.__name__, te - ts))
        return result

    return timed

@timeit
def merge_2_vcf(vcf_1, vcf_2, out_dir_prefix) -> str:
    name_1 = get_name(vcf_1)
    map_1 = generate_map_file(vcf_1, renew=False)
    ""
    name_2 = get_name(vcf_2)
    map_2 = generate_map_file(vcf_2, renew=False)
    ""
    out_name = f"{name_1}_{name_2}"
    out_dir = os.path.join(out_dir_prefix, out_name)
    mkdir(out_dir)
    out_prefix = os.path.join(out_dir, out_name)
    [abab] = compare_maps(map_1, map_2, out_prefix + "common_map.txt")
    out_vcf = out_prefix + ".common_merged.vcf"
    common_vcf_gz_1 = prepare_for_merging(vcf_1, abab, out_dir)
    common_vcf_gz_1 = prepare_for_merging(vcf_2, abab, out_dir)
    sh = " ".join([vcf_merge, common_vcf_gz_1, common_vcf_gz_1, ">", out_vcf])
    print (f"{sh}")
    os.system(sh)
    return out_vcf

@timeit
def generate_map_file(in_file, renew=False) -> str:
    map_file = generate_output_filename(in_file, "mapping")
    if not renew and os.path.isfile(map_file):
        return map_file
    else:
        with open(map_file, "w") as f:
            for line in open(in_file):
                if line.startswith("#"):
                    continue
                arr = line.strip("\n").split("\t")
                new_line = "\t".join(arr[:5]) + "\n"
                f.write(new_line)
        return map_file

@timeit
def get_name(file):
    return os.path.basename(file).strip().split(".")[0]

@timeit
def generate_output_filename(in_file, suffix="out") -> str:
    file_path, file_extension = os.path.splitext(in_file)
    return file_path + "." + suffix + file_extension

@timeit
def mkdir(dir_name):
    os.makedirs(dir_name, exist_ok=True)

@timeit
def compare_maps(map_1, map_2, out_map) -> list:
    num_1 = num_of_lines(map_1)
    num_2 = num_of_lines(map_2)
    if num_1 < num_2:
        map_1, map_2 = map_2, map_1
    d2 = load_map(map_2)
    name_abab = generate_output_filename(out_map, "file_abab")
    file_abab = open(name_abab, "w")
    rs_set = set()
    for line in open(map_1):
        arr = line.strip("\n").split("\t")
        if arr:
            if (arr[0], arr[1]) in d2:
                [c, p, r1, a1, b1] = [arr[i] for i in range(5)]
                [r2, a2, b2] = d2[(c, p)]
                r1 = r2 if r1 in ["---", "."] else r1
                r2 = r1 if r2 in ["---", "."] else r2
                r_flag = (r2 == r1) and (r1 not in rs_set)
                if a1 == a2 and b1 == b2 and r_flag:
                    file_abab.write("\t".join([c, p, r1, a1, b1]) + "\n")
                    rs_set.add(r1)
    file_abab.close()
    return [name_abab]

@timeit
def num_of_lines(in_file) -> int:
    i_ = -1
    for i_, _ in enumerate(open(in_file)):
        pass
    return i_ + 1

@timeit
def load_map(map_file) -> dict:
    d = dict()
    for line in open(map_file):
        arr = line.strip("\n").split("\t")
        if arr:
            d[(arr[0], arr[1])] = [arr[2], arr[3], arr[4]]
    return d

@timeit
def prepare_for_merging(in_vcf, abab, out_dir) -> str:
    in_vcf_common = extract_map_from_vcf_to_vcf(in_vcf, abab, out_dir)
    print("# in_vcf_common:", num_of_lines(in_vcf_common), in_vcf_common)
    ""
    in_vcf_common_sorted = generate_output_filename(in_vcf_common, "sorted")
    in_vcf_common_sorted_gz = in_vcf_common_sorted + ".gz"
    cmd = [
        " ".join(["rm", "-f", in_vcf_common_sorted, in_vcf_common_sorted_gz]),
        " ".join([vcf_sort, in_vcf_common, ">", in_vcf_common_sorted]),
        " ".join([bgzip, in_vcf_common_sorted]),
        " ".join([tabix, in_vcf_common_sorted_gz]),
    ]
    sh = "; ".join(cmd)
    print(f"{sh}")
    os.system(sh)
    return in_vcf_common_sorted_gz

@timeit
def extract_map_from_vcf_to_vcf(vcf_file, abab, out_dir) -> str:
    base = os.path.basename(vcf_file)
    out_file = os.path.join(out_dir, generate_output_filename(base, "common_chrom"))
    d_abab = load_map(abab)
    u = set()
    with open(out_file, "w") as f:
        for line in open(vcf_file):
            if "CHROM" in line:
                f.write(line)
                continue
            elif line.startswith("#"):
                continue
            arr = line.strip("\n").split("\t")
            c, p = arr[0], arr[1]
            if (c, p) in u:
                continue
            if (c, p) in d_abab:
                f.write(line)
                u.add((c, p))
    return out_file


################################################################################
################################################################################
################################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_vcf", default=[], required=False, nargs="+")
    parser.add_argument("-o", "--output_dir", required=True)
    args = parser.parse_args()
    vcf_list = args.input_vcf
    out_dir = args.output_dir

    new_vcf = vcf_list[0]
    for vcf in vcf_list[1:]:
        new_vcf = merge_2_vcf(new_vcf, vcf, out_dir)
        print(f"# merge_2_vcf   {new_vcf}   DONE")
