"""
a script to vcf-merge files with common chr pos ref alt
it is adviced to place small ones in the begining and bigger to the end of the input list
python python-vcf-merge.py \
    -i \
    /home/PublicData/kazwg5_hgdp/kazwg5.hgdp.common_merged.vcf \
    /home/PublicData/jordelab/jordelab_2014_11_25/2010-Xing.TowardAmoreUniformSamplingOfHumanGeneticDiversity/vcf/jordelab.chr_gt.gold.vcf \
    /home/PublicData/hgdp/hgdp940_2014_10_15/plink2vcf/hgdp.ALL.gold.vcf \
    /home/PublicData/1000g_phase3/1000G.ALL.chr_gt.gold.vcf \
    -o /home/PublicData/merged/
"""

__VERSION__ = "0.1.0"

import argparse
import os

bgzip = "bgzip"
tabix = "tabix"
vcf_sort = "vcf-sort"
vcf_merge = "vcf-merge"


def timeit(method):
    import time
    def timed(*args, **kw):
        print("#     %r " % (method.__name__))
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print("# +++ %r %2.2f sec" % (method.__name__, te - ts), args[0])
        return result

    return timed

@timeit
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_vcf", default=[], required=False, nargs="+")
    parser.add_argument("-o", "--output_dir", required=True)
    args = parser.parse_args()
    vcf_list = args.input_vcf
    out_dir = args.output_dir
    python_vcf_merge_pipeline(vcf_list, out_dir)


def python_vcf_merge_pipeline(vcf_list, out_dir):
    for _, vcf in enumerate(vcf_list):
        map_file = generate_map_file(vcf, renew=False)
        map = load_map(map_file)
        if _ == 0:
            common_map = map
        else:
            common_map = common_map.intersection(map)
    names_list = []
    gz_list = []
    for _, vcf in enumerate(vcf_list):
        names_list += [get_name(vcf)]
        in_vcf_common = extract_map_from_vcf_to_vcf(vcf, common_map, out_dir)
        in_vcf_common_sorted = generate_output_filename(in_vcf_common, "sorted")
        in_vcf_common_sorted_gz = in_vcf_common_sorted + ".gz"
        cmd = [
            " ".join(["rm", "-f", in_vcf_common_sorted, in_vcf_common_sorted_gz]),
            " ".join([vcf_sort, in_vcf_common, ">", in_vcf_common_sorted]),
            " ".join([bgzip, in_vcf_common_sorted]),
            " ".join([tabix, in_vcf_common_sorted_gz]),
        ]
        sh = "; ".join(cmd)
        print (sh)
        os.system(sh)
        gz_list += [in_vcf_common_sorted_gz]
    out_vcf = "_".join(['merged'] + names_list) + '.vcf'
    cmd = " ".join([vcf_merge, " ".join(gz_list), ">", out_vcf])
    print (cmd)
    os.system(cmd)


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
def load_map(map_file) -> set:
    s = set()
    for line in open(map_file):
        arr = line.strip("\n").split("\t")
        if arr:
            s.add(tuple((arr[0], arr[1], arr[3], arr[4])))
    return s


@timeit
def get_name(file):
    return os.path.basename(file).strip().split(".")[0]


@timeit
def extract_map_from_vcf_to_vcf(vcf_file, map, out_dir) -> str:
    name = get_name(vcf_file)
    out_file = os.path.join(out_dir, name + '.common_map.vcf')
    with open(out_file, "w") as f:
        for line in open(vcf_file):
            if "CHROM" in line:
                f.write(line)
                continue
            elif line.startswith("#"):
                continue
            arr = line.strip("\n").split("\t")
            if tuple((arr[0], arr[1], arr[3], arr[4])) in map:
                f.write(line)
    return out_file


@timeit
def generate_output_filename(in_file, suffix="out") -> str:
    file_path, file_extension = os.path.splitext(in_file)
    return file_path + "." + suffix + file_extension


if __name__ == "__main__":
    main()
