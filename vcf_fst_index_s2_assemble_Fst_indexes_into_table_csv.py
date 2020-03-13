import os
import re
from collections import defaultdict

__fst_pop_dir__ = "fst_pop_dir"
__fst_index_dir__ = "fst_index_dir"
__fst_run_cmd__ = "fst_index_run_vcftools.sh"


def parse_arguments() -> dict:  # dict of string_key to string_value
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_fst_index_dir", default=__fst_index_dir__)
    args = parser.parse_args()
    d = {
        "fst_dir": args.output_fst_index_dir
    }
    return d


def get_log_files_list(d):
    fst_out_dir = d["fst_dir"]
    files_list = os.listdir(fst_out_dir)
    result_list = [fst_out_dir + '/' + file for file in files_list if file.endswith("txt.log")]
    return result_list


def main() -> None:
    d = parse_arguments()
    assert os.path.exists(d["fst_dir"])
    d['log_files_list'] = get_log_files_list(d)
    d['log_info_dict'] = get_log_info_dict(d)
    write_fst_index_table(d, "mean")
    write_fst_index_table(d, "weighted")
    write_fst_index_report(d)
    print("done")


def get_log_info_dict(d):
    log_files_list = d['log_files_list']
    log_info_dict = defaultdict(tuple)
    for file in log_files_list:
        pop1, pop2 = file.split("/")[-1].replace("fst_", "").replace("_.txt.log", "").split("_vs_")
        target, total, mean, weighted, run_time = "0", "0", "0", "0", "0"
        for line in open(file, "r"):
            line = line.strip()
            if line.startswith("After filtering, kept") and line.endswith("Individuals"):
                target, total = extract_digits_from_str(line)[:2]
            elif line.startswith("Weir and Cockerham mean Fst estimate:"):
                mean = extract_digits_from_str(line)[0]
            elif line.startswith("Weir and Cockerham weighted Fst estimate:"):
                weighted = extract_digits_from_str(line)[0]
            elif line.startswith("Run Time = "):
                run_time = extract_digits_from_str(line)[0]
        log_info_dict[(pop1, pop2)] = (pop1, pop2, target, total, mean, weighted, run_time)
    return log_info_dict


def extract_digits_from_str(s):
    result = re.findall(r"\d+\.*\d*", s)
    # print(s, ">>>", result)
    if result:
        return result
    else:
        return ["nan", "nan"]


def write_fst_index_table(d, target):
    log_info_dict = d['log_info_dict']
    assert target == "mean" or target == "weighted"
    index = 4 if target == "mean" else 5
    pop_list = sorted(set([j for i in log_info_dict for j in i]))
    with open("Fst_index_table.{:s}.csv".format(target), "w") as f:
        f.write("\t".join(pop_list) + "\n")
        for i in pop_list:
            arr = [i]
            for j in pop_list:
                if i == j:
                    arr.append("0")
                elif log_info_dict[(i, j)]:
                    arr.append(log_info_dict[(i, j)][index])
                elif log_info_dict[(j, i)]:
                    arr.append(log_info_dict[(j, i)][index])
                else:
                    print(i, j, "???")
            f.write("\t".join(arr) + "\n")


def write_fst_index_report(d):
    log_info_dict = d['log_info_dict']
    with open("Fst_index_table.report.csv", "w") as f:
        for pair, value in log_info_dict.items():
            f.write("\t".join(list(pair) + list(value)) + "\n")


if __name__ == "__main__":
    main()
