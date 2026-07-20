import json
import os

import numpy as np
import matplotlib.pyplot as plt

from parameters import RUN_RESULTS_FOLDER, FIGURES_FOLDER
from simulation_setup import _read_parameter_value


def _count_non_header_lines(csv_path):
    """Returns the number of data (non-header) lines in a *_sequenced_dsbs.csv file."""
    with open(csv_path, "r") as f:
        return sum(1 for _ in f) - 1


def _read_n_dsb_blunted_ends(simulation_data_csv_path):
    """Returns the n_dsb_blunted_ends value recorded in a *_simulation_data.csv file (one header line, one data row)."""
    with open(simulation_data_csv_path, "r") as f:
        f.readline()
        data_line = f.readline()
    return int(data_line.strip().split(",")[1])


def _output_dir_for(sim_dir, parameters_lines):
    output_rel = _read_parameter_value(parameters_lines, "output_directory_path") or "./output"
    return os.path.join(sim_dir, output_rel.lstrip("./"))


def graph_n_detected(n_ssbs):
    """
    For every simulation folder in RUN_RESULTS_FOLDER (as created by set_up_simulation_folder)
    whose description.json records the given n_ssbs, counts the number of sequenced DSBs
    recorded for each cell (the non-header line count of that cell's *_sequenced_dsbs.csv
    output file), then records the mean and standard deviation of that count across all cells
    in the folder, together with the folder's n_dsbs (from description.json). Plots these
    means, with the standard deviations as error bars, against n_dsbs. Folders without a
    description.json or without an output folder yet, or whose n_ssbs does not match, are skipped.
    """
    n_dsbs_list = []
    mean_list = []
    std_list = []

    for entry in sorted(os.listdir(RUN_RESULTS_FOLDER)):
        sim_dir = os.path.join(RUN_RESULTS_FOLDER, entry)
        description_path = os.path.join(sim_dir, "description.json")
        parameters_path = os.path.join(sim_dir, "parameters.txt")
        if not os.path.isdir(sim_dir) or not os.path.isfile(description_path) or not os.path.isfile(parameters_path):
            continue

        with open(description_path, "r") as f:
            description = json.load(f)
        if description["n_ssbs"] != n_ssbs:
            continue
        n_dsbs = description["n_dsbs"]

        with open(parameters_path, "r") as f:
            parameters_lines = f.readlines()
        output_dir = _output_dir_for(sim_dir, parameters_lines)
        if not os.path.isdir(output_dir):
            continue

        counts = [
            _count_non_header_lines(os.path.join(output_dir, filename))
            for filename in sorted(os.listdir(output_dir))
            if filename.endswith("_sequenced_dsbs.csv")
        ]
        if not counts:
            continue

        n_dsbs_list.append(n_dsbs)
        mean_list.append(np.mean(counts))
        std_list.append(np.std(counts))

    order = np.argsort(n_dsbs_list)
    n_dsbs_arr = np.array(n_dsbs_list)[order]
    mean_arr = np.array(mean_list)[order]
    std_arr = np.array(std_list)[order]

    plt.errorbar(n_dsbs_arr, mean_arr/n_dsbs_arr, yerr=std_arr/n_dsbs_arr, fmt="o-", capsize=4)
    plt.xscale("log")
    plt.xlabel("n_dsbs")
    plt.ylabel("Number of reads per DSB (mean ± std across cells)")
    plt.title("number of read per DSB vs. number of dsbs")
    plt.grid(True)

    os.makedirs(FIGURES_FOLDER, exist_ok=True)
    plt.savefig(os.path.join(FIGURES_FOLDER, "reads_per_dsb_vs_n_dsbs.png"))

    plt.show()

def graph_n_reads_vs_n_ssbs(n_dsbs):
    """
    For every simulation folder in RUN_RESULTS_FOLDER (as created by set_up_simulation_folder)
    whose description.json records the given n_dsbs, computes, for each cell, the number of
    sequenced DSBs (the non-header line count of that cell's *_sequenced_dsbs.csv output file)
    divided by that same cell's actual number of blunted DSB ends (n_dsb_blunted_ends, recorded
    in that cell's *_simulation_data.csv file) rather than the nominal n_dsbs, since the two can
    differ (e.g. nearby DSBs merging into a single blunted end). Records the mean and standard
    deviation of that per-cell ratio across all cells in the folder, together with the mean and
    standard deviation of n_dsb_blunted_ends itself, and the folder's n_ssbs (from
    description.json). Plots the reads-per-blunted-end means (with std as error bars) against
    n_ssbs, and separately plots the n_dsb_blunted_ends means (with std as error bars) against
    n_ssbs. Folders without a description.json or without an output folder yet, or whose n_dsbs
    does not match, are skipped.
    """
    n_ssbs_list = []
    mean_list = []
    std_list = []
    blunted_ends_mean_list = []
    blunted_ends_std_list = []

    for entry in sorted(os.listdir(RUN_RESULTS_FOLDER)):
        sim_dir = os.path.join(RUN_RESULTS_FOLDER, entry)
        description_path = os.path.join(sim_dir, "description.json")
        parameters_path = os.path.join(sim_dir, "parameters.txt")
        if not os.path.isdir(sim_dir) or not os.path.isfile(description_path) or not os.path.isfile(parameters_path):
            continue

        with open(description_path, "r") as f:
            description = json.load(f)
        if description["n_dsbs"] != n_dsbs:
            continue
        n_ssbs = description["n_ssbs"]

        with open(parameters_path, "r") as f:
            parameters_lines = f.readlines()
        output_dir = _output_dir_for(sim_dir, parameters_lines)
        if not os.path.isdir(output_dir):
            continue

        ratios = []
        blunted_ends = []
        for filename in sorted(os.listdir(output_dir)):
            if not filename.endswith("_sequenced_dsbs.csv"):
                continue
            simulation_data_path = os.path.join(
                output_dir, filename[: -len("_sequenced_dsbs.csv")] + "_simulation_data.csv"
            )
            if not os.path.isfile(simulation_data_path):
                continue
            count = _count_non_header_lines(os.path.join(output_dir, filename))
            n_dsb_blunted_ends = _read_n_dsb_blunted_ends(simulation_data_path)
            ratios.append(count / n_dsb_blunted_ends)
            blunted_ends.append(n_dsb_blunted_ends)
        if not ratios:
            continue

        n_ssbs_list.append(n_ssbs)
        mean_list.append(np.mean(ratios))
        std_list.append(np.std(ratios))
        blunted_ends_mean_list.append(np.mean(blunted_ends))
        blunted_ends_std_list.append(np.std(blunted_ends))

    order = np.argsort(n_ssbs_list)
    n_ssbs_arr = np.array(n_ssbs_list)[order]
    mean_arr = np.array(mean_list)[order]
    std_arr = np.array(std_list)[order]
    blunted_ends_mean_arr = np.array(blunted_ends_mean_list)[order]
    blunted_ends_std_arr = np.array(blunted_ends_std_list)[order]

    os.makedirs(FIGURES_FOLDER, exist_ok=True)

    plt.figure()
    plt.errorbar(n_ssbs_arr, mean_arr, yerr=std_arr, fmt="o-", capsize=4)
    plt.xscale("symlog")
    plt.xlabel("n_ssbs")
    plt.ylabel("Number of reads per DSB (mean ± std across cells)")
    plt.title("number of read per DSB vs. number of ssbs")
    plt.grid(True)
    plt.savefig(os.path.join(FIGURES_FOLDER, "reads_per_dsb_vs_n_ssbs.png"))
    plt.show()

    plt.figure()
    plt.errorbar(n_ssbs_arr, blunted_ends_mean_arr, yerr=blunted_ends_std_arr, fmt="o-", capsize=4)
    plt.xscale("symlog")
    plt.xlabel("n_ssbs")
    plt.ylabel("Number of DSB blunted ends (mean ± std across cells)")
    plt.title("number of DSBs vs. number of SSBs")
    plt.grid(True)
    plt.savefig(os.path.join(FIGURES_FOLDER, "n_dsb_blunted_ends_vs_n_ssbs.png"))
    plt.show()

if __name__ == "__main__" : graph_n_reads_vs_n_ssbs(62000)