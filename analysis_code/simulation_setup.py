import json
import os
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from parameters import RUN_RESULTS_FOLDER, SDD_TEMPLATE, PARAMETERS_TEMPLATE


def _read_template_header_and_chrom_sizes():
    """
    Reads SDD_TEMPLATE and returns (header_lines, chrom_sizes_bp):
    header_lines: list of lines (with trailing '\n'), from the start of the file up to and
                  including the '***EndOfHeader***;' line, copied verbatim.
    chrom_sizes_bp: list of chromosome sizes in bp (one per chromosome slot declared in the
                    'Chromosome sizes' field), in the same order as that field, rounded the
                    same way RadiSeq itself does (Mbp value * 1e6, rounded to the nearest bp).
    """
    with open(SDD_TEMPLATE, "r") as f:
        lines = f.readlines()

    header_lines = []
    chrom_sizes_bp = None
    for line in lines:
        header_lines.append(line)
        stripped = line.strip()
        if stripped.startswith("Chromosome sizes"):
            content = stripped.rstrip(";").rstrip(",")
            parts = [p.strip() for p in content.split(",")]
            n_chrom = int(parts[1])
            sizes_mbp = [float(x) for x in parts[2:2 + n_chrom]]
            chrom_sizes_bp = [round(size * 1e6) for size in sizes_mbp]
        if stripped.startswith("***EndOfHeader***"):
            break

    if chrom_sizes_bp is None:
        raise ValueError(f"Could not find a 'Chromosome sizes' field in template SDD file: {SDD_TEMPLATE}")
    if not header_lines or not header_lines[-1].strip().startswith("***EndOfHeader***"):
        raise ValueError(f"Could not find '***EndOfHeader***;' in template SDD file: {SDD_TEMPLATE}")

    return header_lines, chrom_sizes_bp


def _generate_exposure_lines(n_dsbs, n_ssbs, chrom_sizes_bp, n_chroms, dsb_threshold):
    """
    Returns the list of SDD data lines (each ending in '\n') for a single exposure containing
    n_dsbs DSBs and n_ssbs SSBs, placed at random locations across the genome described by
    chrom_sizes_bp (a numpy array of per-chromosome-slot lengths in bp, in 'Chromosome sizes'
    field order). Chromosomes are sampled with probability proportional to their length, so
    damages are uniformly distributed across the whole genome rather than uniformly across
    chromosomes. The exposure's first line is marked as a new exposure (leading '2'); every
    other line is marked '1' (see countExposuresSDD in fileio.cpp).
    """
    chrom_probabilities = chrom_sizes_bp / chrom_sizes_bp.sum()

    data_lines = []

    # DSBs: one line per DSB, with a backbone1 damage and a backbone2 damage no more than
    # dsb_threshold bp apart (as required by find_DSBs in induce_seq.cpp)
    if n_dsbs > 0:
        dsb_chrom_idx = np.random.choice(n_chroms, size=n_dsbs, p=chrom_probabilities)
        dsb_chrom_sizes = chrom_sizes_bp[dsb_chrom_idx]
        # Leave room for the second break, which lands up to dsb_threshold bp further along the chromosome
        max_location = np.maximum(dsb_chrom_sizes - dsb_threshold, 1)
        dsb_locations = (np.random.rand(n_dsbs) * max_location).astype(np.int64) + 1
        dsb_deltas = np.random.randint(0, dsb_threshold + 1, size=n_dsbs)

        for chrom_idx, location, delta in zip(dsb_chrom_idx, dsb_locations, dsb_deltas):
            chrom_n = chrom_idx + 1
            data_lines.append(
                f"1,0;\t 1,0;\t 1,{chrom_n},1,0;\t{location};\tx;\tz;\t1,1,1/4,{delta + 1},1;\n"
            )

    # SSBs: one line per SSB, with a single damage on a randomly chosen backbone
    if n_ssbs > 0:
        ssb_chrom_idx = np.random.choice(n_chroms, size=n_ssbs, p=chrom_probabilities)
        ssb_chrom_sizes = chrom_sizes_bp[ssb_chrom_idx]
        ssb_locations = (np.random.rand(n_ssbs) * ssb_chrom_sizes).astype(np.int64) + 1
        ssb_backbones = np.random.choice([1, 4], size=n_ssbs)

        for chrom_idx, location, backbone_n in zip(ssb_chrom_idx, ssb_locations, ssb_backbones):
            chrom_n = chrom_idx + 1
            data_lines.append(
                f"1,0;\t 1,0;\t 1,{chrom_n},1,0;\t{location};\tx;\tz;\t{backbone_n},1,1;\n"
            )

    # The very first data line of an exposure (and only that one) marks the start of a new exposure
    if data_lines:
        data_lines[0] = "2" + data_lines[0][1:]

    return data_lines


def generate_sdd_file(n_dsbs, n_ssbs, file_path, dsb_threshold):
    """
    sdd dna damage line format:
    [irradiation_number],0;    1,0;   1,[chrom_n],1,0;   [location_in_chrom]; x; z; [backbone_n],[offset],1;
    irradiation_number: first line should be 2, all others should be 1
    chrom_n: chromosome number.
    location_in_chrom: position of base pair in the chromosome
    backbone_n: 1 for 5' to 3' backbone, 4 for the other backbone
    offset: damage location relative to location_in_chrom. offset=1 means at that location, offset=3 means 2 base pairs further

    example of dsb line:
    2,0;	 1,0;	 1,6,1,0;	1920;	x;	z;	1,1,1/4,3,1;

    Generates an SDD file at file_path with a single exposure containing n_dsbs DSBs and
    n_ssbs SSBs. The new file's header is copied verbatim from SDD_TEMPLATE.
    """
    header_lines, chrom_sizes_bp = _read_template_header_and_chrom_sizes()
    n_chroms = len(chrom_sizes_bp)
    chrom_sizes_bp = np.array(chrom_sizes_bp, dtype=np.int64)

    data_lines = _generate_exposure_lines(n_dsbs, n_ssbs, chrom_sizes_bp, n_chroms, dsb_threshold)

    with open(file_path, "w") as f:
        f.writelines(header_lines)
        f.write("\n")
        f.writelines(data_lines)


def generate_multi_exposure_sdd_file(n_dsbs, n_ssbs, n_exposures, file_path, dsb_threshold):
    """
    Generates a single SDD file at file_path containing n_exposures separate exposures, each
    with n_dsbs DSBs and n_ssbs SSBs (see _generate_exposure_lines). RadiSeq processes each
    exposure in an SDD file as its own independent damaged cell, so with
    number_of_cells_in_sample/number_of_cells_to_sequence set to n_exposures, a single RadiSeq
    run against this file performs n_exposures independent trials.
    """
    header_lines, chrom_sizes_bp = _read_template_header_and_chrom_sizes()
    n_chroms = len(chrom_sizes_bp)
    chrom_sizes_bp = np.array(chrom_sizes_bp, dtype=np.int64)

    with open(file_path, "w") as f:
        f.writelines(header_lines)
        f.write("\n")
        for _ in range(n_exposures):
            f.writelines(_generate_exposure_lines(n_dsbs, n_ssbs, chrom_sizes_bp, n_chroms, dsb_threshold))

def _read_parameter_value(lines, key):
    """Returns the (string) value assigned to key in a RadiSeq parameters-file's lines, or None if absent."""
    for line in lines:
        before_comment = line.split("#", 1)[0]
        if "=" not in before_comment:
            continue
        line_key, line_value = before_comment.split("=", 1)
        if line_key.strip() == key:
            return line_value.strip()
    return None


def _set_parameter_lines(lines, updates):
    """
    Returns a new list of parameters-file lines with each key in updates set to its
    new value. Existing 'key = ...' lines are replaced in place (any trailing comment
    is preserved); keys not already present in lines are appended as new lines.
    """
    updates_remaining = dict(updates)
    new_lines = []
    for line in lines:
        before_comment, sep, comment = line.partition("#")
        if sep and "=" in before_comment:
            line_key = before_comment.split("=", 1)[0].strip()
            if line_key in updates_remaining:
                new_lines.append(f"{line_key} = {updates_remaining.pop(line_key)} #{comment}")
                continue
        elif "=" in before_comment:
            line_key = before_comment.split("=", 1)[0].strip()
            if line_key in updates_remaining:
                new_lines.append(f"{line_key} = {updates_remaining.pop(line_key)}\n")
                continue
        new_lines.append(line)

    for key, value in updates_remaining.items():
        new_lines.append(f"{key} = {value}\n")

    return new_lines


def set_up_simulation_folder(n_dsbs, n_ssbs, n_trials):
    """
    Sets up a folder (RUN_RESULTS_FOLDER/simulation_{n_dsbs}_{n_ssbs}) containing everything
    needed to run n_trials independent RadiSeq induce_seq trials with a single RadiSeq
    invocation:
      - ssd.sdd: a single SDD file containing n_trials separate exposures (generated via
        generate_multi_exposure_sdd_file), each with n_dsbs DSBs and n_ssbs SSBs, using the
        DSB threshold declared in PARAMETERS_TEMPLATE so the generated DSBs are actually
        detected as DSBs when RadiSeq is run on them.
      - parameters.txt: a copy of PARAMETERS_TEMPLATE, pointing sddFilePath at ssd.sdd with
        number_of_cells_in_sample and number_of_cells_to_sequence set to n_trials, so RadiSeq
        processes every exposure in ssd.sdd as its own independent cell/trial in one run.
      - description.json: records n_dsbs and n_ssbs.
    Returns the path to the created simulation folder.
    """
    sim_folder = os.path.join(RUN_RESULTS_FOLDER, f"simulation_{n_dsbs}_{n_ssbs}")
    os.makedirs(sim_folder, exist_ok=False)

    with open(PARAMETERS_TEMPLATE, "r") as f:
        template_lines = f.readlines()

    dsb_threshold_str = _read_parameter_value(template_lines, "DSB_threshold_in_bp")
    if dsb_threshold_str is None:
        raise ValueError(f"Could not find 'DSB_threshold_in_bp' in parameters template: {PARAMETERS_TEMPLATE}")
    dsb_threshold = int(dsb_threshold_str)

    sdd_filename = "ssd.sdd"
    generate_multi_exposure_sdd_file(
        n_dsbs, n_ssbs, n_trials, os.path.join(sim_folder, sdd_filename), dsb_threshold
    )

    parameters_lines = _set_parameter_lines(template_lines, {
        "sddFilePath": f"./{sdd_filename}",
        "number_of_cells_in_sample": n_trials,
        "number_of_cells_to_sequence": n_trials,
    })
    with open(os.path.join(sim_folder, "parameters.txt"), "w") as f:
        f.writelines(parameters_lines)

    description = {"n_dsbs": n_dsbs, "n_ssbs": n_ssbs, "time_of_simulation_setup": datetime.now().isoformat()}
    with open(os.path.join(sim_folder, "description.json"), "w") as f:
        json.dump(description, f, indent=4)

    return sim_folder

# if __name__ == "__main__":
    # density = 1/50000
    # for ssb_density in [1/250, 1/500, 1/1000, 1/10000, 1/100000, 1/50000, 1/1000000]:
    #     genome_length = 3.1e9
    #     n_dsbs = int(density*genome_length)
    #     n_ssbs = int(ssb_density*genome_length)
    #     set_up_simulation_folder(n_dsbs, n_ssbs, 3)
