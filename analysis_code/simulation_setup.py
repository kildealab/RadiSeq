import json
import os
import re
import shutil
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from parameters import (
    RUN_RESULTS_FOLDER,
    SDD_TEMPLATE,
    PARAMETERS_TEMPLATE,
    SDD_FILES,
    FINAL_SDDS_RUN_RESULTS_FOLDER,
    FINAL_SDDS_PARAMETERS_TEMPLATE,
)

_SDD_OUTPUT_NUMBER_RE = re.compile(r"SDDOutput_(\d+)\.txt$")


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


def _sdd_numbers_in_folder(folder):
    """Returns {sdd_number: file_path} for every 'SDDOutput_<N>.txt' file found directly in folder."""
    numbers = {}
    for name in os.listdir(folder):
        match = _SDD_OUTPUT_NUMBER_RE.match(name)
        if match:
            numbers[int(match.group(1))] = os.path.join(folder, name)
    return numbers


def _write_final_sdds_parameters(template_lines, dest_folder, sdd_file_paths, merge, particle_names, genome_fasta_path):
    """
    Creates dest_folder and writes a parameters.txt in it (based on template_lines), with
    sddFilePath set to sdd_file_paths (a comma-separated list of absolute paths), merge flags set
    according to merge/len(sdd_file_paths), primary_particles_simulated set to particle_names, and
    induce_seq_genome_fasta_path set to genome_fasta_path (an absolute path, since the template's
    relative path assumes a shallower folder depth than these generated folders sit at).

    If dest_folder already has an output folder (its path taken from the new parameters.txt's own
    output_directory_path, so this stays correct even if the template changes it), that folder --
    along with any results from a previous run -- is removed, so the next run starts clean instead
    of mixing with stale output.
    """
    os.makedirs(dest_folder, exist_ok=True)
    parameters_lines = _set_parameter_lines(template_lines, {
        "sddFilePath": ", ".join(sdd_file_paths),
        "merge_damages_from_multiple_particles": "true" if merge else "false",
        "number_of_particles_to_merge": len(sdd_file_paths),
        "primary_particles_simulated": ",".join(particle_names),
        "induce_seq_genome_fasta_path": genome_fasta_path,
    })
    with open(os.path.join(dest_folder, "parameters.txt"), "w") as f:
        f.writelines(parameters_lines)

    output_rel = _read_parameter_value(parameters_lines, "output_directory_path") or "./output"
    output_dir = os.path.join(dest_folder, output_rel.lstrip("./"))
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)


def set_up_final_sdds_simulations():
    """
    Mirrors Final_SDDs' Neutron/Photon folder structure (skipping any folder with 'nico' in its
    name, e.g. '1MeV_outer_nico') into FINAL_SDDS_RUN_RESULTS_FOLDER, creating one simulation
    folder per SDD number found, each containing a parameters.txt based on
    FINAL_SDDS_PARAMETERS_TEMPLATE:
      - Neutron/<energy>/<dose>/<number>/parameters.txt: combines that number's proton and electron
        SDD files (merge_damages_from_multiple_particles=true, number_of_particles_to_merge=2).
        Only SDD numbers present in *both* the proton and electron folders are used, since a
        combined-damage simulation needs both.
      - Photon/<energy>/<dose>/<number>/parameters.txt: uses that number's single SDD file directly
        (merge_damages_from_multiple_particles=false), since Photon folders aren't split by particle.
    Re-running this (e.g. after Final_SDDs changes) overwrites each folder's parameters.txt and
    deletes any output folder already there (see _write_final_sdds_parameters), so any results from
    a previous run of that simulation are removed along with it.
    Returns the number of simulation folders created.
    """
    with open(FINAL_SDDS_PARAMETERS_TEMPLATE, "r") as f:
        template_lines = f.readlines()

    genome_fasta_path = os.path.abspath(os.path.join("radiSeqData", "induce_seq_human_genome.fa"))

    n_created = 0

    # ----- Neutron: combine the proton + electron SDDs sharing the same number -----
    neutron_root = os.path.join(SDD_FILES, "Neutron")
    for energy_name in sorted(os.listdir(neutron_root)):
        energy_folder = os.path.join(neutron_root, energy_name)
        if not os.path.isdir(energy_folder) or "nico" in energy_name.lower():
            continue

        for dose_name in sorted(os.listdir(energy_folder)):
            dose_folder = os.path.join(energy_folder, dose_name)
            if not os.path.isdir(dose_folder):
                continue
            electron_folder = os.path.join(dose_folder, "electron")
            proton_folder = os.path.join(dose_folder, "proton")
            if not (os.path.isdir(electron_folder) and os.path.isdir(proton_folder)):
                continue

            electron_files = _sdd_numbers_in_folder(electron_folder)
            proton_files = _sdd_numbers_in_folder(proton_folder)
            shared_numbers = sorted(set(electron_files) & set(proton_files))

            for number in shared_numbers:
                dest_folder = os.path.join(FINAL_SDDS_RUN_RESULTS_FOLDER, "Neutron", energy_name, dose_name, str(number))
                _write_final_sdds_parameters(
                    template_lines, dest_folder,
                    [os.path.abspath(proton_files[number]), os.path.abspath(electron_files[number])],
                    merge=True, particle_names=["proton", "electron"], genome_fasta_path=genome_fasta_path,
                )
                n_created += 1

    # ----- Photon: a single SDD file per simulation, no particle split -----
    photon_root = os.path.join(SDD_FILES, "Photon")
    for energy_name in sorted(os.listdir(photon_root)):
        energy_folder = os.path.join(photon_root, energy_name)
        if not os.path.isdir(energy_folder) or "nico" in energy_name.lower():
            continue

        for dose_name in sorted(os.listdir(energy_folder)):
            dose_folder = os.path.join(energy_folder, dose_name)
            if not os.path.isdir(dose_folder):
                continue

            sdd_files = _sdd_numbers_in_folder(dose_folder)
            for number in sorted(sdd_files):
                dest_folder = os.path.join(FINAL_SDDS_RUN_RESULTS_FOLDER, "Photon", energy_name, dose_name, str(number))
                _write_final_sdds_parameters(
                    template_lines, dest_folder,
                    [os.path.abspath(sdd_files[number])],
                    merge=False, particle_names=["photon"], genome_fasta_path=genome_fasta_path,
                )
                n_created += 1

    return n_created


# GRCh37 primary-assembly chromosome lengths in bp, from the NCBI GRC data page:
# https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
GRCH37_CHROM_LENGTHS_BP = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260,
    "6": 171115067, "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747,
    "11": 135006516, "12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392,
    "16": 90354753, "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
}


def _read_chrom_sizes_from_sdd_header(sdd_path):
    """
    Reads sdd_path's header and returns the chromosome sizes declared in its 'Chromosome
    sizes' field, in bp (Mbp value * 1e6, rounded to the nearest bp, matching how RadiSeq
    itself interprets this field -- see NGSsdd::set_chrom_size_bp). Returns None if the file
    has no such field before '***EndOfHeader***;' (or no such line at all).
    """
    with open(sdd_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("Chromosome sizes"):
                content = stripped.rstrip(";").rstrip(",")
                parts = [p.strip() for p in content.split(",")]
                n_chrom = int(parts[1])
                sizes_mbp = [float(x) for x in parts[2:2 + n_chrom]]
                return [round(size * 1e6) for size in sizes_mbp]
            if stripped.startswith("***EndOfHeader***"):
                break
    return None


def _slot_reference_chrom_names(n_chrom):
    """
    Returns the GRCH37_CHROM_LENGTHS_BP key for each of the n_chrom chromosome slots declared
    in an SDD file's 'Chromosome sizes' field, assuming slots are listed in the order used
    throughout this project's SDD files: autosomes 1..22 (copy 1), then autosomes 1..22 again
    (copy 2), then Y, then X.
    """
    autosomes = [str(i) for i in range(1, 23)]
    names = autosomes + autosomes + ["Y", "X"]
    if n_chrom != len(names):
        raise ValueError(f"Expected {len(names)} chromosome slots (autosomes x2 + Y + X), got {n_chrom}")
    return names


def _sdd_file_paths_referenced_by(results_folder):
    """
    Returns the set of SDD file paths used by simulation folders under results_folder: every
    path assigned to 'sddFilePath' in any parameters.txt found in the tree (comma-separated
    lists are split into individual paths), plus any '*.sdd' file found directly in the tree.
    """
    sdd_paths = set()
    for root, _dirs, files in os.walk(results_folder):
        for name in files:
            if name == "parameters.txt":
                with open(os.path.join(root, name), "r") as f:
                    lines = f.readlines()
                value = _read_parameter_value(lines, "sddFilePath")
                if value:
                    sdd_paths.update(p.strip() for p in value.split(",") if p.strip())
            elif name.endswith(".sdd"):
                sdd_paths.add(os.path.join(root, name))
    return sdd_paths


def check_sdd_files_for_out_of_bounds_damage(results_folder=FINAL_SDDS_RUN_RESULTS_FOLDER):
    """
    Goes through every SDD file referenced from results_folder (see
    _sdd_file_paths_referenced_by), reads the chromosome sizes declared in its own header, and
    checks every DNA damage line's position-in-chromosome (SDD field 4) against the real
    GRCh37 length of the chromosome that slot corresponds to (GRCH37_CHROM_LENGTHS_BP, ordered
    per _slot_reference_chrom_names). Prints one entry per SDD file that contains at least one
    damage whose position exceeds the real chromosome length, summarizing -- per affected
    chromosome slot -- how many such damages were found and the worst overshoot in bp.
    """
    sdd_paths = sorted(_sdd_file_paths_referenced_by(results_folder))
    print(f"Checking {len(sdd_paths)} SDD file(s) referenced from {results_folder}...")

    n_files_with_overflow = 0
    for sdd_path in sdd_paths:
        chrom_sizes_bp = _read_chrom_sizes_from_sdd_header(sdd_path)
        if chrom_sizes_bp is None:
            print(f"  {sdd_path}: no 'Chromosome sizes' header field found, skipping")
            continue
        chrom_names = _slot_reference_chrom_names(len(chrom_sizes_bp))

        overflow = {}  # slot_idx -> (count_of_damages_beyond_real_length, max_overshoot_bp)
        with open(sdd_path, "r") as f:
            past_header = False
            for line in f:
                if not past_header:
                    if line.strip().startswith("***EndOfHeader***"):
                        past_header = True
                    continue
                fields = line.split(";")
                if len(fields) < 4:
                    continue
                chrom_id = int(fields[2].split(",")[1].strip())
                position_in_chrom = int(fields[3].strip())
                slot_idx = chrom_id - 1
                real_length = GRCH37_CHROM_LENGTHS_BP[chrom_names[slot_idx]]
                if position_in_chrom > real_length:
                    count, max_overshoot = overflow.get(slot_idx, (0, 0))
                    overshoot = position_in_chrom - real_length
                    overflow[slot_idx] = (count + 1, max(max_overshoot, overshoot))

        if overflow:
            n_files_with_overflow += 1
            print(f"  {sdd_path}:")
            for slot_idx in sorted(overflow):
                count, max_overshoot = overflow[slot_idx]
                chrom_name = chrom_names[slot_idx]
                print(
                    f"    chr{chrom_name} (slot {slot_idx + 1}): {count} damage(s) beyond the real "
                    f"GRCh37 length ({GRCH37_CHROM_LENGTHS_BP[chrom_name]:,} bp; SDD declares "
                    f"{chrom_sizes_bp[slot_idx]:,} bp) -- worst overshoot {max_overshoot:,} bp"
                )

    print(
        f"\n{n_files_with_overflow} of {len(sdd_paths)} SDD file(s) contain damage positions "
        "beyond the real GRCh37 chromosome length."
    )


if __name__ == "__main__":
    # density = 1/50000
    # for ssb_density in [1/250, 1/500, 1/1000, 1/10000, 1/100000, 1/50000, 1/1000000]:
    #     genome_length = 3.1e9
    #     n_dsbs = int(density*genome_length)
    #     n_ssbs = int(ssb_density*genome_length)
    #     set_up_simulation_folder(n_dsbs, n_ssbs, 3)

    # for density in [1/500, 1/1000, 1/1e4, 1/1e5, 1/5e5, 1/1e6]:
    #     genome_length = 3.1e9 * 2
    #     n_dsbs = int(density*genome_length)
    #     # n_ssbs = int(ssb_density*genome_length)
    #     set_up_simulation_folder(n_dsbs, 0, 3)

    set_up_final_sdds_simulations()