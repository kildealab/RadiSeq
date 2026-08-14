RUN_RESULTS_FOLDER = "run_results"

SDD_TEMPLATE = "analysis_code/template_files/sdd_template.sdd"

PARAMETERS_TEMPLATE = "analysis_code/template_files/parameters_template.txt"

FIGURES_FOLDER = "analysis_code/figures"

SDD_FILES = "Final_SDDs"

FINAL_SDDS_RUN_RESULTS_FOLDER = "Final_SDDs_runs_results"

FINAL_SDDS_PARAMETERS_TEMPLATE = "Final_SDDs_runs_results/parameters_template.txt"

FINAL_SDDS_GRAPHS_FOLDER = "Final_SDDs_graphs"

FINAL_SDDS_ANALYSIS_FOLDER = "Final_SDDs_analysis"

FINAL_SDDS_ALIGNMENT_DATA_FOLDER = "Final_SDDs_alignment_data"

ALIGNING_TEST_CHROM_SIZES_PATH = "aligning_test/accessory_files/hg19.chrom.sizes.bed"

# The experimental setups (particle, energy, dose) simulated under Final_SDDs_runs_results.
# Every folder listed here holds many numbered simulation-run folders (one per SDD file,
# e.g. ".../161/") for that setup, plus (once generate_dsb_location_graphs has been run) a
# "dsb_graphs" folder. Update this if a new setup (particle/energy/dose combination) is added.
FINAL_SDDS_SETUPS = {
    "neutron_1MeV_0.5Gy": f"{FINAL_SDDS_RUN_RESULTS_FOLDER}/Neutron/1MeV_outer/0_5Gy",
    "neutron_1MeV_1.5Gy": f"{FINAL_SDDS_RUN_RESULTS_FOLDER}/Neutron/1MeV_outer/1_5Gy",
    "neutron_1MeV_3.0Gy": f"{FINAL_SDDS_RUN_RESULTS_FOLDER}/Neutron/1MeV_outer/3_0Gy",
    "photon_250keV_0.5Gy": f"{FINAL_SDDS_RUN_RESULTS_FOLDER}/Photon/250keV_outer/0_5Gy",
    "photon_250keV_1.5Gy": f"{FINAL_SDDS_RUN_RESULTS_FOLDER}/Photon/250keV_outer/1_5Gy",
    "photon_250keV_3.0Gy": f"{FINAL_SDDS_RUN_RESULTS_FOLDER}/Photon/250keV_outer/3_0Gy",
    "photon_6MeV_0.5Gy": f"{FINAL_SDDS_RUN_RESULTS_FOLDER}/Photon/6MeV_outer/0_5Gy",
    "photon_6MeV_1.5Gy": f"{FINAL_SDDS_RUN_RESULTS_FOLDER}/Photon/6MeV_outer/1_5Gy",
    "photon_6MeV_3.0Gy": f"{FINAL_SDDS_RUN_RESULTS_FOLDER}/Photon/6MeV_outer/3_0Gy",
}