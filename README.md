# radiSeq_Simulator

This is a tool to simulate Next Generation Sequencing (NGS) of radiation-exposed cells using a Standard DNA Damage (SDD) data file from Monte Carlo simulations of cell irradiations. radiSeq Simulator can perform both bulk- and single-cell whole genome DNA sequencing.


## Table of Contents

* [Authors](#authors)
* [Description](#description)
* [Features](#features)
* [Installation](#installation)
* [Input Parameters](#input-parameters)
* [Acknowledgements](#acknowledgements)

## Authors

Felix Mathew and John Kildea

Contact email: felix.mathew@mail.mcgill.ca

Website: [www.kildealab.com](https://kildealab.com/software/radiseq_simulator/)

## Description
Use radiSeq Simulator to computationally simulate whole genome DNA sequencing of radiation-exposed cells in a sample. The complete working logic is shown in the flowchart below. <br>

![Logo](https://github.com/felixmat/radiSeq_Simulator/blob/main/radiSeq%20Simulator.svg)

## Features

* All input values can be specified in a single parameter file
* Compatible with any nuclear DNA model that can output an SDD
* Multi-threading enabled 
* Option to switch between single-cell and bulk-cell sequencing
* Option to perform paired-end sequencing
* 6 built-in Illumina sequencers to choose from
* Users can use any reference genome of their choosing
* Option to generate a detailed run summary output file

## Installation
### Prerequisites

* Compiler: Should support C++17 or above
* OS: Unix-like systems (e.g.: Linux, macOS)

**Note**: This application was developed on macOS Sierra 10.12.6.

### Getting started

1. Download the latest version of radiSeq_Simulator
2. Unzip the downloaded radiSeq_Simulator file
3. Download the generic human reference genome [(click here to download)](https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_rm.primary_assembly.fa.gz)<br>
   Unzip the downloaded file and save it in the 'radiSeqData' directory under the name `'Human_reference_genome.fa'`.
5. Compile radiSeq_Simulator:
   * `cd path/to/radiSeq_Simulator`
   * `make`<br>
   Ignore any warnings that you might see.
6. Set up the environment variable 'RADISEQ_DATA_DIR'
   * `export RADISEQ_DATA_DIR=path/to/radiSeq_Simulator/radiSeqData`<br>
   You will need to do this step every time you open a fresh Terminal window. Alternatively, you may choose to add this into one of your statup files (eg: .bashrc) if you are comfortable doing so.<br>
**Note**: Replace 'path/to/radiSeq_Simulator' in step 5 and 6 with the actual full path to the directory radiSeq_Simulator in your system

### Running a test
Run the test sequencing to check if the radiSeq_Simulator is working fine. You can run the test using the command:<br> 
* `cd ./example_test`
* `../radiSeq ./TestParameters.txt`

A successful test run will generate FASTQ output files and a run summary file in a folder called 'output' inside the example_test directory.

**Note**: User can create their own [input parameter files](#input-parameters) to run specific simulations 

## Input parameters

The user should specify all the input parameters for the simulation in a parameter text(.txt) file. Users can make any text file a parameter file with a filename
of their choosing, as long as the contents of the file are formatted in a specified parameter name and value pair. Input values can be given in the format:
`Parameter Name = Parameter Value   #comment`. The list of all acceptable parameter names is given in the table below. All the parameters are optional except the path to the SDD file. 

| Parameter name | Description | Parameter value |
|----------|------------|------------|
| sddFilePath | Complete path to the SDD file(s) | Comma separated list of paths (string) |
| merge_damages_from_multiple_particles | Flag to indicate if the user wishes to define a single genome by combining damages from multiple SDD files | 'True' or 'False' |
| number_of_particles_to_merge | Number of SDD files from individual primary simulations to combine | Number (integer) |
| primary_particles_simulated | Names of primary particles that introduced the damages that are going to be combined into a single genome | Comma-separated list of names  (strings) |
| relative_dose_contributions | Relative dose contributions of each primary particle in order of the SDD files | Comma-separated list of numbers between 0 and 1 (float(s)) |
| adjust_damages_with_actual_dose | Flag to indicate if the user wishes to scale the number of damages with the actual dose delivered and it is different from the expected dose | 'True' or 'False' |
| actual_dose_delivered_data | Complete path to the file containing the actual dose delivered in each run. One file is expected for each SDD file specified | Comma-separated list of paths (string) |
| reference_genome_FASTAfile | Complete path to the reference genome file | path to file (sting) |
| acceptable_difference_in_seq_length_percent | Acceptable difference in the lengths of the reference genome provided and the genome length of the Monte Carlo model | Percentage (double) |
| number_of_cells_in_sample | Total number of cells the user assumes to have in your sample. This is different from the number of cells to sequence | Number (integer) |
| number_of_cells_to_sequence | Number of cells (damaged and undamaged) to be sequenced. These many cells will be randomly selected from the number of cells in the sample | Number (integer) |
| illumina_sequencer | Name of the Illumina sequencer to be used for sequencing from the in-built list | 'HiSeq1000', 'HiSeq2000', 'HiSeq2500_v125', 'HiSeq2500_v150', 'HiSeqX' and 'NovaSeq6000' |
| single_or_bulk_sequencing | Flag to specify if single-cell or bulk-cell sequencing is to be performed | 'single' or 'bulk' |
| do_paired_end_sequencing | Flag to indicate if paired-end sequencing to be performed | 'True' or 'False' |
| mean_DNA_fragment_length | Mean DNA fragment length (in bp) if paired-end sequencing. This is different from the read length | Number (integer) |
| std_dev_DNA_fragment_length | Standard deviation in the mean DNA fragment length (in bp) | Number (integer) |
| total_read_coverage | Total read coverage the user wants to get from this sequencing. If single-cell sequencing the read coverage will get distributed over the total number of cells sequenced | Number (integer) |
| maximum_errors_in_reads | Maximum number of indel errors to have in a read. A default value of -1 indicates that there is no limit or threshold to the number of errors in a read | Number (integer) |
| read1_insertion_error_rate | Insertion error rate for read 1. Expects a double value between 0 and 1 | Number (double) |
| read1_deletion_error_rate | Deletion error rate for read 1. Expects a double value between 0 and 1 | Number (double) |
| read2_insertion_error_rate | Insertion error rate for read 2. Expects a double value between 0 and 1 | Number (double) |
| read2_deletion_error_rate | Deletion error rate for read 2. Expects a double value between 0 and 1 | Number (double) |
| output_directory_path | Complete path to the directory where the output fastq files and the run summary file should be stored | Path to directory (string) |
| output_FASTQ_filename_prefix | Prefix for the sequenced output FASTQ file (omit file extension) | String |
| make_summary_report | Flag to indicate if the user wishes to generate a summary report file at the end of run | 'True' or 'False' |
| random_seed | Seed number for the random number generator to be initialized with a fixed seed. A default value of 0 indicates that the system will automatically generate random seeds completely random | Number (integer) |
| number_of_threads | Number of threads to be used for a multithreaded run. Default value is 1 | Number (integer) |


## Acknowledgements

The sequencing segment of the radiSeq Simulator was inspired by the ART Illumina sequencing tool published open-source by the US National Institute of Health.
Some of the in-built Illumina sequencer error profiles and the indel error rates used in the radiSeq Simulator were adopted from the ART sequencing tool. 
We sincerely thank the authors Weichun Huang et al. for making the ART toolkit available open-source. 

To learn more about the ART toolkit:

* Publication: Weichun Huang and others, ART: a next-generation sequencing read simulator, Bioinformatics, Volume 28, Issue 4, February 2012, Pages 593–594, https://doi.org/10.1093/bioinformatics/btr708
* ART toolkit: https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm
* GitHub: https://github.com/scchess/Art/tree/master






