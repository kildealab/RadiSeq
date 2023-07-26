# SDD-NGS_Simulator




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

Website: www.kildealab.com

## Description


## Features

* All input values can be specified in a single parameter file
* Compatible with any nuclear DNA model that can output an SDD
* Option to switch between single-cell and bulk-cell sequencing
* Option to perform paired-end sequencing
* 5 built-in Illumina sequencers to choose from
* User can use any reference genome of their choosing
* Option to generate a detailed run summary output file

## Installation
### Prerequisites

* Compiler: Should support C++11 or above
* OS: Unix-like systems (e.g.: Linux, macOS)

**Note**: This application was developed on macOs Sierra 10.12.6.

### Getting started

1. Download the latest version of SDD-NGS_Simulator
2. Unzip the downloaded SDD-NGS_Simulator file
   * Also unzip `Human_reference_genome.fa.zip` in NGSSimData directory
3. Compile SDD-NGS_Simulator:
   * `cd path/to/SDD-NGS_Simulator`
   * `make`

Ignore any warnings that you might see.

### Running a test
Run the test sequencing to check if the SDD-NGS_Simulator is working fine. You can run the test using the command:
  * `./SDD_NGSSimulator ./inputfiles/TestParameters.txt`

A successful test run will generate FASTQ output files and a run summary file in a folder called 'output'.

**Note**: User can create their own [input parameter files](#input-parameters) to run specific simulations 

## Input parameters

The user should specify all the input parameters for the simulation in a parameter text(.txt) file. Users can make any text file a parameter file with a filename
of their choosing, as long as the contents of the file are formatted in a specified parameter name and value pair. Input values can be given in the format:
`_Parameter Name = Parameter Value   #comment_`. The list of all acceptable parameter names is given in the table below. All the parameters are optional except the path to the SDD file. 

| Parameter name | Description | Parameter value |
|----------|------------|------------|
| sddFilePath | Complete path to the SDD file(s) | Comma seperated list of paths (string) |
| merge_damages_from_multiple_particles | Flag to indicate if the user wishes to define a single genome by combiining damages from multiple SDD files | 'True' or 'False' |
| number_of_particles_to_merge | Number of SDD files from individual primary simulations to combine | Number (integer) |
| primary_particles_simulated | Names of primary particle that introduced the damages that are going to be combined into a single genome | Comma seperated list of names  (strings) |
| relative_dose_contributions | Relative dose contributions of each primary particle in order of the SDD files | Comma seperated list of numbers between 0 and 1 (float(s)) |
| adjust_damages_with_actual_dose | Flag to indicate if the user wishes to scale the number of damages with the actual dose delivered and it is different from expected dose | 'True' or 'False' |
| actual_dose_delivered_data | Complete path to the file containing actual dose delivered in each run. One file is expected for each SDD file specified | Comma seperated list of paths (string) |
| reference_genome_FASTAfile | Complete path to the reference genome file | path to file (sting) |
| acceptable_difference_in_seq_length_percent | Acceptable difference in the lengths of the reference genome provided and the genome length of the Monte Carlo model | Percentage (double) |
| number_of_cells_in_sample | Total number of cells user assumes to have in your sample. This is different from the number of cells to sequence | Number (integer) |
| number_of_cells_to_sequence | Number of cells (damaged and undamaged) to be sequenced. These many cells will be randomly selected from the number of cells in sample | Number (integer) |
| illumina_sequencer | Name of the Illumina sequencer to be used for sequencing from the in-built list | 'HiSeq1000', 'HiSeq2000', 'HiSeq2500_v125', 'HiSeq2500_v150' and 'HiSeqX' |
| single_or_bulk_sequencing | Flag to specify if single-cell or bulk-cell sequencing to be performed | 'single' or 'bulk' |
| do_paired_end_sequencing | Flag to indicate if paired-end sequencing to be performed | 'True' or 'False' |
| mean_DNA_fragment_length | Mean DNA fragment length (in bp) if paired-end sequencing. This is different from the read length | Number (integer) |
| std_dev_DNA_fragment_length | Standard deviation in the mean DNA fragment length (in bp) | Number (integer) |
| total_read_coverage | Total read coverage user wants to get from this sequencing. If single-cell sequencing the read coverage will get distributed over the total number of cells seqeunced | Number (integer) |
| maximum_errors_in_reads | Maximum number of indel errors to have in a read. Default value of -1 indicates that there is no limit or threshold to the number of errors in a read | Number (integer) |
| read1_insertion_error_rate | Insertion error rate for read 1. Expects a double value between 0 and 1 | Number (double) |
| read1_deletion_error_rate | Deletion error rate for read 1. Expects a double value between 0 and 1 | Number (double) |
| read2_insertion_error_rate | Insertion error rate for read 2. Expects a double value between 0 and 1 | Number (double) |
| read2_deletion_error_rate | Deletion error rate for read 2. Expects a double value between 0 and 1 | Number (double) |
| output_directory_path | Complete path to the directory where the output fastq files and the run summary file should be stored | Path to directory (string) |
| output_FASTQ_filename_prefix | Prefix for the sequenced output FASTQ file (omit file extension) | String |
| make_summary_report | Flag to indicate if the user wishes to generate a summary report file at the end of run | 'True' or 'False' |
| random_seed | Seed number for the random number generator to be initialized with a fixed seed. Default value of 0 indicates that the system will automatically generate random seeds completely random | Number (integer) |


## Acknowledgements

The sequencing segment of the SDD_NGS Simulator was inspired by the ART Illumina sequencing tool published open-source by the US National Institute of Health.
Some of the in-built Illumina sequencer error profiles and the indel error rates used in SDD_NGS Simulator were adopted from the ART sequencing tool. 
We sincerely thank the authors Weichun Huang et al. for making the ART toolkit available open-source. 

To learn more about the ART toolkit:

* Publication: Weichun Huang and others, ART: a next-generation sequencing read simulator, Bioinformatics, Volume 28, Issue 4, February 2012, Pages 593–594, https://doi.org/10.1093/bioinformatics/btr708
* ART toolkit: https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm
* GitHub: https://github.com/scchess/Art/tree/master






