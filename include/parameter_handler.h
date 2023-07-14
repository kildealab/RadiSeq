#ifndef PARAMETER_HANDLER_H
#define PARAMETER_HANDLER_H

#include <iostream>
#include <string>
#include <vector>


class NGSParameters{
    std::string default_parameter_file;                                          // Object holding defualt parameter file path. Hard-coded value
    // Following variable objects hold respective parameter value(s)
    unsigned int random_seed;                                                    // RNG seed value. If user provides, fixed value will be stored, else random
    bool is_merge_damages;                                                       // Assigns 1 if parameter value is "True" or "true". Else 0. If damages from different SDD to be merged
    int num_particles_to_merge;                                                  // Intiger value indicating the number of particle-simulations to merge
    std::vector<std::string> names_of_particles_to_merge;                        // Vector to hold the list of particles for which 'is_merge_damages' applicable
    std::vector<double> relative_dose_contributions;                             // Vector to hold the relative dose contribution from each particle
    std::vector<std::string> sddfile_path;                                       // Vector to hold the paths to SDD files from different particles to construct a genome
    bool is_adjust_dose;                                                         // If damages needs to be adjusted with actual dose delivered
    std::vector<std::string> actual_dosefile_path;                               // Vector to hold the paths to files containing actual dose deposited in each run by each particle
    std::string reference_genome_file;                                           // Variable holding reference genome fasta file path
    std::string output_directory;                                                // Path to the output directory where all the fastq files generated should be stored
    //int dsb_threshold;                                                           // This will hold the threshold distance between two opposite SSBs to form a DSB, in bp.
    std::string sequencer;                                                       // Name of the illumina sequencer to be used. Must be one with the build-in error profiles
    std::vector<std::string> list_sequencers;                                    // The list of names of all the built-in sequencer profiles we have
    int read_length;                                                             // Automatically set the read length based on the user provided sequencer name.
    std::vector<int> list_read_lengths;                                          // The list of read lengths of all the built-in sequencers in order of their names
    std::string r1_quality_profile;                                              // File name of the read 1 quality profile of the sequencer user-selected 
    std::vector<std::string> list_r1_quality_profiles;                           // The list of file names of read 1 quality profiles of all the built-in sequencers in order of their names
    std::string r2_quality_profile;                                              // File name of the read 2 quality profile of the sequencer user-selected 
    std::vector<std::string> list_r2_quality_profiles;                           // The list of file names of read 2 quality profiles of all the built-in sequencers in order of their names
    std::string sequencing_mode;                                                 // This value sets if user wants single-cell or bulk-cell sequencing
    int num_of_cells_in_sample;                                                  // Total number of cells we assume to have in a sample. Cells to sequence will be randomly sampled from this sample
    int num_of_cells_to_sequence;                                                // Total number of cells to be sequenced (damaged or not) at the end of a single run of this program
    int total_read_coverage;                                                     // Total read coverage: collectively from all cells sequenced
    bool is_paired_end_seq;                                                      // True if user wishes to have paired-end sequencing
    int mean_DNA_fragment_length;                                                // Mean DNA fragment length in bp (Mandatory input if paired-end)
    int std_dev_DNA_fragment_length;                                             // Standard deviation in mean DNA fragment length in bp (Mandatory input if paired-end)
    int max_read_errors;                                                         // Maximum number of errors the program can incorporate in a generated read
    double r1_insError_rate;                                                     // Insertion error rate in read 1 
    double r1_delError_rate;                                                     // Deletion error rate in read 1 
    double r2_insError_rate;                                                     // Insertion error rate in read 2 
    double r2_delError_rate;                                                     // Deletion error rate in read 2 
    std::string fastq_filename_prefix;                                           // String to hold the user-specified fastq output filename prefix
    bool is_summary_report;                                                      // True if user wishes to generate a summary report at the end of the run

public:
    NGSParameters();                                                             // Default constructor
    void process_parameterFile(const std::string*, NGSParameters&);              // Function to set all the user-defined sequencing parameters
    void set_parameters(std::string*, std::string*);                             // Mother function, when called will initiate all set functions

    void set_random_seed(std::string*);                                          // function to set the RNG seed if there is one provided
    unsigned int get_random_seed();                                              // functio to get the random seed

    void set_merge_damages_from_particles(std::string*, std::string*);           // function to set 'is_merge_damages' 
    bool get_merge_damages_from_particles();                                     // function to get 'is_merge_damages' 

    void set_num_of_particles_to_merge(std::string*);                            // function to set 'num_particles_to_merge'
    int get_num_of_particles_to_merge();                                         // function to get 'num_particles_to_merge'

    void set_names_of_particles_to_merge(std::string*);                          // function to set 'names_of_particles_to_merge' 
    const std::vector<std::string>* get_names_of_particles_to_merge();           // function to get 'names_of_particles_to_merge'

    void set_relative_dose_contributions(std::string*);                          // function to set 'relative_dose_contributions'
    std::vector<double>& get_relative_dose_contributions();                      // function to get 'relative_dose_contributions'

    void set_sddfile_path(std::string*);                                         // function to set 'sddfile_path'
    std::vector<std::string>& get_sddfile_path();                                // function to get 'sddfile_path'

    void set_adjust_damages_with_actual_dose(std::string*, std::string*);        // function to set 'adjust_damages_with_actual_dose'
    bool get_adjust_damages_with_actual_dose();                                  // function to get 'adjust_damages_with_actual_dose'

    void set_actual_dosefile_path(std::string*);                                 // function to set 'actual_dosefile_path'
    const std::vector<std::string>* get_actual_dosefile_path();                  // function to get 'actual_dosefile_path'

    void set_reference_genome(std::string*, std::string*);                       // function to set the reference genome path
    const std::string* get_reference_genome();                                   // function to get the reference genome path

    void set_output_directory(std::string*);                                     // function to set the path to output directory
    const std::string* get_output_directory();                                   // function to get the output directory path

    //void set_dsb_threshold(std::string*);                                        // function to set the DSB threshold value
    //int get_dsb_threshold();                                                     // function to get the DSB threshold values

    void set_sequencer(std::string*);                                            // function to set the name of the illumina sequencer
    const std::string* get_sequencer();                                          // function to get the sequencer name

    void set_read_length(int);                                                   // function to automatically set the read length 
    int get_read_length();                                                       // function to get the read length 

    void set_read_quality_profiles(int);
    const std::string* get_r1_quality_profile();
    const std::string* get_r2_quality_profile();

    void set_sequencing_mode(std::string*, std::string*);                        // function to set the sequencing mode: single-cell or bulk-cell
    const std::string* get_sequencing_mode();                                    // function to get the sequencing mode

    void set_num_of_cells_in_sample(std::string*, std::string*);                 // function to set the 'num_of_cells_in_sample'
    int get_num_of_cells_in_sample();                                            // function to get the 'num_of_cells_in_sample'

    void set_num_of_cells_to_sequence(std::string*, std::string*);               // function to set the 'num_of_cells_to_sequence'
    int get_num_of_cells_to_sequence();                                          // function to get the 'num_of_cells_to_sequence'

    void set_total_read_coverage(std::string*, std::string*);                    // function to set the read coverage value
    int get_total_read_coverage();                                               // function to get the total read coverage

    void set_paired_end_sequencing(std::string*, std::string*);                  // function to set 'is_paired_end_seq'
    bool get_paired_end_sequencing();                                            // function to get 'is_paired_end_seq'

    void set_mean_DNA_fragment_length(std::string*);                             // function to set the mean DNA fragment length. Different from read length
    int get_mean_DNA_fragment_length();                                          // function to get the mean DNA fragment length

    void set_std_dev_DNA_fragment_length(std::string*);                          // function to set the standard deviation in mean DNA fragment length
    int get_std_dev_DNA_fragment_length();                                       // function to get the standard deviation in mean DNA fragment length

    void set_max_errors_in_read(std::string*);                                   // function to set the maximum number of errors in a read
    int get_max_errors_in_read();                                                // function to get the maximum number of errors in a read

    void set_insertion_error_rate_read1(std::string*);                           // function to set the insertion error rate in read 1 
    double get_insertion_error_rate_read1();                                     // function to get the insertion error rate in read 1 
    
    void set_deletion_error_rate_read1(std::string*);                            // function to set the deletion error rate in read 1 
    double get_deletion_error_rate_read1();                                      // function to get the deletion error rate in read 1 
    
    void set_insertion_error_rate_read2(std::string*);                           // function to set the insertion error rate in read 2 
    double get_insertion_error_rate_read2();                                     // function to get the insertion error rate in read 2
    
    void set_deletion_error_rate_read2(std::string*);                            // function to set the deletion error rate in read 2 
    double get_deletion_error_rate_read2();                                      // function to get the deletion error rate in read 2 

    void set_output_fastq_filename_prefix(std::string*);                         // function to set the fastq output filename prefix
    const std::string* get_output_fastq_filename_prefix();                       // function to get the fastq output filename prefix

    void set_summary_report(std::string*, std::string*);                                       // function to set "is_summary_report"
    bool get_summary_report();                                                   // function to get "is_summary_report"

    void help_parameter(std::string*);                                           // Function to print help message for every parameter
    void success_parameter();                                                    // Function to check the appropriateness of all parameters

};





#endif