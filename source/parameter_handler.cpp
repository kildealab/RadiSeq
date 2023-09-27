#include "parameter_handler.h"
#include "fileio.h"
#include "support_functions.h"
#include "random_generator.h"

#include <iostream>
#include <string>
#include <algorithm>

// Default constructor
NGSParameters::NGSParameters(){
    default_parameter_file = "./NGSSimData/NGSDefaultParameters.txt";
    reference_genome_file = "./NGSSimData/Human_reference_genome.fa";
    list_sequencers = {"HiSeq1000","HiSeq2000","HiSeq2500_v125","HiSeq2500_v150","HiSeqX","NovaSeq6000","test"};
    list_read_lengths = {100, 100, 125, 150, 150, 150, 5};
    list_r1_quality_profiles = {"HiSeq1000_R1.txt","HiSeq2000_R1.txt","HiSeq2500_v125_R1.txt","HiSeq2500_v150_R1.txt","HiSeqX_R1.txt","NovaSeq6000_R1.txt","test_R1.txt"};
    list_r2_quality_profiles = {"HiSeq1000_R2.txt","HiSeq2000_R2.txt","HiSeq2500_v125_R2.txt","HiSeq2500_v150_R2.txt","HiSeqX_R2.txt","NovaSeq6000_R2.txt","test_R2.txt"};
}



//--------------------------------------------------------------------------------------------
// This function takes user-specified input file as argument and process the file
//--------------------------------------------------------------------------------------------
void NGSParameters::process_parameterFile(const std::string* parameterfile, NGSParameters& parameter){
    // Parameter file is a MANDATORY requirement since we need to know SDD file path. If not found, exit with an error.
    if (!checkFileExists(parameterfile)){
        std::cerr<<"\n ERROR: Unable to read the parameter file : "<< *parameterfile<<'\n';
        exit(EXIT_FAILURE);
    }else{
        std::cout<<"\n Successfully read the parameter file: "<< *parameterfile<<'\n';
    }
    
    // Continue if parameter file is present and set all parameter values
    readParameterFile(&default_parameter_file, parameter);         // Read default parameter file first to set default parameter values
    readParameterFile(parameterfile, parameter);                   // Read user-specified parameter file to overwrite default parameter values

    success_parameter();                                           // Check if provided parameters meet all conditions and if it is appropriate

    if(get_random_seed()){                                         // If user provided a non-zero seed value for the RNG
        rng::initFixedSeed(get_random_seed());                     // Initiate the RNG with the fixed seed
    }else{
        rng::initRandomSeed();                                     // Otherwise, initiate the RNG with a random seed
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// Global set_parameter function that calls the respective set_parameter functions depending
// on the parameterName passed
//--------------------------------------------------------------------------------------------
void NGSParameters::set_parameters(std::string* paramName, std::string* paramValue){
    if (*paramName == "random_seed"){
        set_random_seed(paramValue);
    }
    else if (*paramName == "merge_damages_from_multiple_particles"){
        set_merge_damages_from_particles(paramName, paramValue);
    } 
    else if (*paramName == "number_of_particles_to_merge"){
        set_num_of_particles_to_merge(paramValue);
    } 
    else if (*paramName == "primary_particles_simulated"){
        set_names_of_particles_to_merge(paramValue);
    } 
    else if (*paramName == "relative_dose_contributions"){
        set_relative_dose_contributions(paramValue);
    } 
    else if (*paramName == "sddFilePath"){
        set_sddfile_path(paramValue);
    }
    else if (*paramName == "adjust_damages_with_actual_dose"){
        set_adjust_damages_with_actual_dose(paramName, paramValue);
    }
    else if (*paramName == "actual_dose_delivered_data"){
        set_actual_dosefile_path(paramValue);
    }
    else if (*paramName == "reference_genome_FASTAfile"){
        set_reference_genome(paramName, paramValue);
    }
    else if(*paramName == "acceptable_difference_in_seq_length_percent"){
        set_max_acceptable_seq_length_difference(paramValue);
    }
    else if (*paramName == "output_directory_path"){
        set_output_directory(paramValue);
    }
    //else if (*paramName == "DSB_threshold_in_bp"){
    //    set_dsb_threshold(paramValue);
    //}
    else if (*paramName == "illumina_sequencer"){
        set_sequencer(paramValue);
    }
    else if (*paramName == "single_or_bulk_sequencing"){
        set_sequencing_mode(paramName, paramValue);
    }
    else if (*paramName == "number_of_cells_in_sample"){
        set_num_of_cells_in_sample(paramName, paramValue);
    }
    else if (*paramName == "number_of_cells_to_sequence"){
        set_num_of_cells_to_sequence(paramName, paramValue);
    }
    else if (*paramName == "total_read_coverage"){
        set_total_read_coverage(paramName, paramValue);
    }
    else if (*paramName == "do_paired_end_sequencing"){
        set_paired_end_sequencing(paramName, paramValue);
    }
    else if (*paramName == "mean_DNA_fragment_length"){
        set_mean_DNA_fragment_length(paramValue);
    }
    else if (*paramName == "std_dev_DNA_fragment_length"){
        set_std_dev_DNA_fragment_length(paramValue);
    }
    else if (*paramName == "maximum_errors_in_reads"){
        set_max_errors_in_read(paramValue);
    }
    else if (*paramName == "read1_insertion_error_rate"){
        set_insertion_error_rate_read1(paramValue);
    }
    else if (*paramName == "read1_deletion_error_rate"){
        set_deletion_error_rate_read1(paramValue);
    }
    else if (*paramName == "read2_insertion_error_rate"){
        set_insertion_error_rate_read2(paramValue);
    }
    else if (*paramName == "read2_deletion_error_rate"){
        set_deletion_error_rate_read2(paramValue);
    }
    else if (*paramName == "output_FASTQ_filename_prefix"){
        set_output_fastq_filename_prefix(paramValue);
    }
    else if (*paramName == "make_summary_report"){
        set_summary_report(paramName, paramValue);
    }
    else{
        std::cerr<<"\n WARNING: Unrecognized parameter specified : \""<<*paramName<<"\"\n"
        <<" ----- This parameter will be ignored -----\n";
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// set_parameter functions : They assign parameter values. User-defined values override default
//--------------------------------------------------------------------------------------------
void NGSParameters::set_random_seed(std::string* paramValue){
    random_seed = static_cast<unsigned int>(std::stoi(*paramValue));
}
void NGSParameters::set_merge_damages_from_particles(std::string* paramName,std::string* paramValue){
    if(lowercaseString(paramValue) == "true"||lowercaseString(paramValue) == "false"){
        is_merge_damages = (lowercaseString(paramValue) == "true");                           // Gets 1 if "True" or "true"; else 0
    }else{                                                                                    // If user provided unrecognized value format, print help
        help_parameter(paramName);
        std::cerr<<" ----- Setting \""<<*paramName<<"\" to its default value: \""<<std::boolalpha<<get_merge_damages_from_particles()<<"\" -----\n";
    }
}
void NGSParameters::set_num_of_particles_to_merge(std::string* paramValue){
    if(std::stoi(*paramValue) > 0){                                                           // set this parameter only if user specified a non-zero value. Else default
        num_particles_to_merge = atoi(paramValue->c_str());                                   // Converting the parameterValue to an int
    }
}
void NGSParameters::set_names_of_particles_to_merge(std::string* paramValue){
    stringToVec(',', paramValue, names_of_particles_to_merge);                                // Converting comma seperated list of names to a vector
}
void NGSParameters::set_relative_dose_contributions(std::string* paramValue){
    stringToVec(',', paramValue, relative_dose_contributions);                                // Registering comma seperated double values as a vector
}
void NGSParameters::set_sddfile_path(std::string* paramValue){
    stringToVec(',', paramValue, sddfile_path);                                               // Registering comma seperated SDD paths to a vector
}
void NGSParameters::set_adjust_damages_with_actual_dose(std::string* paramName, std::string* paramValue){
   if(lowercaseString(paramValue) == "true"||lowercaseString(paramValue) == "false"){
        is_adjust_dose = (lowercaseString(paramValue) == "true");                              // Gets 1 if "True" or "true"; else 0
    }else{                                                                                     // If user provided unrecognized value format, print help
        help_parameter(paramName);
        std::cerr<<" ----- Setting \""<<*paramName<<"\" to its default value: \""<<std::boolalpha<<get_adjust_damages_with_actual_dose()<<"\" -----\n";
    } 
}
void NGSParameters::set_actual_dosefile_path(std::string* paramValue){
    stringToVec(',', paramValue, actual_dosefile_path);                                       // Registering comma seperated actual_dosefile paths to a vector
}
void NGSParameters::set_reference_genome(std::string* paramName, std::string* paramValue){
    int param_len = paramValue->size();
    if(param_len>3 && paramValue->substr(param_len-3) == ".fa"){                              // If only filepath value has minimum three chars, and that three chars are '.fa' at the end
        reference_genome_file = *paramValue;
    }else if(*paramValue == "default"){                                                       // If reading default parameter file, ignore
        return;
    }else{
        help_parameter(paramName);
        std::cerr<<" ----- Setting \""<<*paramName<<"\" to its default value: Using Ashkenazi human reference genome -----\n";
    }
}
void NGSParameters::set_max_acceptable_seq_length_difference(std::string* paramValue){
    max_diff_model_Vs_reference = std::stod(*paramValue);
}
void NGSParameters::set_output_directory(std::string* paramValue){
    output_directory = *paramValue;
}
//void NGSParameters::set_dsb_threshold(std::string* paramValue){
//    dsb_threshold = std::stoi(*paramValue);
//}
void NGSParameters::set_sequencer(std::string* paramValue){
    sequencer = *paramValue;
}
void NGSParameters::set_read_length(int sequencer_index){                                     // This parameter is automatically set based on the user provided sequencer
    read_length = list_read_lengths[sequencer_index];
}
void NGSParameters::set_read_quality_profiles(int sequencer_index){
    r1_quality_profile = "./NGSSimData/"+list_r1_quality_profiles[sequencer_index];
    r2_quality_profile = "./NGSSimData/"+list_r2_quality_profiles[sequencer_index];
}
void NGSParameters::set_sequencing_mode(std::string* paramName, std::string* paramValue){
    if(lowercaseString(paramValue) == "single"||lowercaseString(paramValue) == "bulk"){
        sequencing_mode = lowercaseString(paramValue);                                        // Gets 'single' or 'bulk' depending on the user provided sequencing mode
    }else{                                                                                    // If user provided unrecognized value format, print help
        help_parameter(paramName);
        std::cerr<<" ----- Setting \""<<*paramName<<"\" to its default value: \"single\" -----\n";
    }
}
void NGSParameters::set_num_of_cells_in_sample(std::string* paramName, std::string* paramValue){
    if(std::stoi(*paramValue) > 0){
        num_of_cells_in_sample = std::stoi(*paramValue);
    }else{
        help_parameter(paramName);
    }
}
void NGSParameters::set_num_of_cells_to_sequence(std::string* paramName, std::string* paramValue){
    if(std::stoi(*paramValue) > 0){
        num_of_cells_to_sequence = std::stoi(*paramValue);
    }else{
        help_parameter(paramName);
    }
}
void NGSParameters::set_total_read_coverage(std::string* paramName, std::string* paramValue){
    if(std::stoi(*paramValue) > 0){
        total_read_coverage = std::stoi(*paramValue);
    }else{
        help_parameter(paramName);
    }
}
void NGSParameters::set_paired_end_sequencing(std::string* paramName, std::string* paramValue){
    if(lowercaseString(paramValue) == "true"||lowercaseString(paramValue) == "false"){
        is_paired_end_seq = (lowercaseString(paramValue) == "true");                           // Gets 1 if "True" or "true"; else 0
    }else{                                                                                     // If user provided unrecognized value format, print help
        help_parameter(paramName);
        std::cerr<<" ----- Setting \""<<*paramName<<"\" to its default value: \"True\" -----\n";
    }
}
void NGSParameters::set_mean_DNA_fragment_length(std::string* paramValue){
    mean_DNA_fragment_length = std::stoi(*paramValue);
}
void NGSParameters::set_std_dev_DNA_fragment_length(std::string* paramValue){
    std_dev_DNA_fragment_length = std::stoi(*paramValue);
}
void NGSParameters::set_max_errors_in_read(std::string* paramValue){
    max_read_errors = std::stoi(*paramValue);
}
void NGSParameters::set_insertion_error_rate_read1(std::string* paramValue){
    r1_insError_rate = std::stod(*paramValue);
}
void NGSParameters::set_deletion_error_rate_read1(std::string* paramValue){
    r1_delError_rate = std::stod(*paramValue);
}
void NGSParameters::set_insertion_error_rate_read2(std::string* paramValue){
    r2_insError_rate = std::stod(*paramValue);
}
void NGSParameters::set_deletion_error_rate_read2(std::string* paramValue){
    r2_delError_rate = std::stod(*paramValue);
}
void NGSParameters::set_output_fastq_filename_prefix(std::string* paramValue){
    fastq_filename_prefix = *paramValue;
}
void NGSParameters::set_summary_report(std::string* paramName, std::string* paramValue){
    if(lowercaseString(paramValue) == "true"||lowercaseString(paramValue) == "false"){
        is_summary_report = (lowercaseString(paramValue) == "true");                           // Gets 1 if "True" or "true"; else 0
    }else{                                                                                     // If user provided unrecognized value format, print help
        help_parameter(paramName);
        std::cerr<<" ----- Setting \""<<*paramName<<"\" to its default value: \"True\" -----\n";
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// get_parameter functions : They return the parameter values that are set
//--------------------------------------------------------------------------------------------
unsigned int NGSParameters::get_random_seed(){
    return(random_seed);
}
bool NGSParameters::get_merge_damages_from_particles(){
    return(is_merge_damages);
}
int NGSParameters::get_num_of_particles_to_merge(){
    return(num_particles_to_merge);
} 
const std::vector<std::string>* NGSParameters::get_names_of_particles_to_merge(){
    return(&names_of_particles_to_merge);
}
std::vector<double>&  NGSParameters::get_relative_dose_contributions(){
    return(relative_dose_contributions);
}
std::vector<std::string>& NGSParameters::get_sddfile_path(){
    return(sddfile_path);
}
bool NGSParameters::get_adjust_damages_with_actual_dose(){
    return(is_adjust_dose);
}
const std::vector<std::string>* NGSParameters::get_actual_dosefile_path(){
    return(&actual_dosefile_path);
}
const std::string* NGSParameters::get_reference_genome(){
    return(&reference_genome_file);
}
double NGSParameters::get_max_acceptable_seq_length_difference(){
    return(max_diff_model_Vs_reference);
}
const std::string* NGSParameters::get_output_directory(){
    return(&output_directory);
}
//int NGSParameters::get_dsb_threshold(){
//    return(dsb_threshold);
//}
const std::string* NGSParameters::get_sequencer(){
    return(&sequencer);
}
int NGSParameters::get_read_length(){
    return(read_length);
}
const std::string* NGSParameters::get_r1_quality_profile(){
    return(&r1_quality_profile);
}
const std::string* NGSParameters::get_r2_quality_profile(){
    return(&r2_quality_profile);
}
const std::string* NGSParameters::get_sequencing_mode(){
    return(&sequencing_mode);
}
int NGSParameters::get_num_of_cells_in_sample(){
    return(num_of_cells_in_sample);
}
int NGSParameters::get_num_of_cells_to_sequence(){
    return(num_of_cells_to_sequence);
}
int NGSParameters::get_total_read_coverage(){
    return(total_read_coverage);
}
bool NGSParameters::get_paired_end_sequencing(){
    return(is_paired_end_seq);
}
int NGSParameters::get_mean_DNA_fragment_length(){
    return(mean_DNA_fragment_length);
}
int NGSParameters::get_std_dev_DNA_fragment_length(){
    return(std_dev_DNA_fragment_length);
}
int NGSParameters::get_max_errors_in_read(){
    return(max_read_errors);
}
double NGSParameters::get_insertion_error_rate_read1(){
    return(r1_insError_rate);
}
double NGSParameters::get_deletion_error_rate_read1(){
    return(r1_delError_rate);
}
double NGSParameters::get_insertion_error_rate_read2(){
    return(r2_insError_rate);
}
double NGSParameters::get_deletion_error_rate_read2(){
    return(r2_delError_rate);
}
const std::string* NGSParameters::get_output_fastq_filename_prefix(){
    return(&fastq_filename_prefix);
}
bool NGSParameters::get_summary_report(){
    return(is_summary_report);
}
//--------------------------------------------------------------------------------------------




//--------------------------------------------------------------------------------------------
// This function will print help message for every parameter
//--------------------------------------------------------------------------------------------
void NGSParameters::help_parameter(std::string* paramName){
    std::cerr<<"\n WARNING: Unrecognizable/inappropriate value specified for parameter : \""<<*paramName<<"\"\n";
    if (*paramName == "merge_damages_from_multiple_particles"){
        std::cerr<<" This parameter should be set \"True\" or \"False\" "
        <<"to specify whether or not you wish to combine damages from multiple primary particles into a single genome\n";
    }
    else if (*paramName == "number_of_particles_to_merge"){
        std::cerr<<" This parameter should be set to a non-zero positive integer number "
        <<"if you wish to merge damages from multiple primary particles (from multiple SDD files) into a single genome\n";
    }
    else if (*paramName == "primary_particles_simulated"){
        std::cerr<<" Comma seperated list of names of primary particles (strings) that introduced  "
        <<"the damages that are going to be combined into a single genome.\n Specify one name, if there is only one primary (optional)\n";
    }
    else if (*paramName == "relative_dose_contributions"){
        std::cerr<<" Comma seperated list of relative dose contributions (floats) for each primary particle in order. \n"
        <<" Same order should be kept when specifying respective SDD files. Sum of these values should be 1. \n";
    }
    else if (*paramName == "sddFilePath"){
        std::cerr<<" You MUST specify the path to the SDD file. If you wish to combine damages from multiple SDD file into a single genome, \n"
        <<" then provide a comma seperated list of paths, in the same order in which their relative_dose_contributions are specified \n";
    }
    else if (*paramName == "adjust_damages_with_actual_dose"){
        std::cerr<<" This parameter should be set \"True\" or \"False\" "
        <<"to specify whether or not you wish to adjust the number of DNA damages with the actual dose delivered \n"
        <<" If True, you must provide the actual dose delivered data in each exposure in a single file in order\n";
    }
    else if (*paramName == "actual_dose_delivered_data"){
        std::cerr<<" Specify the path to the file containing actual dose delivered (in Gy) in each run/exposure. \n"
        <<" If there are more than one SDD file corresponding to different primary particle simulations, \n"
        <<" then, provice a comma seperated list of paths in the same order as the respective SDD file \n";
    }
    else if (*paramName == "reference_genome_FASTAfile"){
        std::cerr<<" Specify the path to the FASTA file with the reference genome. \n"
        <<" Specified reference genome must be in a valid FASTA file ending with '.fa' extension. \n"
        <<" If not specified, the default Ashkenazi human reference genome will be used\n";
    }
    else if (*paramName == "illumina_sequencer"){
        std::cerr<<" Specify the name of the Illumina sequencer to be used for sequencing. \n"
        <<" The name of the sequencer has to match one of the buil-in sequencer profiles. \n"
        <<" Built-in sequencers: HiSeq1000, HiSeq2000, HiSeq2500_v125, HiSeq2500_v150, HiSeqX and NovaSeq6000\n";
    }
    else if (*paramName == "single_or_bulk_sequencing"){
        std::cerr<<" Specify whether you want to perform single-cell or bulk-cell whole-genome sequencing. \n"
        <<" This parameter can take only values 'single' and 'bulk' for two different sequencing types. \n";
    }
    else if (*paramName == "number_of_cells_in_sample"){
        std::cerr<<" Specify the total number of cells you assume to have in your sample. \n"
        <<" Cells to sequence will be randomly sampled from a pool of these many cells. \n";
    }
    else if (*paramName == "number_of_cells_to_sequence"){
        std::cerr<<" Specify the total number of cells (damaged + undamaged) you want to sequence. \n"
        <<" These many cells will be randomly sampled from the total number of cells in the sample. \n"
        <<" This value is expected to be less than or equal to the total number of cells in the sample \n";
    }
    else if (*paramName == "total_read_coverage"){
        std::cerr<<" Specify the total read coverage you want to get from this sequencing (Integer value expected). \n"
        <<" If single-cell sequencing, the coverage gets distributed over the number of cells to be sequenced. \n"
        <<" i.e, if 100 cells sequenced with a coverage of 100x each cell will have 1x coverage in single-cell sequencing \n";
    }
    else if (*paramName == "do_paired_end_sequencing"){
        std::cerr<<" This parameter should be set \"True\" or \"False\" "
        <<"to specify whether or not you wish to perform paired-end sequencing \n";
    }
    else if (*paramName == "mean_DNA_fragment_length"){
        std::cerr<<" Specify the mean DNA fragment length in bp. This is different from the read length. \n"
        <<" A read is read from a DNA fragment, therefore, the mean DNA fragment length must be longer than the read length\n";
    }
    else if (*paramName == "std_dev_DNA_fragment_length"){
        std::cerr<<" Specify the standard deviation of the mean DNA fragment length in bp. \n";
    }
    else if (*paramName == "make_summary_report"){
        std::cerr<<" This parameter should be set \"True\" or \"False\" "
        <<"to specify whether or not you wish to generate a summary report at the end of the run \n";
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function will check the appropriateness of the parameters provided. 
// If there are critical mistakes, then exit with error.
//--------------------------------------------------------------------------------------------
void NGSParameters::success_parameter(){
    //------------- Check if the provided parameter values are appropriate -----------------//
    std::string temp_str;                                                                       // A temporary strin object to hold parameter name

    // If user wants to merge damages from multiple particles: 
    if(get_merge_damages_from_particles()){
        int reqSize = get_num_of_particles_to_merge();                                          // Get the number of particles to merge
        //check if an appropriate number of SDD files are provided
        if (reqSize != (get_sddfile_path()).size()){                                            // Exit with an error if number of SDD paths given not matches with reqired number
            std::cerr<<"\n ERROR: Mismatch between the number of particles to merge("<<reqSize<<") and the number of SDD files provided("<<(get_sddfile_path()).size()<<")"<<'\n';
            temp_str= "number_of_particles_to_merge"; help_parameter(&temp_str);                // Print help for number of particles to merge parameter
            temp_str= "sddFilePath"; help_parameter(&temp_str);                                 // Print help for SDD file path parameter
            exit(EXIT_FAILURE);
        }
        //check if an appropriate number of relative dose contributions are provided, if yes, check if their sum is 1
        if (reqSize != (get_relative_dose_contributions()).size()){                             // Exit with an error if number rel_dose_contributions given not matches with reqired number
            std::cerr<<"\n ERROR: Mismatch between the number of particles to merge("<<reqSize<<") and the number of relative contibutions provided("<<(get_relative_dose_contributions()).size()<<")"<<'\n';
            temp_str= "relative_dose_contributions"; help_parameter(&temp_str);                 // Print help for relative_dose_contributions parameter
            exit(EXIT_FAILURE);
        }else{
            double total_rel_dose{0};                                                           // Temporary variable to hold the total relative dose contribution
            for (int i=0;i<reqSize;i++){
                total_rel_dose += (get_relative_dose_contributions())[i];                       // Calculate the total relative dose contribution given
            }
            if (total_rel_dose != 1.0){                                                         // If total relative dose is not equal to 1, print error and exit
                std::cerr<<"\n ERROR: The sum of all relative dose contributions should be equal to 1 \n";
                temp_str= "relative_dose_contributions"; help_parameter(&temp_str);             // Print help for relative_dose_contributions parameter
                exit(EXIT_FAILURE);
            }  
        }
    }else{
        //If there are more than one SDD file path given, but the merge damage flag is off:
        if ((get_sddfile_path()).size()>1){
            std::cerr<<"\n ERROR: More than one SDD file paths are provided\n"
            <<" If you wish to merge damages from multiple SDD files, specify using \'merge_damages_from_multiple_particles parameter\'\n";
            exit(EXIT_FAILURE);
        }
        // If there are more than one relative dose contribution values given, but the merge damage flag is off:
        if ((get_relative_dose_contributions()).size()>1){
            std::string reset_rel_dose ="1";
            set_relative_dose_contributions(&reset_rel_dose);
            std::cerr<<"\n WARNING: More than 1 relative dose contribution values are provided for a single SDD file\n"
            <<"  ----- Re-setting \"relative_dose_contributions\" to 1 ----- \n";
        } 
    }

    // User MUST provide a valid path(s) to SDD file(s). If not, exit with error
    for (auto element : get_sddfile_path()){                                                     // Iterate through each element of the vector sddfile_path
        if (!checkFileExists(&element)){
            std::cerr<<"\n ERROR: Unable to read the SDD file : "<< element<<'\n';
            temp_str= "sddFilePath"; help_parameter(&temp_str);
            exit(EXIT_FAILURE);
        }else{
            std::cout<<"\n Successfully read the SDD file : "<< element<<'\n';
        }
    }

    // If user wants to adjust the damages according to the actual dose delivered in each exposure
    if(get_adjust_damages_with_actual_dose()){
        int reqSize = (get_sddfile_path()).size();                                              // Get the number of SDD files provided
        //check if an appropriate number of actual dose delivered files are provided
        if (reqSize != (*get_actual_dosefile_path()).size()){                                   // Exit with an error if number actual dose delivered files paths given not matches with reqired number
            std::cerr<<"\n ERROR: Mismatch between the number of particles to merge("<<reqSize<<") and the number of actual dose delivered files provided("<<(*get_actual_dosefile_path()).size()<<")"<<'\n';
            temp_str= "actual_dose_delivered_data"; help_parameter(&temp_str);                  // Print help for actual_dose_delivered_data parameter
            exit(EXIT_FAILURE);
        }
        // User MUST provide a valid path(s) to Actual_dose_delivered file(s). If not, exit with error
         for (auto element : *(get_actual_dosefile_path())){                                    // Iterate through each element of the vector actual_dosefile_path
            if (!checkFileExists(&element)){
                std::cerr<<"\n ERROR: Unable to read the Actual_dose_delivered file : "<< element<<'\n';
                temp_str= "actual_dose_delivered_data"; help_parameter(&temp_str);
                exit(EXIT_FAILURE);
            }else{
                std::cout<<"\n Successfully read the Actual_dose_delivered file : "<< element<<'\n';
            }
        }
    }

    // If the reference genome is user specified, then the filename should be valid. Else, exit with error
    if (!checkFileExists(get_reference_genome())){
        std::cerr<<"\n ERROR: Unable to read the reference genome file provided : "<< *get_reference_genome()<<'\n';
        temp_str= "reference_genome_FASTAfile"; help_parameter(&temp_str);
        exit(EXIT_FAILURE);
    }else{
        std::cout<<"\n Successfully read the reference genome file : "<< *get_reference_genome()<<'\n';
    }
    
    // Check if the user specified Illumina sequencer profile is in the list of built-in sequencer profiles
    auto it = std::find(list_sequencers.begin(), list_sequencers.end(), *get_sequencer());      // Find the location of the sequencer in the list
    if ( it == list_sequencers.end()) {
        std::cerr<<"\n ERROR: Illumina sequencer name provided : "<< *get_sequencer()<<" is not compatible \n";
        temp_str= "illumina_sequencer"; help_parameter(&temp_str);
        exit(EXIT_FAILURE);
    }else{
        std::cout<<"\n Successfully set the Illumina sequencer : "<< *get_sequencer()<<'\n';
        int index = std::distance(list_sequencers.begin(), it);                                 // Calculate the index of the sequencer name in the list
        set_read_length(index);                                                                 // Set the read length accordingly
        set_read_quality_profiles(index);                                                       // Set the filenames of read quality profiles accordingly
    }

    // Check if the number of cells to sequence is > number of cells in the sample. If true, exit with error
    if (get_num_of_cells_to_sequence() > get_num_of_cells_in_sample()){
        std::cerr<<"\n ERROR: The number of cells you wish to sequence ("<<get_num_of_cells_to_sequence()<<") is more than the total number of cells in the sample ("<<get_num_of_cells_in_sample()<<")\n";
        temp_str= "number_of_cells_to_sequence"; help_parameter(&temp_str);
        exit(EXIT_FAILURE);
    }

    // If user wants to do paired-end sequencing, make sure mean DNA fragment length and std dev are also provided. Else exit with error
    if (get_paired_end_sequencing()){
        if (get_mean_DNA_fragment_length()<get_read_length()){
            std::cerr<<"\n ERROR: The mean DNA fragment length ("<<get_mean_DNA_fragment_length()<<") provided is smaller than the read length ("<<get_read_length()<<")\n";
            temp_str= "mean_DNA_fragment_length"; help_parameter(&temp_str);
            exit(EXIT_FAILURE);
        }else if (!(get_std_dev_DNA_fragment_length() > 0)){
            std::cerr<<"\n ERROR: User must provide a non-zero positive standard deviation value\n";
            temp_str= "std_dev_DNA_fragment_length"; help_parameter(&temp_str);
            exit(EXIT_FAILURE);
        }
    }


}
//--------------------------------------------------------------------------------------------