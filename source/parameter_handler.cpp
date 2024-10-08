#include "parameter_handler.h"
#include "fileio.h"
#include "support_functions.h"
#include "random_generator.h"

#include <iostream>
#include <string>
#include <algorithm>

// Default constructor
NGSParameters::NGSParameters(){
    default_parameter_file_name = "/NGSDefaultParameters.txt";
    reference_genome_file_name = "/Human_reference_genome.fa";
    list_sequencers = {"HiSeq1000","HiSeq2000","HiSeq2500_v125","HiSeq2500_v150","HiSeqX","NovaSeq6000","Test","Custom"};
    list_read_lengths = {100, 100, 125, 150, 150, 151, 5};
    list_r1_quality_profiles = {"HiSeq1000_R1.txt","HiSeq2000_R1.txt","HiSeq2500_v125_R1.txt","HiSeq2500_v150_R1.txt","HiSeqX_R1.txt","NovaSeq6000_R1.txt","test_R1.txt"};
    list_r2_quality_profiles = {"HiSeq1000_R2.txt","HiSeq2000_R2.txt","HiSeq2500_v125_R2.txt","HiSeq2500_v150_R2.txt","HiSeqX_R2.txt","NovaSeq6000_R2.txt","test_R2.txt"};
}



//--------------------------------------------------------------------------------------------
// This function takes user-specified input file as argument and process the file
//--------------------------------------------------------------------------------------------
void NGSParameters::process_parameterFile(const std::string* parameterfile, NGSParameters& parameter, std::string* dataPath){
    // Parameter file is a MANDATORY requirement since we need to know SDD file path. If not found, exit with an error.
    if (!checkFileExists(parameterfile)){
        std::cerr<<"\n ERROR: Unable to read the parameter file : "<< *parameterfile<<'\n';
        exit(EXIT_FAILURE);
    }else{
        std::cout<<"\n Successfully read the parameter file: "<< *parameterfile<<'\n';
    }
    
    dataFolderPath = *dataPath;
    if(!checkFolderExists(dataFolderPath.c_str())){                // Check if the environment variable set by the user is pointing to the right directory
        std::cerr<<"\n ERROR: Unable to find the environment variable \"RADISEQ_DATA_DIR\" path correctly. The path provided is "<<dataFolderPath<<"\n";
        exit(EXIT_FAILURE);
    }
    // Defining data filenames w.r.t the variable dataFolderPath
    default_parameter_file = dataFolderPath+default_parameter_file_name;
    reference_genome_file = dataFolderPath+reference_genome_file_name;

    // Continue if parameter file is present and set all parameter values
    readParameterFile(&default_parameter_file, parameter);         // Read default parameter file first to set default parameter values
    readParameterFile(parameterfile, parameter);                   // Read user-specified parameter file to overwrite default parameter values

    success_parameter();                                           // Check if provided parameters meet all conditions and if it is appropriate

    int nThreads = get_number_of_threads();                        // Get the number of threads user is asking (default is 1)
    if(get_random_seed()){                                         // If user provided a non-zero seed value for the RNG (default is 0)
        rng::initThreadLocalFixedSeed(get_random_seed(), nThreads);// Initiate local instances of the RNG with the fixed seed
    }else{
        rng::initThreadLocalRandomSeed(nThreads);                  // Otherwise, initiate the RNG with a random seed
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
    else if (*paramName == "number_of_threads"){
        set_number_of_threads(paramValue);
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
    else if (*paramName == "custom_read1_quality_profile_path"){
        set_custom_r1_quality_profile_path(paramValue);
    }
    else if (*paramName == "custom_read2_quality_profile_path"){
        set_custom_r2_quality_profile_path(paramValue);
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
    else if (*paramName == "read_length"){
        set_read_length(paramValue);
    }
    else if (*paramName == "total_read_coverage"){
        set_total_read_coverage(paramName, paramValue);
    }
    else if (*paramName == "coverage_distribution"){
        set_coverage_distribution(paramName, paramValue);
    }
    else if (*paramName == "degree_of_GC_bias"){
        set_degree_of_GC_bias(paramName, paramValue);
    }
    else if (*paramName == "bin_size_for_GC_bias_estimation"){
        set_GC_binSize(paramValue);
    }
    else if (*paramName == "do_paired_end_sequencing"){
        set_paired_end_sequencing(paramName, paramValue);
    }
    else if (*paramName == "fragment_size_distribution_path"){
        set_fragment_size_distribution_path(paramName, paramValue);
    }
    else if (*paramName == "min_DNA_fragment_length"){
        set_min_DNA_fragment_length(paramValue);
    }
    else if (*paramName == "max_DNA_fragment_length"){
        set_max_DNA_fragment_length(paramValue);
    }
    else if (*paramName == "mode_DNA_fragment_length"){
        set_mode_DNA_fragment_length(paramValue);
    }
    else if (*paramName == "beta_of_beta_distribution"){
        set_beta_of_beta_distribution(paramName, paramValue);
    }
    else if (*paramName == "maximum_errors_in_reads"){
        set_max_errors_in_read(paramValue);
    }
    else if (*paramName == "max_fraction_unknown_bases_in_reads"){
        set_max_fraction_unknown_bases_in_reads(paramName, paramValue);
    }
    else if (*paramName == "fraction_of_other_oriented_read_pairs"){
        set_fraction_nonFR_read_pairs(paramName, paramValue);
    }
    else if (*paramName == "read_artifacts_rate"){
        set_read_artifacts_rate(paramName, paramValue);
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
void NGSParameters::set_number_of_threads(std::string* paramValue){
    if(std::stoi(*paramValue) > 0){                                                           // Set this parameter only if user specified a non-zero +ve value. Else default
        number_of_threads = atoi(paramValue->c_str());                                        // Converting the parameterValue to an int
    }
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
void NGSParameters::set_sddfile_path(std::string* paramValue){
    stringToVec(',', paramValue, sddfile_path);                                               // Registering comma seperated SDD paths to a vector
}
void NGSParameters::set_adjust_damages_with_actual_dose(std::string* paramName, std::string* paramValue){
   if(lowercaseString(paramValue) == "true"||lowercaseString(paramValue) == "false"){
        is_adjust_dose = (lowercaseString(paramValue) == "true");                             // Gets 1 if "True" or "true"; else 0
    }else{                                                                                    // If user provided unrecognized value format, print help
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
void NGSParameters::set_read_length(std::string* paramValue){                                 // This parameter is automatically set based on the user provided sequencer
    int tmp = std::stoi(*paramValue);                                                         // Converting the parameterValue to an int
    if(tmp > 0){                                                                              // set this parameter only if user specified a non-zero value. Else default
        read_length = tmp;
    }else{
        read_length = 0;                                                                      // Set to default
    }
}
void NGSParameters::set_read_length(int paramValue){                                          // Function overload for integer argument
    read_length = paramValue;
}
void NGSParameters::set_read_quality_profiles(int sequencer_index){
    r1_quality_profile = dataFolderPath+"/"+list_r1_quality_profiles[sequencer_index];
    r2_quality_profile = dataFolderPath+"/"+list_r2_quality_profiles[sequencer_index];
    if(!checkFileExists(&r1_quality_profile)){
        std::cerr<<"\n ERROR: Invalid RADISEQ_DATA_DIR or Missing files\n"
                 <<" The read quality files cannot be found at "<<dataFolderPath<<"\n";
        exit(EXIT_FAILURE);
    }
}
void NGSParameters::set_custom_r1_quality_profile_path(std::string* paramValue){
    path_to_custom_r1_quality_profile = *paramValue;
}
void NGSParameters::set_custom_r2_quality_profile_path(std::string* paramValue){
    path_to_custom_r2_quality_profile = *paramValue;
}
void NGSParameters::set_custom_read_quality_profiles(){
    r1_quality_profile = path_to_custom_r1_quality_profile;
    r2_quality_profile = path_to_custom_r2_quality_profile;
    if(!checkFileExists(&r1_quality_profile)){
        std::cerr<<"\n ERROR: Invalid file path or Missing files\n"
                 <<" The read quality file cannot be found at "<<path_to_custom_r1_quality_profile<<"\n";
        exit(EXIT_FAILURE);
    }
    if(!r2_quality_profile.empty() && !checkFileExists(&r2_quality_profile)){
        std::cerr<<"\n ERROR: Invalid file path or Missing files\n"
                 <<" The read quality file cannot be found at "<<path_to_custom_r2_quality_profile<<"\n";
        exit(EXIT_FAILURE);
    }
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
    if(std::stod(*paramValue) > 0){
        total_read_coverage = std::stod(*paramValue);
    }else{
        help_parameter(paramName);
    }
}
void NGSParameters::set_coverage_distribution(std::string* paramName, std::string* paramValue){
    if(lowercaseString(paramValue) == "uniform"||lowercaseString(paramValue) == "mda"){
        coverage_distribution = lowercaseString(paramValue);                                   // Gets 'uniform' or 'mda' depending on the user provided coverage distribution
    }else{                                                                                     // If user provided unrecognized value format, print help
        help_parameter(paramName);
        std::cerr<<" ----- Setting \""<<*paramName<<"\" to its default value: \"Uniform\" -----\n";
    }
}
void NGSParameters::set_degree_of_GC_bias(std::string* paramName, std::string* paramValue){
    try{
        degree_of_GC_bias = std::abs(std::stod(*paramValue));
    }catch(...){                                                                               // Error if provided value is not a number
        help_parameter(paramName);
    }
}
void NGSParameters::set_GC_binSize(std::string* paramValue){
    GC_binSize = std::stoi(*paramValue);
}
void NGSParameters::set_paired_end_sequencing(std::string* paramName, std::string* paramValue){
    if(lowercaseString(paramValue) == "true"||lowercaseString(paramValue) == "false"){
        is_paired_end_seq = (lowercaseString(paramValue) == "true");                           // Gets 1 if "True" or "true"; else 0
    }else{                                                                                     // If user provided unrecognized value format, print help
        help_parameter(paramName);
        std::cerr<<" ----- Setting \""<<*paramName<<"\" to its default value: \"True\" -----\n";
    }
}
void NGSParameters::set_fragment_size_distribution_path(std::string* paramName, std::string* paramValue){
    fragment_size_distribution_path = *paramValue;
}
void NGSParameters::set_is_fragment_distribution_from_file(){
    is_fragment_size_distribution = true;
}
void NGSParameters::set_min_DNA_fragment_length(std::string* paramValue){
    min_DNA_fragment_length = std::stoi(*paramValue);
}
void NGSParameters::set_max_DNA_fragment_length(std::string* paramValue){
    max_DNA_fragment_length = std::stoi(*paramValue);
}
void NGSParameters::set_mode_DNA_fragment_length(std::string* paramValue){
    mode_DNA_fragment_length = std::stod(*paramValue);
}
void NGSParameters::set_beta_of_beta_distribution(std::string* paramName, std::string* paramValue){
    if(std::stod(*paramValue) > 1.0){
        beta_of_beta_distribution = std::stod(*paramValue);
    }else{
        help_parameter(paramName);
        std::cerr<<" ----- Setting \""<<*paramName<<"\" to its default value: \"4.5\" -----\n";
    }
}
void NGSParameters::set_max_errors_in_read(std::string* paramValue){
    max_errors_in_read = std::stoi(*paramValue);
}
void NGSParameters::set_max_fraction_unknown_bases_in_reads(std::string* paramName, std::string* paramValue){
    if(std::stod(*paramValue) >= 0.0 && std::stod(*paramValue) <= 1.0){
        max_fraction_unknown_bases_in_reads = std::stod(*paramValue);
    }else{
        help_parameter(paramName);
        std::cerr<<" ----- Setting \""<<*paramName<<"\" to its default value: \"1.0\" -----\n";
    }
}
void NGSParameters::set_N_threshold_in_reads(int read_length, double max_fraction_unknown_bases_in_reads){
    N_threshold_in_reads = round(read_length*max_fraction_unknown_bases_in_reads);
}
void NGSParameters::set_fraction_nonFR_read_pairs(std::string* paramName, std::string* paramValue){
    if(std::stod(*paramValue) >= 0.0 && std::stod(*paramValue) <= 1.0){
        fraction_nonFR_read_pairs = std::stod(*paramValue);
    }else{
        help_parameter(paramName);
        std::cerr<<" ----- Setting \""<<*paramName<<"\" to its default value: \"0.0\" -----\n";
    }
}
void NGSParameters::set_read_artifacts_rate(std::string* paramName, std::string* paramValue){
    if(std::stod(*paramValue) >= 0.0 && std::stod(*paramValue) <= 1.0){
        read_artifacts_rate = std::stod(*paramValue);
    }else{
        help_parameter(paramName);
        std::cerr<<" ----- Setting \""<<*paramName<<"\" to its default value: \"0.0\" -----\n";
    }
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
int NGSParameters::get_number_of_threads(){
    return(number_of_threads);
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
const std::string* NGSParameters::get_custom_r1_quality_profile_path(){
    return(&path_to_custom_r1_quality_profile);
}
const std::string* NGSParameters::get_custom_r2_quality_profile_path(){
    return(&path_to_custom_r2_quality_profile);
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
double NGSParameters::get_total_read_coverage(){
    return(total_read_coverage);
}
const std::string* NGSParameters::get_coverage_distribution(){
    return(&coverage_distribution);
}
double NGSParameters::get_degree_of_GC_bias(){
    return(degree_of_GC_bias);
}
int NGSParameters::get_GC_binSize(){
    return(GC_binSize);
}
bool NGSParameters::get_paired_end_sequencing(){
    return(is_paired_end_seq);
}
const std::string* NGSParameters::get_fragment_size_distribution_path(){
    return(&fragment_size_distribution_path);
}
bool NGSParameters::get_is_fragment_distribution_from_file(){
    return(is_fragment_size_distribution);
}
int NGSParameters::get_min_DNA_fragment_length(){
    return(min_DNA_fragment_length);
}
int NGSParameters::get_max_DNA_fragment_length(){
    return(max_DNA_fragment_length);
}
double NGSParameters::get_mode_DNA_fragment_length(){
    return(mode_DNA_fragment_length);
}
double NGSParameters::get_beta_of_beta_distribution(){
    return(beta_of_beta_distribution);
}
int NGSParameters::get_max_errors_in_read(){
    return(max_errors_in_read);
}
double NGSParameters::get_max_fraction_unknown_bases_in_reads(){
    return(max_fraction_unknown_bases_in_reads);
}
int NGSParameters::get_N_threshold_in_reads(){
    return(N_threshold_in_reads);
}
double NGSParameters::get_fraction_nonFR_read_pairs(){
    return(fraction_nonFR_read_pairs);
}
double NGSParameters::get_read_artifacts_rate(){
    return(read_artifacts_rate);
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
        <<" Built-in sequencers: HiSeq1000, HiSeq2000, HiSeq2500_v125, HiSeq2500_v150, HiSeqX, NovaSeq6000 and Custom\n";
    }
    else if (*paramName == "custom_r1_quality_profile"){
        std::cerr<<" Specify the path to the file containing read 1 quality profile. \n"
        <<" This is a required parameter when the custom sequencer is chosen. \n";
    }
    else if (*paramName == "custom_r2_quality_profile"){
        std::cerr<<" Specify the path to the file containing read 2 quality profile. \n"
        <<" When paired-end sequencing is needed with a custom sequencer, both read 1 and read 2 quality profiles MUST also be specified. \n";
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
        std::cerr<<" Specify the total read coverage you want to get from this sequencing (Double value expected). \n"
        <<" If single-cell sequencing, the coverage gets distributed over the number of cells to be sequenced. \n"
        <<" i.e, if 100 cells sequenced with a coverage of 100x each cell will have 1x coverage in single-cell sequencing \n";
    }
    else if (*paramName == "coverage_distribution"){
        std::cerr<<" Specify whether you want to use the read coverage distribution of MDA or not for single-cell sequencing. \n"
        <<" This parameter can take only values 'uniform' and 'mda' for two different coverage distribution types. \n"
        <<" This option to choose the WGA distribution is only available for single-cell sequencing \n";
    }
    else if (*paramName == "degree_of_GC_bias"){
        std::cerr<<" Specify the slope of the lines that make a triangular function for GC bias model. \n"
        <<" A double value is expected. And the same slope will be used for the +ve and -ve lines. \n";
    }
    else if (*paramName =="bin_size_for_GC_bias_estimation"){
        std::cerr<<" Specify a non-zero integer value for the bin size. GC bias will be estimated over these bins.\n"
        <<" Note that the smaller the bin size, the longer the simulation will take to run\n";
    }
    else if (*paramName == "do_paired_end_sequencing"){
        std::cerr<<" This parameter should be set \"True\" or \"False\" "
        <<"to specify whether or not you wish to perform paired-end sequencing \n";
    }
    else if (*paramName == "fragment_size_distribution_path"){
        std::cerr<<" Specify the path to the file containing the DNA fragment size distribution\n"
        <<" Each row in this file is expected to be space-seperated value pairs corresponding to the DNA fragment size and the normalized count of that fragment respectively\n"
        <<" Alternatively, the DNA fragment size distribution can be described using a beta function parameters. However, if fragment_size_distribution_path is set, \n"
        <<" then the beta function will be overriden with the provided fragment size distribution.\n";
    }
    else if (*paramName == "DNA_fragment_length"){
        std::cerr<<" Specify the minimum and maximum lengths of DNA fragments (in bp) obtained after size selection. This is different from the read length. \n"
        <<" Make sure the minimum and maximum values correspond to the lower and upper limit of the desired distribution respectively\n"
        <<" Alternatively, a file containing fragment size distribution can be provided using the parameter \'fragment_size_distribution_path\' \n";
    }
    else if (*paramName == "mode_DNA_fragment_length"){
        std::cerr<<" Specify the mode DNA fragment size for the desired DNA fragment size distribution. \n"
        <<" This value should be between the minimum and maximum DNA fragment sizes \n";
    }
    else if (*paramName == "beta_of_beta_distribution"){
        std::cerr<<" The beta parameter for the beta distribution needs to be greater than 1 \n";
    }
    /* else if (*paramName == "maximum_errors_in_reads"){
        std::cerr<<" The maximum number of indel errors you can have in a read (Integer value expected)\n"
        <<" If not specified, it will default to no restrictions on the number of indels\n";
    } */
    else if (*paramName == "max_fraction_unknown_bases_in_reads"){
        std::cerr<<" Specify the allowable fraction of bases in a read to be unknown (N's). \n"
        <<" This value should be a double in the range [0,1] \n";
    }
    else if (*paramName == "fraction_of_other_oriented_read_pairs"){
        std::cerr<<" Specify what fraction of the total read pairs needs to be not in inward-orientation in paired-end sequencing. \n"
        <<" This value should be a double in the range [0,1]\n";
    }
    else if (*paramName == "read_artifacts_rate"){
        std::cerr<<" Specify the read artifacts formation rate for both read 1 and 2 combined. \n"
        <<" This value should be a double in the range [0,1]\n";
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
        size_t reqSize = static_cast<size_t>(get_num_of_particles_to_merge());                  // Get the number of particles to merge
        //check if an appropriate number of SDD files are provided
        if (reqSize != (get_sddfile_path()).size()){                                            // Exit with an error if number of SDD paths given not matches with reqired number
            std::cerr<<"\n ERROR: Mismatch between the number of particles to merge("<<reqSize<<") and the number of SDD files provided("<<(get_sddfile_path()).size()<<")"<<'\n';
            temp_str= "number_of_particles_to_merge"; help_parameter(&temp_str);                // Print help for number of particles to merge parameter
            temp_str= "sddFilePath"; help_parameter(&temp_str);                                 // Print help for SDD file path parameter
            exit(EXIT_FAILURE);
        }
    }else{
        //If there are more than one SDD file path given, but the merge damage flag is off:
        if ((get_sddfile_path()).size()>1){
            std::cerr<<"\n ERROR: More than one SDD file paths are provided\n"
            <<" If you wish to merge damages from multiple SDD files, specify using \'merge_damages_from_multiple_particles parameter\'\n";
            exit(EXIT_FAILURE);
        }
    }

    // User MUST provide a valid path(s) to SDD file(s). If not, exit with error
    for (auto element : get_sddfile_path()){                                                    // Iterate through each element of the vector sddfile_path
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
        size_t reqSize = (get_sddfile_path()).size();                                           // Get the number of SDD files provided
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
    
    // Check if the user provided value for the maximum errors in a read is more than the read length or if it is less than -1 (default value)
    /* if (get_max_errors_in_read()>get_read_length() || get_max_errors_in_read()<-2){
        std::cerr<<"\n WARNING: The value provided for the maximum number of indel errors in a read ("<<get_max_errors_in_read()<<") is invalid. Make sure the value is not exceeding the read length\n"
        <<"  ----- Re-setting \"maximum_errors_in_reads\" to default: no limits ----- \n";
        std::string reset_max_errors ="-1";
        set_max_errors_in_read(&reset_max_errors);
    } */
    
    // Check if the user specified Illumina sequencer profile is in the list of built-in sequencer profiles
    auto it = std::find_if(list_sequencers.begin(),list_sequencers.end(), [this](const std::string& s){return compareStrings(&s, get_sequencer());});
    //auto it = std::find(list_sequencers.begin(), list_sequencers.end(), *get_sequencer());      // Find the location of the sequencer in the list
    std::vector<int>::size_type index;
    if ( it == list_sequencers.end()){
        std::cerr<<"\n ERROR: Illumina sequencer name provided : "<< *get_sequencer()<<" is not compatible \n";
        temp_str = "illumina_sequencer"; help_parameter(&temp_str);
        exit(EXIT_FAILURE);
    }else{
        std::cout<<"\n Successfully set the Illumina sequencer : "<< *get_sequencer()<<'\n';
        index = std::distance(list_sequencers.begin(), it);                                     // Calculate the index of the sequencer name in the list
    }

    // Check if the user provided read length is appropriate
    if (index < list_read_lengths.size()){                                                      // If a built-in sequencer is chosen
        if (get_read_length() == 0 || get_read_length() > list_read_lengths[index]){
            std::cerr<<"\n WARNING: The maximum allowable read length is "<<list_read_lengths[index]<<" bp and the value provided is inappropirate \n"
            <<" Therefore, the maximum length "<<list_read_lengths[index]<<" will be used (default) in this run\n";
            set_read_length(list_read_lengths[index]);
        }
        set_read_quality_profiles(index);                                                       // Set the filenames of read quality profiles accordingly
    }else{                                                                                      // If the custom sequencer is selected
        std::string empty_path = "\"\"";                                                        // Temporary string to hold empty path, which is the default value
        if(get_paired_end_sequencing() && *get_custom_r2_quality_profile_path()== empty_path){  // If paired-end sequencing is needed, read 2 quality profile must also be provided 
            std::cerr<<"\n ERROR: The path to the read 2 quality profile MUST be provided with the custom sequencer \n";
            temp_str = "custom_r2_quality_profile"; help_parameter(&temp_str);
            exit(EXIT_FAILURE);
        }
        if(*get_custom_r1_quality_profile_path() == empty_path){                                // User must provide the read 1 quality profile for the custom sequencer
            std::cerr<<"\n ERROR: The path to the read 1 quality profile MUST be provided with the custom sequencer \n";
            temp_str = "custom_r1_quality_profile"; help_parameter(&temp_str);
            exit(EXIT_FAILURE);
        }
        set_custom_read_quality_profiles();
        if(get_read_length() == 0){
            std::cerr<<"\n ERROR: A non-zero read length must be provided with the custom sequencer \n";
            exit(EXIT_FAILURE);
        }
    }
    //set_N_threshold_in_reads(get_read_length(), get_max_fraction_unknown_bases_in_reads());     // Set the threshold on the number of unknown bases in a read

    // Check if the number of cells to sequence is > number of cells in the sample. If true, exit with error
    if (get_num_of_cells_to_sequence() > get_num_of_cells_in_sample()){
        std::cerr<<"\n ERROR: The number of cells you wish to sequence ("<<get_num_of_cells_to_sequence()<<") is more than the total number of cells in the sample ("<<get_num_of_cells_in_sample()<<")\n";
        temp_str= "number_of_cells_to_sequence"; help_parameter(&temp_str);
        exit(EXIT_FAILURE);
    }

    // Check if MDA read coverage is set for bulk-cell sequencing and not single-cell sequencing
    if (*get_coverage_distribution()=="mda"){
        if(*get_sequencing_mode()=="bulk"){
            std::cerr<<"\n WARNING: MDA distribution is not available for bulk-cell sequencing\n"
            <<"  ----- Re-setting \"coverage_ditribution\" to \'Uniform\' ----- \n";
            std::string paramValue = "uniform"; std::string paramName = "coverage_distribution"; 
            set_coverage_distribution(&paramName, &paramValue);
        }
    }

    // Check if GC bias is set for single-cell sequencing and not bulk-cell sequencing
    /* if (get_degree_of_GC_bias()!=0.0){
        if(*get_sequencing_mode()=="single"){
            std::cerr<<"\n WARNING: GC bias option is not available for single-cell sequencing\n"
            <<"  ----- Re-setting \"degree of GC bias\" to \'0.0\' ----- \n";
            std::string paramValue = "0.0"; std::string paramName = "degree_of_GC_bias"; 
            set_degree_of_GC_bias(&paramName, &paramValue);
        }
    } */
    // Check if the bin size provided for calculating GC bias is zero. 
    if (get_degree_of_GC_bias()!=0.0){
        if(get_GC_binSize() == 0){
            temp_str="bin_size_for_GC_bias_estimation"; help_parameter(&temp_str);
            exit(EXIT_FAILURE);
        }
    }


    // If user wants to do paired-end sequencing, make sure the bounds of the DNA fragment length distributions are appropriate. Else exit with error
    if (get_paired_end_sequencing()){
        std::string empty_path = "\"\"";                                                        // Temporary string to hold empty path, which is the default value
        if(*get_fragment_size_distribution_path()!=empty_path){                                 // If the fragment size distribution from a file should be used, then check the file
            if (!checkFileExists(get_fragment_size_distribution_path())){
                std::cerr<<"\n ERROR: Unable to read the fragment size distribution file provided : "<< *get_fragment_size_distribution_path()<<'\n';
                temp_str= "fragment_size_distribution_path"; help_parameter(&temp_str);
                exit(EXIT_FAILURE);
            }else{
                std::cout<<"\n Successfully read the fragment size distribution file : "<< *get_fragment_size_distribution_path()<<'\n';
                set_is_fragment_distribution_from_file();
            }
        }else{                                                                                 // If the beta function should be used to generate the distribution, then check the beta parameters. 
            if (get_min_DNA_fragment_length()>=get_max_DNA_fragment_length()){
                std::cerr<<"\n ERROR: The maximum DNA fragment length ("<<get_max_DNA_fragment_length()<<" bp) provided is equal/smaller than the minimum DNA fragment length ("<<get_min_DNA_fragment_length()<<" bp)\n";
                temp_str= "DNA_fragment_length"; help_parameter(&temp_str);
                exit(EXIT_FAILURE);
            }
            if (get_min_DNA_fragment_length()>=get_mode_DNA_fragment_length() || get_max_DNA_fragment_length()<=get_mode_DNA_fragment_length()){
                std::cerr<<"\n ERROR: The mode DNA fragment length ("<<get_max_DNA_fragment_length()<<") provided is not between the lower and upper bounds given ( ["<<get_min_DNA_fragment_length()<<","<<get_max_DNA_fragment_length()<<"] )\n";
                temp_str= "mode_DNA_fragment_length"; help_parameter(&temp_str);
                exit(EXIT_FAILURE);
            }
        }
    }


}
//--------------------------------------------------------------------------------------------