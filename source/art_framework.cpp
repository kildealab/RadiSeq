#include "art_framework.h"
#include "random_generator.h"
#include "fileio.h"
#include "fastafile_handler.h"
#include "support_functions.h"

#include <vector>
#include <iostream>
#include <map>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <bitset>



//--------------------------------------------------------------------------------------------
// Definition of all static variables that can be accessed by all ART class objects
//--------------------------------------------------------------------------------------------
std::vector<double> ART::read_pair_orientation_dist;                                                    // Vector that stores the fraction of read pair orientation for paired end sequencing
std::vector<double> ART::read_artifacts_rate_dist;                                                      // Vector that stores the rate of read artifacts and no artifact for a read generation
std::vector<std::map<unsigned int, unsigned short>> ART::read1_quality_distribution_vec;                // Quality distribution vector for read 1. This vector will not change for each cell 
std::vector<std::map<unsigned int, unsigned short>> ART::read2_quality_distribution_vec;                // Quality distribution vector for read 2. This vector will not change for each cell 
double ART::baseCall_error_probability[80];                                                             // Array that will hold 80 values of base call error probabilities.
int ART::read_length{0};
std::string* ART::chromSegSeq;
long ART::valid_region{0};
long ART::num_initial_primer_sites{0};                                                                  // Variable to hold the number of initial primer sites in a chromSegment for MDA
std::vector<long> ART::initial_MDAsites;                                                                // Temporary vector to hold the initial primer binding locations for MDA
std::vector<double> ART::initial_MDAsites_bias;                                                         // Temporary vector to hold the bias associated with each of the initial primer binding location for MDA
std::string ART::wga_technique{"uniform"};                                                              // Variable to hold the value that indicate the WGA technique chosen. Default uniform amplification
int ART::GC_binSize{10000};                                                                             // Variable to hold the bin size over which the GC fraction and bias needds to be calculated. Default is 10K.
//std::unordered_set<long> ART::forbidden_locations;                                                      // A set that stores all the locations in the sequence segment where a read should not be created    
std::unordered_set<long> ART::forbidden_bins;                                                           // A set that stores all the GC bins in the sequence segment where a read should not be created
//std::vector<double> ART::GCsites_bias;                                                                  // Vector to hold the GC bias of each location in a chromosome segment
std::vector<double> ART::GCbins_bias;                                                                   // Vector to hold the GC bias of each bin in a chromosome segment
//std::vector<double> ART::cum_GCsites_bias;                                                              // Vector to hold the cumilative bias upto each location in the chromosome segment
std::vector<double> ART::cum_GCbins_bias;                                                               // Vector to hold the cumilative bias upto each bin in the chromosome segment
//--------------------------------------------------------------------------------------------



// Default constructor
ART::ART(){}
    


//--------------------------------------------------------------------------------------------
// This function will initialize all the constant parameters that dictates how the read 
// generation engine should work. These parameters are read length, GC bin size, fraction of 
// outward-oriented (RF) read pairs and base call error probability. This function takes in
// all these parameters as arguments except for base call error probability, which is generated
// with a function call within.
//--------------------------------------------------------------------------------------------
void ART::initiate_read_generation(int read_len, int User_GC_binSize, double fraction_nonFR_readPairs, double read_artifacts_rate){
    read_length = read_len;
    GC_binSize = User_GC_binSize;
    read_pair_orientation_dist.clear();
    read_artifacts_rate_dist.clear();
    double fraction_inward_readPairs = 1.0-fraction_nonFR_readPairs;                                    // Fraction of inward-oriented(FR) read pairs is assumed to be all pairs other than outward-oriented(RF) ones
    double fraction_outward_readPairs = fraction_nonFR_readPairs * 2/3;                                 // 2/3 rd of the other read orientation is expected to be RF (outward orientation)
    double fraction_other_readPairs = fraction_nonFR_readPairs - fraction_outward_readPairs;            // Whatever is not outward and inward will be used to generate FF and RR orientations
    for(double val : {fraction_inward_readPairs, fraction_outward_readPairs, fraction_other_readPairs}){// Push these read orientation fractions to the vector for later use. 
        read_pair_orientation_dist.push_back(val);
    }
    for(double val : {read_artifacts_rate, 1.0-read_artifacts_rate}){
        read_artifacts_rate_dist.push_back(val);
    }
    set_baseCall_error_probability();                                                                   // Call the base call error probability setting function
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function sets the base call error probability,that the program will later use to add
// base call errors to the read. The quality scores Q are logarithmically related to the 
// base-calling error probabilities P with the equation: P = pow(10, -Q/10). However, at this
// point, we do not know the value of Q. Therefore we generate an array which can be used later
// to find the base call error probability for a specific Q. 
//--------------------------------------------------------------------------------------------
void ART::set_baseCall_error_probability(){
    for(int i=0; i<80; i++){                                                                            // 80 is the size of the baseCall_error_probability array
        baseCall_error_probability[i] = pow(10,-i/10.0);
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function will set the read quality distribution for read 1 and read 2 based on the 
// read quality profile files passed. 
//--------------------------------------------------------------------------------------------
void ART::set_read_quality_distribution(const std::string& r1_quality_filename,  const std::string& r2_quality_filename){
    std::ifstream r1_quality_file(r1_quality_filename.c_str());
    make_quality_distribution(r1_quality_file, read1_quality_distribution_vec);
    
    std::ifstream r2_quality_file(r2_quality_filename.c_str());
    make_quality_distribution(r2_quality_file, read2_quality_distribution_vec);
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function calculates the probability of having atleast x errors (insertions/deletions) 
// in a read of length read_length if the corresponding sequencing error rate is specified.
// x can take all values from 0 to read_length. The calculated probabilities are stored as 
// an array into the vector object: error_probability_vec that is passed.
//--------------------------------------------------------------------------------------------
void ART::set_read_error_probability(int read_len, double error_rate, std::vector<double>& error_probability_vec, int max_errors){
    error_probability_vec.clear();                                                                      // Empty the vector 
    double cdf_cutoff = 0.999999;                                                                       // Threshold value for total probability
    
    if(max_errors == 0) return;                                                                         // If max errors is set to 0 (no errors), then return

    if(error_rate<1e-30) return;                                                                        // If user-specified error rate is < 10^-30, set it 0
    
    // Probability that there are no errors in a read of length read_length for a sequencing error rate given
    double probability = rng::binomial_cdf(0, error_rate, read_length);
    double total_probability = probability;
    
    // Determine the probability that there are atleast i errors in a read of length read_length for the sequencing error_rate
    for(int i=1; i<read_length+1; i++){
        probability = rng::binomial_cdf(i, error_rate, read_length);
        error_probability_vec.push_back(probability);
        if(max_errors>0 && (i>=max_errors)) break;                                                      // Return when user-specified limit of errors reached
        total_probability += probability; 
        if(total_probability >= cdf_cutoff) break;                                                      // If total probability exceeded the threshold, return
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
void ART::set_read_error_rates(double ins_rate, double del_rate){
    insertion_rate = ins_rate;
    deletion_rate = del_rate;
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function will take in read length and the sequence of the current chromosome segment
// and sets the static variables read_lenth and valid_region along with the chromSegmentSeq.
// Returns true if the chromosome segment is big enough to make reads.
// Function takes two optional arguments, which are essential if the WGA technique is MDA
// to generate phi29 binding sites in the chromosome segment
//--------------------------------------------------------------------------------------------
bool ART::int_set(std::string& chromSeq, int nThreads){
    valid_region = chromSeq.size() - read_length;
    if(valid_region<1) return false;                                                                    // If read length is >= chromSegmentSeq length, then do not proceed with the chromSegmentSeq
    chromSegSeq = &chromSeq;
    indel_map_vec.clear();                                                                              // Clear any pre-existing data
    read_seq_vec.clear();
    indel_map_vec.resize(nThreads);                                                                     // Resize the vector to accomodate all threads
    read_seq_vec.resize(nThreads);
    return true;
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function returns the chromosome segment sequence that is currently being processed. 
//--------------------------------------------------------------------------------------------
/* const std::string* ART::get_chromSegmentSeq(){
    return (chromSegSeq);
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function returns the read length that we want to generate
//--------------------------------------------------------------------------------------------
int ART::get_read_length(){
    return (read_length);
} */
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function identifies all the locations in a chromosome segment from where a read should
// not be created. This analysis is performed after considering the user input on the maximum
// allowed unknown (N's) bases in a read and all the locations that would exceed this limit
// if used to create a read is added to the list of forbidden locations to avoid. Function 
// returns TRUE if succesfull (if the segment has atleast one good location)
// This function also calculates the GC bias of each location in a chromosome segment and
// populate the vector GCsites_bias. 
//--------------------------------------------------------------------------------------------
/* bool ART::get_N_mask_and_GCbias(int N_threshold, int minFragSize){
    forbidden_locations.clear();                                                                        // Reset the list of forbidden locations
    GCsites_bias.clear();                                                                               // Reset the GC bias vector for each chrom segment
    cum_GCsites_bias.clear();                                                                           // Reset the cumulative GC sites bias vector for each chrom segment
    const size_t chromSeqLength = chromSegSeq->length();                                                // Pre-calculate the chromosome length and store it to avoid repeated calculations
    if(chromSeqLength < static_cast<size_t>(read_length)) return false;                                 // Return false if chromosome is smaller than read length
    double GCslope = GCBias::get_GCbias_slope();                                                        // Temporary variable to hold the GC biase slope to avoid calling the function again and again
    if(N_threshold == read_length && GCslope == 0.0) return true;                                       // If number of acceptable N's is equal to the read length, then there is no need of masking. Continue with empty forbidden list
    std::bitset<256> isN;                                                                               // Define bitset to do bit manipulation of ASCII characters
    std::bitset<256> isGC;
    isN.set('N');                                                                                       // Set the character for bitset to be 'N'
    isGC.set('G');                                                                                      // Set the character bit to be of 'G'              
    isGC.set('C');                                                                                      // Set the character bit to be also of 'C'
    int num_Ns{0};                                                                                      // A temporary variable to hold the count of N's in a window of read_length
    int GC_count{0};                                                                                    // A temporary variable to hold the count of G and C in a window of read_length
    double total_GCbias{0};                                                                             // Temporary variable to hold the total GC bias value for the chromosome sequence
    for(int i=0; i<read_length-1; i++){       
        char currentChar = (*chromSegSeq)[i];                                                           // Copy the character value to the currentChar
        if(isN.test(static_cast<unsigned char>(currentChar))) num_Ns++;                                 // Count the number of Ns in the first (read-2) base locations
        else if(isGC.test(static_cast<unsigned char>(currentChar))) GC_count++;                         // Count the number of Gs and Cs in the first (read-2) base locations
    }
    
    // Sliding window approach to count the number of N's and GCs in a region that is read_length size long
    std::vector<int> indices_to_insert;                                                                 // Accumulator to collect indices to insert into forbidden_locations
    for(int j=0; j<valid_region+1; j++){                      
        char rightChar = (*chromSegSeq)[j+read_length-1];                                               // Character to the right of the sliding window
        if(isN.test(static_cast<unsigned char>(rightChar))) num_Ns++;                                   // Check if the next character to the right is an N. If yes, increase the count
        else if(isGC.test(static_cast<unsigned char>(rightChar))) GC_count++;                           // Check if the next character to the right is either a G or a C. If yes, increase the count
        
        if(num_Ns >= N_threshold) indices_to_insert.push_back(j);                                       // Check if the number of N's exceeded the threshold. If yes, add the starting location of the window to the list           

        if(GCslope != 0.0){                                                                             // Do the following only if degree of GC bias is not zero
            double GC_bias = GCBias::get_GCbias(static_cast<double>(GC_count)/(read_length-num_Ns));    // Calculate the GC bias value for each location in the chromosome segment. GC fraction must be obtained by counting non-N bases
            if (GC_bias==0.0)indices_to_insert.push_back(j);                                            // If current location's GC bias is zero, make it a forbidden location 
            GCsites_bias.push_back(GC_bias);                                                            // Make an array of all site biases
            total_GCbias += GC_bias;  
        }
        
        char leftChar = (*chromSegSeq)[j];                                                              // First character to the left of the sliding window
        if(isN.test(static_cast<unsigned char>(leftChar))) num_Ns--;                                    // Check if the first character on the left is an N. If yes, reduce its count for the next window
        else if(isGC.test(static_cast<unsigned char>(leftChar))) GC_count--;                            // Check if the first character on the left is either a G or a C. If yes, reduce its count for the next window
    }
    
    for (int index : indices_to_insert){                                                                // Insert collected indices into forbidden_locations. This will avoid the overhead of inserting each element seperately.
        forbidden_locations.insert(index);                                                              // Since forbidden_locations is an unordered set, any duplicates will not be stored. 
    }

    if(total_GCbias != 0){
        double cum_GCbias{0};                                                                           // Temporary variable to hold the cumulative GC bias upto each location of the chromosome sequence
        for (double& siteBias : GCsites_bias){                                                          // Normalize the GC bias vector for the chromosome segment
            siteBias /= total_GCbias;
            cum_GCbias += siteBias;
            cum_GCsites_bias.push_back(cum_GCbias);
        }
    }else if(total_GCbias == 0 && GCslope != 0.0) return false;                                         // If GC bias in the valid region is zero, then the segment is not usable

    if(N_threshold == read_length) return true;                                                         // If number of acceptable N's is equal to the read length, then there is no need of masking. Continue with empty forbidden list
    if(minFragSize >= (static_cast<int>(chromSeqLength))){                                              // When chromSeg size is less than the minFragment size, the generated fragment will have chromSeq size. Which means the read will start always from location 0
        if(forbidden_locations.find(0)!=forbidden_locations.end()) return false;                        // If location 0 is forbidden, then the chromosome segment is not usable
    }else{
        if(cum_GCsites_bias.size() !=0 && cum_GCsites_bias[static_cast<int>(chromSeqLength)-minFragSize] == 0.0) return false;// If cumulative bias vector is not empty but the cumulative bias upto the available region is zero, segment is not usable
    }
    return static_cast<long int>(forbidden_locations.size()) != valid_region;                           // If the entire chromosome segment is not good for generating any reads, exit false else true
}
//-------------------------------------------------------------------------------------------- */



//--------------------------------------------------------------------------------------------
// This function calculates the GC bias in a chromosome segment for a bin size specified and
// populate the vector GCbinss_bias. 
//--------------------------------------------------------------------------------------------
bool ART::GCbias_maker(int minFragSize){
    forbidden_bins.clear();                                                                             // Reset the list of forbidden bins
    GCbins_bias.clear();                                                                                // Reset the GC bias vector for each chrom segment
    cum_GCbins_bias.clear();                                                                            // Reset the cumulative GC sites bias vector for each chrom segment
    const size_t chromSeqLength = chromSegSeq->length();                                                // Pre-calculate the chromosome length and store it to avoid repeated calculations
    if(chromSeqLength < static_cast<size_t>(read_length)) return false;                                 // Return false if chromosome is smaller than read length
    double GCslope = GCBias::get_GCbias_slope();                                                        // Temporary variable to hold the GC biase slope to avoid calling the function again and again
    if(GCslope == 0.0) return true;                                                                     // If the degree of GC bias is zero, no need to calculate further. Continue with empty forbidden list
    std::bitset<256> isN;
    std::bitset<256> isGC;
    isN.set('N');                                                                                       // Set the character bit to be of 'N'  
    isGC.set('G');                                                                                      // Set the character bit to be of 'G'              
    isGC.set('C');                                                                                      // Set the character bit to be also of 'C'
    int GC_count{0};                                                                                    // A temporary variable to hold the count of G and C in a bin
    int N_count{0};                                                                                     // A temporary variable to hold the count of N in a bin
    double total_GCbias{0};                                                                             // Temporary variable to hold the total GC bias value for the chromosome sequence
    int num_GC_bins{0};                                                                                 // Temporary variable to hold the number of GC bins in the segment. (absolute count not indices)
    std::vector<int> binID_to_insert;                                                                   // Temporary vector that stores the ID of bins that needs to be added to the forbidden bins list

    for(size_t i=0; i<chromSeqLength; i++){
        char currentChar = (*chromSegSeq)[i]; 
        if(isGC.test(static_cast<unsigned char>(currentChar))){GC_count++;}                             // Check if the current character is either G or C. If yes, count them
        else if(isN.test(static_cast<unsigned char>(currentChar))){N_count++;}                          // Check if the current character is N. If yes, count them

        if((i+1)%GC_binSize==0){                                                                        // If i+1 index is a multiple of the GC_binSize, then it is the end of a bin and perform the following calculations
            double GC_binBias{0.0};
            if((GC_binSize-N_count)!=0){
                GC_binBias = GCBias::get_GCbias(static_cast<double>(GC_count)/(GC_binSize-N_count));    // Calculate the GC bias value for each bin in the chromosome segment. GC fraction must be obtained by counting non-N bases
            }
            GCbins_bias.push_back(GC_binBias);                                                          // Make an array of all bin biases
            total_GCbias += GC_binBias;                                                                 // Calculate the total bin biases for normalization
            num_GC_bins++; GC_count = 0; N_count = 0;                                                   // Count the bin and reset other counters
            if (GC_binBias==0.0)binID_to_insert.push_back(num_GC_bins);                                 // If current bin's GC bias is zero, make it a forbidden bin 
        }
    }
    int unprocessed_segment = chromSeqLength-(num_GC_bins*GC_binSize);                                  // Find if there are xtra sequence after the full bins, or if the sequence is smaller than a bin size
    if (unprocessed_segment>=read_length){                                                              // Calculate the bias for this xtra bit only if it's length is equal to or greater than the read length
        double GC_bias{0.0};
        if((unprocessed_segment-N_count)!=0){                                                           // Calculate the GC bias value for the extra sequence. GC fraction must be obtained by counting non-N bases. Else the bias will be 0
            GC_bias = GCBias::get_GCbias(static_cast<double>(GC_count)/(unprocessed_segment-N_count));
        }
        GCbins_bias.push_back(GC_bias);                                                                 // Add this to the bin bias
        total_GCbias += GC_bias;
        num_GC_bins++;                                                                                  // Count this also as a bin since the GC fraction is independent of bin size
        if (GC_bias==0.0)binID_to_insert.push_back(num_GC_bins);                                        // Add the bin to forbidden bin list if bias is zero
    }
    if(total_GCbias == 0) return false;                                                                 // Return if this segment is not useful

    for (int binID : binID_to_insert){                                                                  // Insert collected bin IDs into forbidden_bins list. This will avoid the overhead of inserting each element seperately.
        forbidden_bins.insert(binID);                                                                   // Since forbidden_bins is an unordered set, any duplicates will not be stored. 
    }

    double cum_GCbias{0};                                                                               // Temporary variable to hold the cumulative GC bias upto each GC bin in the chromosome sequence
    for (double& binBias : GCbins_bias){                                                                // Normalize the GC bias vector for the chromosome segment
        binBias /= total_GCbias;
        cum_GCbias += binBias;                                                                          // Calcuate the cumulative GC bias of each bin for later use
        cum_GCbins_bias.push_back(cum_GCbias);
    }
    if(minFragSize >= (static_cast<int>(chromSeqLength))){                                              // When chromSeg size is less than the minFragment size, the generated fragment will have chromSeq size. Which means the read will start always from bin 1
        if(forbidden_bins.find(1)!=forbidden_bins.end()) return false;                                  // If bin 1 is forbidden, then the chromosome segment is not usable
    }else{
        int available_bins = (floor(static_cast<int>(chromSeqLength)-minFragSize)/GC_binSize)+1;        // Total number of bins in the available region if the bin size if GC-binSize. +1 to get the bin ID (not index)
        if(cum_GCbins_bias.size() !=0 && cum_GCbins_bias[available_bins-1] == 0.0) return false;        // If cumulative bias vector is not empty but the cumulative bias upto the available region is zero, segment is not usable. Bin index is used
    }
    return static_cast<int>(forbidden_bins.size()) != num_GC_bins;                                      // If all the bins are forbidden in the segment, return false coz the segment is not usable. Else return true.  
}
//--------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------
// This function determines where should be the initial primer binding sites for MDA and what
// should be the bias of each location to have a read generated from there. It is essentially
// generating a non-uniformity for read distribution according to MDA kinetics. Function returns
// TRUE if successful (If reads need to be generated from this segment)
//--------------------------------------------------------------------------------------------
bool ART::MDA_distribution_maker(long num_reads_per_segment, double coverage_per_cell){
    wga_technique = "mda";                                                                              // Set wga_technique to MDA if this function is accessed
    double GCslope = GCBias::get_GCbias_slope();                                                        // Temporary variable to hold the GC biase slope to avoid calling the function again and again
    initial_MDAsites.clear();
    if(num_reads_per_segment == 0 || coverage_per_cell == 0.0) return false;                            // Return false if no reads can be generated from the segment
    if(coverage_per_cell>1.0){                                                                          // If coverage is greater than 1, the initial sites is obtained according to the coverage needed
        num_initial_primer_sites = num_reads_per_segment/coverage_per_cell;
    }else{                                                                                              // If coverage is =< 1, then, use the number of reads per segment as the number of initial sites.
        num_initial_primer_sites = num_reads_per_segment;
    }
    initial_MDAsites_bias.assign(num_initial_primer_sites, 0.0);                                        // Filling the vector with zeros equal to the initial number of sites
    
    for (int i=0; i<num_initial_primer_sites; i++){
        long initial_site{0};
        if(GCslope == 0.0){                                                                             // If GC bias is zero, then 
            initial_site = static_cast<long>(floor(rng::rand_double(0,1)*valid_region));                // Randomly obtain where in the valid region should the random primer be
        }else{                                                                                          // If GC bias is non-zero, the region where the random primer should be depended on the GC bias
            long initial_bin = rng::weighted_rand_int(GCbins_bias);                                     // Obtain in which bin of the valid region should the random primer be according to the GC bias
            long binStart = (initial_bin-1)*GC_binSize;                                                 // Starting location of the selected bin is the ending location of the previous bin
            long binEnd = std::min(valid_region, (initial_bin*GC_binSize));                             // Find the end locaiton of the bin. If valid region ends before the end of the full bin, then use that as the end
            initial_site = rng::rand_int(binStart,binEnd);                                              // In the selected bin, find where the initial site should be, randomly since the bias is uniform in a single bin
        }
        initial_MDAsites.push_back(initial_site);                                                       // Add this site to the vector of initial MDA sites
    }
    std::sort(initial_MDAsites.begin(), initial_MDAsites.end());                                        // Sort the vector
    long num_biased_sites = rng::rand_int(1L,static_cast<long>(ceil(initial_MDAsites.size())));         // Randomly get the number of highly biased sites among the initial sites.
    double max_bias{0.4};                                                                               // Temporary variable to hold the maximum bias allowed for each site. 
    double cum_probability{0.0};                                                                        // Temporary variable to hold the total cumilative bias (probabilty) already assigned
    for (int j=0; j<num_biased_sites; j++){
        long site_location = rng::rand_int(1L,num_initial_primer_sites)-1;                              // Pick a biased site randomly in the vector of initial primer sites
        double site_bias = rng::rand_double(0,max_bias);                                                // Randomly decide what the bias for the selected site to be
        initial_MDAsites_bias[site_location] += site_bias;                                              // Add the site bias to the respective site postion
        cum_probability += site_bias; max_bias -= site_bias;                                            // Increment cum_probability and decrement max_bias
        if (max_bias <= 0.0){                                                                           // If the max_bias reached zero before all the biased sites got a bias, break
            num_biased_sites = (j+1);                                                                   // Update the number of biased sites
            break;
        }
    }
    double remaining_probability = 1 - cum_probability;                                                 // Probability left for the unbiased sites
    double totalBias{0};                                                                                // Temporary variable to hold the total bias value. This should be 1 ideally, but can very due to round off errors
    for (int k=0; k<num_initial_primer_sites; k++){
        if (remaining_probability != 0.0){                                                              // If there is a non-zero probability remaining
            initial_MDAsites_bias[k] += remaining_probability/(num_initial_primer_sites);               // Divide the remaining probability equally among all the initial sites
        }
        totalBias += initial_MDAsites_bias[k];
    }
    if(totalBias != 0){
        for (double& bias : initial_MDAsites_bias){                                                     // Normalize the MDA bias vector to have a total probability of 1. Needed to do this because of round off errors. 
            bias/=totalBias;                                        
        }
    }
    
    if (initial_MDAsites_bias.size() == 0){                                                             // Return FALSE if there are no sites identified for read generation
        return false;
    }else{
        return true;                                                                                    // Return TRUE if atleast one site is identified to create a read from
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function generates a read with random insertions and/or deletions (indel) errors 
// based on the chrosomsome segment sequence that is being processed. If the WGA technique
// is specified to be MDA, then appropriate read distribution will be generated.
//--------------------------------------------------------------------------------------------
void ART::generate_read_with_indel(int threadID){
    long start_position;                                                                                // Variable to hold the starting position for a read
    long start_bin;
    do{
        if(wga_technique == "mda"){                                                                     // If MDA, then randomly select an amplicon and generate a read from the amplicon randomly (option only for single-cell seq)
            int amplicon = rng::weighted_rand_int(initial_MDAsites_bias,threadID)-1;                    // Choose an amplicon from the initial primer sites according to the MDA bias; -1 to get the index of that site in the initial_MDAsites vector
            long amplicon_start = initial_MDAsites[amplicon];
            long amplicon_end = std::min(valid_region, amplicon_start+12000);                           // Amplicon should have 12kb length in MDA. But if it exeeds valid region, stop at valid region.
            start_position = rng::rand_int(amplicon_start,amplicon_end,threadID);                       // Randomly select a starting position on the MDA amplicons for exponential amplification
            start_bin = (start_position/GC_binSize)+1;                                                  // The bin in which this starting point is 
        }else{                                                                                          // If not MDA, then do uniform distribution of reads
            if (GCBias::get_GCbias_slope() != 0.0){                                                     // If non-zero GC bias to be included
                long start_bin = rng::weighted_rand_int(GCbins_bias,threadID);                          // Find which bin to be selected or the read to start from
                long binStart = (start_bin-1)*GC_binSize;                                               // Starting location of the selected bin
                long binEnd = std::min(valid_region, (start_bin*GC_binSize));                           // End of the selected bin or the end of the valid region; whichever is the smallest
                start_position = rng::rand_int(binStart,binEnd,threadID);                               // Randomly obtain a start position within the selected bin
            }else{
                start_position = static_cast<long>(floor(rng::rand_double(0,1,threadID)*valid_region)); // Randomly selects a starting position on the reference sequence within a valid region 
            }
        }
    }while(forbidden_bins.find(start_bin) != forbidden_bins.end());                                     // Find if the selected bin is a forbidden bin
    
    int length_changed = get_indel_map(threadID);                                                       // Generate an indel map randomly for the read. length_changed = length of insertions - deletions      
    
    // If length_changed is negative, we need to sample more than the read length from the chromSegmentSeq to keep the read length fixed for all generated reads
    if(static_cast<int>(start_position+read_length-length_changed) > static_cast<int>((*chromSegSeq).length())){// Check if the generated read map requires a sequence that extends beyond the chromSegmentSeq size 
        length_changed = get_balanced_indel_map(threadID);                                                      // If it does, generate another indel map, which will have number of deletions <= number of insertions to avoid it happening
    } 
  
    std::string read_template_seq = (*chromSegSeq).substr(start_position, read_length-length_changed);  // Obtain the read template sequence from the chromSegmentSeq 
    int is_read_chimeric = rng::weighted_rand_int(read_artifacts_rate_dist);                            // Check if this read needs to be generated with a chimeric artifact
    if (is_read_chimeric == 1){                                                                         // If chimeric artifact need to be included, do that
        generate_chimeric_read(read_template_seq, threadID);
    }
    read_maker(read_template_seq, threadID);                                                            // Generate the read sequence by incorporating the indels into the read template sequence
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function generate the indel map, which describes where in the read will be an indel
// present. This function uses the deletion_probability_vec and the insertion_probability_vec
// to find the number of deletions and insertions a read can have at random. Then that many
// deletions and insertions are placed in the indel map at randomly obtained locations that
// corresponds to the read positions. Deletions are indicated with '-' and and insertion gets
// a random nitrogen base inserted at the location. Function takes threadID as an argument
//--------------------------------------------------------------------------------------------
int ART::get_indel_map(int threadID){
    std::map<int,char,std::less<int>>& indel_map = indel_map_vec[threadID];                             // Pass the respective map of the thread by reference
    indel_map.clear();
    int insertion_length{0};                                                                            // Variable to hold the size of the insertion made
    int deletion_length{0};                                                                             // Variable to hold the size of the deletion made
   
    // Processing deletions first
    if(deletion_rate != 0.0){
        deletion_length = static_cast<int>(rng::binomial_distribution(deletion_rate, read_length, threadID)); // Find size X
    }
    if(deletion_length>0){
        for(int j=deletion_length; j>0;){                                                               // Make X deletions at random in the read and make a map of the deletions made in a read
            int deletion_position = round((read_length-1)*rng::rand_double(0,1,threadID));              // Randomly generate a deletion location in read. Values can range [0,read_length-1]
            if(indel_map.find(deletion_position) == indel_map.end()){                                   // map.find(pos) will return where in the map we have key 'pos'. Therefore, make a deletion only if the position is not already mutated
                indel_map[deletion_position] = '-';                                                     // Indicate deletion with '-' in the indel map
                j--;
            }
        }
    }

    // Processing insertions 
    if(insertion_rate !=  0.0){
        insertion_length = static_cast<int>(rng::binomial_distribution(insertion_rate, read_length, threadID));
    }
    while((read_length-deletion_length)<insertion_length){                                              // Ensure that there is enough unchanged position for mutation after introducing deletions. Try another length if not
        insertion_length = static_cast<int>(rng::binomial_distribution(insertion_rate, read_length, threadID));
    }
    if(insertion_length>0){                                                                             // If probability of having an insertion of size X is greater than a random probability value generated, then
        for(int j=insertion_length; j>0;){
            int insertion_position = round((read_length-1)*rng::rand_double(0,1,threadID));             // Randomly generate an insertion location in read in the interval [0,read_length-1].
            if(indel_map.find(insertion_position) == indel_map.end()){                                  // Add an insertion mutation if this location does not already have a mutation
                int base = rng::rand_int(1,4,threadID);                                                 // Randomly determine the nitrogen base to add to the insertion_position
                switch(base){
                    case 1:
                        indel_map[insertion_position] = 'A';   break;
                    case 2:
                        indel_map[insertion_position] = 'C';   break;
                    case 3:
                        indel_map[insertion_position] = 'G';   break;
                    case 4:
                        indel_map[insertion_position] = 'T';  
                }
                j--;
            }
        }
    }
/* 
    // Processing deletions first
    std::vector<double> deletion_probability_vec_copy = deletion_probability_vec;
    for(int i=deletion_probability_vec.size(); i>0; i--){                                               // The size of the deletion_probability_vec is the max number of errors a read can have,if specified, else is the read length
        if(deletion_probability_vec[i-1]>=rng::rand_double(0,1,threadID)){                              // If probability of having a deletion of size X is greater than a random probability value generated, then
            deletion_length = i;                                                                        // Find size X
            for(int j=i; j>0;){                                                                         // Make X deletions at random in the read and make a map of the deletions made in a read
                int deletion_position = round((read_length-1)*rng::rand_double(0,1,threadID));          // Randomly generate a deletion location in read. Values can range [0,read_length-1]
                if(indel_map.find(deletion_position) == indel_map.end()){                               // map.find(pos) will return where in the map we have key 'pos'. Therefore, make a deletion only if the position is not already mutated
                    indel_map[deletion_position] = '-';                                                 // Indicate deletion with '-' in the indel map
                    j--;
                }
            }
            break;                                                                                      // Stop processing deletions if deletions of size X is made and mapped once
        }
    }
    // Processing insertions 
    std::vector<double> insertion_probability_vec_copy = insertion_probability_vec;
    for(int i=insertion_probability_vec.size(); i>0; i--){                                              // The size of the insertion_probability_vec is the max number of errors a read can have,if specified, else is the read length
        if((read_length-deletion_length) < i) continue;                                                 // Ensure that there is enough unchanged position for mutation after introducing deletions. Continue if not
        if(insertion_probability_vec[i-1]>=rng::rand_double(0,1,threadID)){                             // If probability of having an insertion of size X is greater than a random probability value generated, then
            insertion_length = i;
            for(int j=i; j>0;){
                int insertion_position = round((read_length-1)*rng::rand_double(0,1,threadID));         // Randomly generate an insertion location in read in the interval [0,read_length-1].
                if(indel_map.find(insertion_position) == indel_map.end()){                              // Add an insertion mutation if this location does not already have a mutation
                    int base = rng::rand_int(1,4,threadID);                                             // Randomly determine the nitrogen base to add to the insertion_position
                    switch(base){
                        case 1:
                            indel_map[insertion_position] = 'A';   break;
                        case 2:
                            indel_map[insertion_position] = 'C';   break;
                        case 3:
                            indel_map[insertion_position] = 'G';   break;
                        case 4:
                            indel_map[insertion_position] = 'T';  
                    }
                    j--;
                }
            }
            break;                                                                                      // Stop processing insertions if insertions of size X is made and mapped once
        }
    } */
    return (insertion_length-deletion_length);
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function also generate an indel map, like the function 'get_indel_map()'. However, 
// this function forces the number of deletions in a read to be less than or equal to the 
// number of insertions. Function takes threadID as an argument
//--------------------------------------------------------------------------------------------
int ART::get_balanced_indel_map(int threadID){
    std::map<int,char,std::less<int>>& indel_map = indel_map_vec[threadID];                             // Pass the respective map of the thread by reference
    indel_map.clear();
    int insertion_length{0};                                                                            // Variable to hold the size of the insertion made
    int deletion_length{0};                                                                             // Variable to hold the size of the deletion made
   
    if(insertion_rate){
        insertion_length = static_cast<int>(rng::binomial_distribution(insertion_rate, read_length, threadID));
    }
    if(insertion_length>0){
        for(int j=insertion_length; j>0;){
            int insertion_position = round((read_length-1)*rng::rand_double(0,1,threadID));             // Randomly generate an insertion location in read in the interval [0,read_length-1].
            if(indel_map.find(insertion_position) == indel_map.end()){                                  // Add an insertion mutation if this location does not already have a mutation
                int base = rng::rand_int(1,4,threadID);                                                 // Randomly determine the nitrogen base to add to the insertion_position
                switch(base){
                    case 1:
                        indel_map[insertion_position] = 'A';   break;
                    case 2:
                        indel_map[insertion_position] = 'C';   break;
                    case 3:
                        indel_map[insertion_position] = 'G';   break;
                    case 4:
                        indel_map[insertion_position] = 'T';  
                }
                j--;
            }
        }
    }


    // Processing deletions and ensures the number of deletions are not more than the number of insertions
    if(deletion_rate){
        int max_possible_deletion = std::min(insertion_length, (read_length-insertion_length));         // Limit the deletion length to the minimum of insertion length or the unchanged read portion
        deletion_length = static_cast<int>(rng::binomial_distribution(deletion_rate, max_possible_deletion, threadID)); // Make sure the maximum number of deletions is <= the number of insertions
    }                                                             
    if(deletion_length>0){                                                
        for(int j=deletion_length; j>0;){                                                               // Make X deletions at random in the read and make a map of the deletions made in a read
            int deletion_position = round((read_length-1)*rng::rand_double(0,1,threadID));              // Randomly generate a deletion location in read. End positions are invalid for deletion. Therefore floor function is used to ignore 'read_length-1'
            if(deletion_position == 0) continue;                                                        // Similarly, position '0' is also ignored
            if(indel_map.find(deletion_position) == indel_map.end()){                                   // map.find(pos) will return where in the map we have key 'pos'. Therefore, make a deletion only if the position is not already mutated
                indel_map[deletion_position] = '-';                                                     // Indicate deletion with '-' in the indel map
                j--;
            }
        }
    }
/* 
    // Processing insertions first
    std::vector<double> insertion_probability_vec_copy = insertion_probability_vec;
    for(int i=insertion_probability_vec.size(); i>0; i--){                                              // The size of the insertion_probability_vec is the max number of errors a read can have,if specified, else is the read length
        if(insertion_probability_vec[i-1]>=rng::rand_double(0,1,threadID)){                             // If probability of having an insertion of size X is greater than a random probability value generated, then
            insertion_length = i;
            for(int j=i; j>0;){
                int insertion_position = round((read_length-1)*rng::rand_double(0,1,threadID));         // Randomly generate an insertion location in read in the interval [0,read_length-1].
                if(indel_map.find(insertion_position) == indel_map.end()){                              // Add an insertion mutation if this location does not already have a mutation
                    int base = rng::rand_int(1,4,threadID);                                             // Randomly determine the nitrogen base to add to the insertion_position
                    switch(base){
                        case 1:
                            indel_map[insertion_position] = 'A';   break;
                        case 2:
                            indel_map[insertion_position] = 'C';   break;
                        case 3:
                            indel_map[insertion_position] = 'G';   break;
                        case 4:
                            indel_map[insertion_position] = 'T';  
                    }
                    j--;
                }
            }
            break;                                                                                      // Stop processing insertions if insertions of size X is made and mapped once
        }
    }
    

    // Processing deletions and ensures the number of deletions are not more than the number of insertions
    std::vector<double> deletion_probability_vec_copy = deletion_probability_vec;
    for(int i=insertion_length; i>0; i--){                                                              // Make sure the maximum number of deletions is <= the number of insertions
        if((read_length-insertion_length) < i) continue;                                                // Ensure that there is enough unchanged position for mutation after introducing insertions. Continue if not
        if(deletion_probability_vec[i-1]>=rng::rand_double(0,1,threadID)){                              // If probability of having a deletion of size X is greater than a random probability value generated, then
            deletion_length = i;                                                                        // Find size X
            for(int j=i; j>0;){                                                                         // Make X deletions at random in the read and make a map of the deletions made in a read
                int deletion_position = round((read_length-1)*rng::rand_double(0,1,threadID));          // Randomly generate a deletion location in read. End positions are invalid for deletion. Therefore floor function is used to ignore 'read_length-1'
                if(deletion_position == 0) continue;                                                    // Similarly, position '0' is also ignored
                if(indel_map.find(deletion_position) == indel_map.end()){                               // map.find(pos) will return where in the map we have key 'pos'. Therefore, make a deletion only if the position is not already mutated
                    indel_map[deletion_position] = '-';                                                 // Indicate deletion with '-' in the indel map
                    j--;
                }
            }
            break;                                                                                      // Stop processing deletions if deletions of size X is made and mapped once
        }
    } */
    return (insertion_length-deletion_length);
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function is used to generate read chimeras that occure due to errors in DNA amplification.
// This is a very crude model of chimera, which is essentially, reverse complementing a small
// segment of the read to to form a chimeric read. A segment of the proper read is randomly obtained
// which then gets replaced with it reverse complement and the modified read is returned. 
//--------------------------------------------------------------------------------------------
void ART::generate_chimeric_read(std::string& read, int threadID){
    int min_chimeric_length = static_cast<int>(read_length*1/4);                                        // Minimum is set to be atleast as long as 1/4th of the read length                         
    int max_chimeric_length = static_cast<int>(read_length*3/4);                                        // Maximum is set to be at most as long as 3/4th of the read length
    int chimeric_length = rng::rand_int(min_chimeric_length, max_chimeric_length, threadID);            // Randomly obtain the length of the chimera between the minimum and maximum lengths
    int chimeric_position = rng::rand_int(0, read_length-chimeric_length, threadID);                    // Randomly determine where in the read the chimera needs to be
    std::string read_segment = read.substr(chimeric_position,chimeric_length);                          // Cut out the portion of the read where the chimera needs to be
    std::string chimeric_insert;                                                                        // Temporary string to hold the chimeric sequence
    getReverseComplementarySeq(read_segment, chimeric_insert);                                          // Obtain the chimeric sequence as the reverse complement of the read segment
    read.replace(chimeric_position,chimeric_length,chimeric_insert);                                    // Replace the read segment with the chimeric sequence in the read to generate a chimeric read
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function converts a read template sequence obtained from a chromosome segment to a read 
// sequence by incorporating indels in specified locations as per the indel map generated. 
// It first checks if there are any indels, and if not, simply copies the templace sequence. 
// Otherwise, it loops through the templace sequence and adds characters to the read sequence 
// based on the indel map. Function takes threadID as an argument to make it thread safe
//--------------------------------------------------------------------------------------------
void ART::read_maker(std::string& read_template_seq, int threadID){
    std::string& read_seq = read_seq_vec[threadID];                                                     // Pass the respective read_seq place holder for each thread by reference
    std::map<int,char,std::less<int>>& indel_map = indel_map_vec[threadID];                             // Pass the respective indel map of the thread by reference
    read_seq.clear();
    if(indel_map.size() == 0){                                                                          // If indel map is empty, i.e, there are no insertions and deletions, then
        read_seq = read_template_seq;                                                                   // Make read sequence the same as the chromsome sequence without alterations
        return;
    }
    int k{0};
    for(size_t i=0; i<read_template_seq.size();){
        if(indel_map.count(k) == 0){                                                                    // For a base location that was unaltered, obtain the base from the chromosome sequence
            read_seq.push_back(read_template_seq[i]); i++; k++; 
        }else if(indel_map[k] == '-'){                                                                  // If the base location corresponds to a deletion, ignore that base
            i++;k++;
        }else{                                                                                          // If the base location corresponds to an insertion, get the base value from the indel map itself
            read_seq.push_back(indel_map[k]); k++;
        }
    }
    while(indel_map.count(k)>0){                                                                        // If there are more indels in the map that is not already processed, then include them as well
        read_seq.push_back(indel_map[k]);
        k++;
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function will get the read quality scores for every position at random in a read for 
// each and every read based on the read_quality_distribution vectors. Function takes threadID
// as an argument to make it thread safe
//--------------------------------------------------------------------------------------------
void ART::get_read_quality(std::vector<short>& read_quality_vec, int read_number, int threadID){
    read_quality_vec.clear();
    std::vector<std::map<unsigned int, unsigned short>>* quality_distribution_ptr = nullptr;            // Temporary vector pointer of the vector that will hold the quality distribution maps
    if(read_number == 1){                                                                               // If it is read 1, then get the quality distribution vector for read 1 
        quality_distribution_ptr = &read1_quality_distribution_vec;
    }else if(read_number == 2){                                                                         // If it is read 2, then get the quality distribution vector for read 2 
        quality_distribution_ptr = &read2_quality_distribution_vec;
    }
    std::vector<std::map<unsigned int, unsigned short>>& quality_distribution = *quality_distribution_ptr;
    
    unsigned int cumCC;                                                                                 // Temporary integer to hold the map key
    std::map<unsigned int, unsigned short>::iterator it;                                                // Declares an iterator 'it' for the map that is used to search for quality values.
    for(int i=0; i<read_length; i++){                                                                   // For each position of the read
        cumCC = static_cast<int>(ceil(rng::rand_double(0.0,1.0,threadID)*10000000))+1;                  // Randomly get a map key. Key is in the range [1,10000001]
        it = quality_distribution[i].lower_bound(cumCC);                                                // Find the iterator pointing to the first element in the map whose key is greater than or equal to cumCC
        read_quality_vec.push_back(it->second);                                                         // Get the corresponding quality score and store it in the read_quality_vec
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function adds base errors to the read after considering the quality scores in each 
// read positions, which are specified in the read_quality_vec that is passed to the funciton.
// If the base is 'N', the quality score is changed to 1 (low). For any other position, a substitution
// is randomly made based on the read quality score.
//--------------------------------------------------------------------------------------------
void ART::add_baseCall_error(std::vector<short>& read_quality_vec, int threadID){
    std::string& read_seq = read_seq_vec[threadID];                                                     // Pass the respective read_seq place holder for each thread by reference
    for(size_t i=0; i<read_quality_vec.size(); i++){                                                    // read_quality_vec has same size as the read length
        if(read_seq[i]=='N'){                                                                           // If the base in a read is 'N', then change the quality score to 1 (low)
            read_quality_vec[i]= static_cast<short>(1);
            continue;
        }
        // Sustitution of the existing base with an error
        if(rng::rand_double(0,1,threadID) < baseCall_error_probability[read_quality_vec[i]]){           // If the base call error probability is greater than a random probability, then proceed. 
            char achar = read_seq[i];                                                                   // Get the current base at that location
            while(read_seq[i] == achar){                                                                // Substitute a base that is not the original base the read had
                int base = rng::rand_int(1,4,threadID);                                                 // Randomly determine the nitrogen base to add to the base call error position
                switch(base){
                    case 1:
                        achar = 'A';   break;
                    case 2:
                        achar = 'C';   break;
                    case 3:
                        achar = 'G';   break;
                    case 4:
                        achar = 'T';  
                } 
            }
            read_seq[i] = achar;                                                                        // Make the substitution in the read sequence
        }
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function generates two paired-reads with random insertions and/or deletions (indel) errors 
// based on the chrosomsome segment sequence that is being processed. First, this function will
// randomly generate a DNA fragment based on the mean and stdDev values user-specified. Then,
// two reads will be generated from opposite ends of the same fragment.
//--------------------------------------------------------------------------------------------
bool ART::generate_paired_reads_with_indel(ART& read_1, ART& read_2, int minFragSize, std::vector<double>& fragment_weights, int threadID){
    const std::string& chromSegmentSeq = *chromSegSeq;                                                  // Variable holding chromosome segment sequence
    
    // Determine the fragment length to be simulated from the chromSegmentSequence. 
    long chromSegment_length = chromSegmentSeq.length();                                                // Variable to hold the length of the chromosome segment being processed                               
    long fragment_length{0};                                                                            // Variable to hold the fragment length that needs to be simulated
    while(true){                                                                                        // Repeat fragmentation until a successful fragment is obtained
        if(minFragSize >= chromSegment_length){                                                         // If chromSegment length <= minimumFragment length, then set fragment length to be chromSegmentSeq length
            fragment_length = chromSegment_length;
            break;
        }else{                                                                                          // Else randomly sample a fragment length from the beta distribution weights
            int minRandomFragSize = minFragSize + rng::weighted_rand_int(fragment_weights, threadID)-1;
            fragment_length = std::max(minRandomFragSize, read_length);                                 // Ensure fragment length is at least read_length
            if(fragment_length <= chromSegment_length){                                                 // Break the loop when a fragment length that is <= chromosome segment length is obtained
                break;
            }
        }
    }

    int length_changed_1 = read_1.get_indel_map(threadID);                                              // Generate an indel map randomly for read 1. length_changed = length of insertions - deletions
    int length_changed_2 = read_2.get_indel_map(threadID);                                              // Generate an indel map randomly for read 2

    // If length_changed is negative, we need to sample more than the read length from the fragment to keep the read length fixed for all generated reads
    if((read_length-length_changed_1) > fragment_length){                                               // Check if the generated read map requires a sequence that extends beyond the fragment size 
        length_changed_1 = read_1.get_balanced_indel_map(threadID);                                     // If it does, generate another indel map, which will have number of deletions <= number of insertions to avoid it happening
    } 
    if((read_length-length_changed_2) > fragment_length){                                               // Check if the generated read map appropriate for read 2 
        length_changed_2 = read_2.get_balanced_indel_map(threadID);                                     // Re-generate the indel map for read 2 if needed
    } 
    long start_position_1;                                                                              // Variable to hold the starting location of the fragment
    long start_position_2;

    if(wga_technique == "mda"){
        int amplicon = rng::weighted_rand_int(initial_MDAsites_bias,threadID)-1;                        // Choose an amplicon from the initial primer sites according to the MDA bias; -1 to get the index of that site in the initial_MDAsites vector
        long amplicon_start = initial_MDAsites[amplicon];
        long amplicon_end = std::min(chromSegment_length-fragment_length, amplicon_start+12000);        // Amplicon should have 12kb length in MDA. But if it exeeds valid region, stop at valid region.
        if(amplicon_start>amplicon_end){                                                                // Incase if the amplicon starting point is beyong the valid region for this fragment, choose a different (new) amplicon position
            start_position_1 = static_cast<long>(floor((chromSegment_length-fragment_length)*rng::rand_double(0,1, threadID)));
        }else{
            start_position_1 = rng::rand_int(amplicon_start,amplicon_end,threadID);                     // Randomly select a starting position on the MDA amplicons for exponential amplification
        }
    }else{                                                                                              // If not MDA, then do uniform distribution of reads
        if(GCBias::get_GCbias_slope() != 0.0){                                                          // If non-zero GC bias to be included
            long available_region = chromSegment_length-fragment_length;                                // Find the length of the available region in this segment for sequencing
            long num_available_fullBins = static_cast<long>(floor(available_region/GC_binSize));
            long remaining_chromLength = chromSegment_length-(num_available_fullBins*GC_binSize);       // Find the length of remaining or extra chromosome after the full bins or if the segment is smaller than one bin size
            int bins_in_remaining_chrmLength = static_cast<int>(floor(remaining_chromLength/GC_binSize));// How many full bins are remaining in the extra chrom length
            long length_of_xtraBits = remaining_chromLength-fragment_length;                            // The length of the sequence that is needed to make the available region after the available number of full GC bins
            double bias_in_xtraBits{0.0};
            if (length_of_xtraBits>=read_length){                                                       // Bias of extra bit needs to be calculated only if that segment is bigger than or equal to the read length
                if(bins_in_remaining_chrmLength==0){                                                    // The the remaining chromosome length is not big enough to be one full bin
                    bias_in_xtraBits = (GCbins_bias[num_available_fullBins]/remaining_chromLength)*length_of_xtraBits;// Caculate the bias of the xtra bits as the total bias in that non-full bin region/ the size of the region * the xtra bit size
                }else{
                    bias_in_xtraBits = (GCbins_bias[num_available_fullBins]/GC_binSize)*length_of_xtraBits;// If there are atleast one full bin in the remaining chromosome legth, then total/binSize will get a unit's bias
                }
            }
            double max_cumGCbinBias = bias_in_xtraBits;                                                 // Variable to hold the cumulative bias upto the end of the available region for the fragment to use it to normalize the new vector
            if(num_available_fullBins!=0){                                                              // If the number of full bins in available region is not zero, then add the cumulative bias up to that fina full bin with the bias of the xtra bits
                max_cumGCbinBias += cum_GCbins_bias[num_available_fullBins-1];
            }else{
                if (available_region == 0){                                                             // Special case when chromSegment_length == fragment length
                    max_cumGCbinBias = cum_GCbins_bias[0];
                }
            }
            double random_binBias = rng::rand_double(0.0,max_cumGCbinBias);                             // Randomly pick a value of cumulative bias in the available region (This approach is identical to using rng::weighted_rand_int() function)
            auto it = std::lower_bound(cum_GCbins_bias.begin(), (cum_GCbins_bias.end()), random_binBias); // Find the first position in the available region of the cum_GCsites_bias vector where the cumulative sum is greater than or equal to a randomly generated number (random_number)
            if(it != cum_GCbins_bias.end()){                                                            // If there is a bin that meets the condition
                long selected_bin = std::distance(cum_GCbins_bias.begin(), it);                         // Get the index of that bin
                long binStart = selected_bin*GC_binSize;                                                // Starting position of the selected bin    
                long binEnd = std::min(available_region, ((selected_bin+1)*GC_binSize));                // Ending position of the selected bin or the end of available region; whichever is the smallest
                start_position_1  = rng::rand_int(binStart,binEnd);                                     // Find a position in this bin to generate the read randomly, since the bias in a bin is uniform
            }else{                                                                                      // If  the condition is not met for some reason, return false. 
                return false;
            }
        }else{                                                                                          // If GC bias is zero with uniform read distribution
            start_position_1 = static_cast<long>(floor((chromSegment_length-fragment_length)*rng::rand_double(0,1, threadID))); // Starting position for read 1
        }
    }
    start_position_2 = (start_position_1+fragment_length)-(read_length-length_changed_2);               // This would be the starting position for read 2 in forward strand
    long bin_1 = (start_position_1/GC_binSize)+1; long bin_2 = (start_position_2/GC_binSize)+1;         // Find the bins in which starting position 1 and 2 belong
    if((forbidden_bins.find(bin_1) != forbidden_bins.end()) || (forbidden_bins.find(bin_2) != forbidden_bins.end())) return false; // Make sure both starting bins are not forbidden 
    //long start_position_2 = chromSegment_length-fragment_length-start_position_1;                       // This would be the starting position for read 2. Must be on the reverse-complementary strand
    // If length_changed is negative, we need to sample more than the read length from the chromSegmentSeq to keep the read length fixed for all generated reads
    /*if((start_position_1+read_length-length_changed_1) > chromSegment_length){                          // Check if the generated read map requires a sequence that extends beyond the chromSegmentSeq size 
        length_changed_1 = read_1.get_balanced_indel_map(threadID);                                     // If it does, generate another indel map, which will have number of deletions <= number of insertions to avoid it happening
    } 
    if((start_position_2+read_length-length_changed_2) > chromSegment_length){                          // Check if the generated read map appropriate for read 2 
        length_changed_2 = read_2.get_balanced_indel_map(threadID);                                     // Re-generate the indel map for read 2 if needed
    } */
    //std::string chromSegmentSeq_revComp;                                                                // Temporary string to hold the reverse complementary strand with same damages to make the paired-red from opposite end of the fragment
    //getReverseComplementarySeq(chromSegmentSeq, chromSegmentSeq_revComp);                               // Make the reverse complementary sequence strand
    std::string read1_forward_seq = chromSegmentSeq.substr(start_position_1, read_length-length_changed_1); // Obtain the read template sequence from the chromSegmentSeq for read 1 (in forward strand)
    std::string read2_forward_seq = chromSegmentSeq.substr(start_position_2, read_length-length_changed_2); // Obtain the read template sequence from the chromSegmentSeq for read 2 (in forward strand)
    std::string read1_template_seq;                                                                     // Temporary string to hold the correctly oriented read 1
    std::string read2_template_seq;                                                                     // Temporary string to hold the correctly oriented read 2                                                                     
    
    int read_orientation{1};                                                                            // Temporaray variable to hold the read pair orientation. Default is inward-orientation
    if(read_pair_orientation_dist[0]!=1.0){                                                             // If all the read pairs needs to be inward-oriented, skip the random step to save time
        read_orientation = rng::weighted_rand_int(read_pair_orientation_dist, threadID);                // Else, randomly obtain the orientation of the pair to be generated based of the ratio provided
    }
    
    int is_read1_chimeric = rng::weighted_rand_int(read_artifacts_rate_dist, threadID);                 // Check if read 1 needs to be chimeric based on the read_artifacts_rate
    int is_read2_chimeric = rng::weighted_rand_int(read_artifacts_rate_dist, threadID);                 // Check if read 2 needs to be chimeric based on the read_artifacts_rate
    
    switch (read_orientation){                                                                          // Switch condition is used instead of an if-else is to include other read orientations in the future if needed
        case 1:                                                                                         // If inward-orientated read pairs are needed (----> <----)
            read1_template_seq = read1_forward_seq;
            getReverseComplementarySeq(read2_forward_seq, read2_template_seq);                          // Make the reverse complementary sequence read
            if (is_read1_chimeric == 1){
                read_1.generate_chimeric_read(read1_template_seq, threadID);
            }
            if (is_read2_chimeric == 1){
                read_2.generate_chimeric_read(read2_template_seq, threadID);
            }
            break;
        case 2:                                                                                         // If outward-oriented read pairs are needed (<---- ---->)
            getReverseComplementarySeq(read1_forward_seq, read1_template_seq);
            read2_template_seq = read2_forward_seq;
            break;
        case 3:                                                                                         // If other orientation, generate forward-forward reads (----> ---->)
            read1_template_seq = read1_forward_seq;
            read2_template_seq = read2_forward_seq;
            break;
    }
    read_1.read_maker(read1_template_seq, threadID);                                                    // Generate the read sequence by incorporating the indels into the template sequence of read 1
    read_2.read_maker(read2_template_seq, threadID);                                                    // Generate the read sequence by incorporating the indels into the template sequence of read 2

    return true;                                                                                        // Return true if successful
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This is a function that can return the value of read_seq variable. If this function is 
// called at the end of all read operations, it will return the final read sequence that has
// indel errors and base call errors according to the quality scores. 
//--------------------------------------------------------------------------------------------
const std::string* ART::get_final_read_sequence(int threadID){
    std::string& read_seq = read_seq_vec[threadID];                                                     // Pass the respective read_seq place holder for each thread by reference
    return(&read_seq);
}