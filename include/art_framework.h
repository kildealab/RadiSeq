#ifndef ART_FRAMEWORK_H
#define ART_FRAMEWORK_H

#include <vector>
#include <map>
#include <string>
#include <unordered_set>

class ART{
    static std::vector<std::map<unsigned int, unsigned short>> read1_quality_distribution_vec;          // Static vector to hold the quality score distribution map of read 1 
    static std::vector<std::map<unsigned int, unsigned short>> read2_quality_distribution_vec;          // Static vector to hold the quality score distribution map of read 2
    static double baseCall_error_probability[80];                                                       // Static array to store the 80 values of base calling error probabilities
    static int read_length;                                                                             // Variable to hold the read_length. Static variable will be accessible by all class objects created
    static std::string* chromSegSeq;                                                                    // Variable to hold the chromosome segment sequence currently being processed
    static long valid_region;                                                                           // Variable to hold the value chromSegmentSeq.size()-read_length
    static long num_initial_primer_sites;                                                               // Variable to hold the number of initial primer sites in a chromSegment for MDA
    static std::vector<long> initial_MDAsites;                                                          // Temporary vector to hold the initial primer binding locations for MDA
    static std::vector<double> initial_MDAsites_bias;                                                   // Temporary vector to hold the bias associated with each of the initial primer binding location for MDA
    static std::string wga_technique;                                                                   // Variable to hold the value that indicate the WGA technique chosen. Default uniform amplification
    static std::unordered_set<long> forbidden_locations;                                                // A set that stores all the locations in the sequence segment where a read should not be created    
    static std::vector<double> GCsites_bias;                                                            // Vector to hold the GC bias for each location in a chromosome segment
    static std::vector<double> cum_GCsites_bias;                                                        // Vector to hold the cumilative bias upto each location in the chromosome segment

    // Following variables are private to each ART class object being generated
    double insertion_rate;
    double deletion_rate;
    std::vector<std::map<int,char,std::less<int>>> indel_map_vec;                                       // Vector that holds the maps that store the indel distribution of a read for each thread
    std::vector<std::string> read_seq_vec;                                                              // Vector to hold the final read sequence after introducing damages for each thread

public:
    
    std::vector<double> insertion_probability_vec;                                                      // Vector to hold the probability of having atleast x insertions in a read
    std::vector<double> deletion_probability_vec;                                                       // Vector to hold the probability of having atleast x deletions in a read
    
    ART();
    static void set_read_quality_distribution(const std::string&,  const std::string&);                 // Function to set the read quality distributions for read1 and read2
    static void set_baseCall_error_probability();                                                       // Function to set the base calling error probabilites

    //void set_read_error_probability(int, double, std::vector<double>&, int);                            // Function to generate the probability distribution of x indels occuring in a read
    void set_read_error_rates(int, double, double);

    bool int_set(std::string&, int);                                                                    // Function to set the chromsome segment Sequence, read lengh, indel map vector
    //const std::string* get_chromSegmentSeq();                                                           // Function to get the chrom segment sequence currently processing
    //int get_read_length();                                                                              // Function to get the read length that we want to generate
    
    static bool MDA_distribution_maker(long, double);                                                   // Function to create the coverage distribution profile for MDA for each cell
    static bool get_N_mask_and_GCbias(int, int minFragSize=1);                                          // Function to determine the forbidden locations in a chromosome fragment where reads cannot be created. minFragSize is optional

    void generate_read_with_indel(int);                                                                 // Function that that ultimately call all subfunctions to generate a read with indel damages
    int get_indel_map(int);                                                                             // Function to generate an indel distribution map for a read
    int get_balanced_indel_map(int);                                                                    // Function to generate an indel distribution map that will have deletions <= insertions
    void read_maker(std::string&, int);                                                                 // Function to make read from template using indel map

    void get_read_quality(std::vector<short>&, int, int);                                               // Function to generate a read-specific quality score profile for every read 
    void add_baseCall_error(std::vector<short>&, int);                                                  // Function to add base call errors to the read based on the quality scores of each read positions
    
    static bool generate_paired_reads_with_indel(ART&, ART&, int, std::vector<double>&, int);           // Function to generate two paired-reads from the opposite ends of the same DNA fragment
    
    const std::string* get_final_read_sequence(int);
};



#endif