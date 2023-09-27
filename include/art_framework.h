#ifndef ART_FRAMEWORK_H
#define ART_FRAMEWORK_H

#include <vector>
#include <map>
#include <string>

class ART{
    static std::vector<std::map<unsigned int, unsigned short>> read1_quality_distribution_vec;          // Static vector to hold the quality score distribution map of read 1 
    static std::vector<std::map<unsigned int, unsigned short>> read2_quality_distribution_vec;          // Static vector to hold the quality score distribution map of read 2
    static double baseCall_error_probability[80];                                                       // Static array to store the 80 values of base calling error probabilities
    static int read_length;                                                                             // Variable to hold the read_length. Static variable will be accessible by all class objects created
    static std::string& chromSegmentSeq;                                                                // Static variable to hold the chromosome segment sequence currently being processed
    static int valid_region;                                                                            // Static variable to hold the value chromSegmentSeq.size()-read_length
    static std::string read_seq;                                                                        // Static variable to hold the final read sequence after introducing damages
    std::map<int,char,std::less<int>> indel_map;                                                        // Map that stores the indel distribution of a read

public:
    
    std::vector<double> insertion_probability_vec;                                                      // Vector to hold the probability of having atleast x insertions in a read
    std::vector<double> deletion_probability_vec;                                                       // Vector to hold the probability of having atleast x deletions in a read
    
    ART();
    static void set_read_quality_distribution(const std::string&,  const std::string&);                 // Function to set the read quality distributions for read1 and read2
    static void set_baseCall_error_probability();                                                       // Function to set the base calling error probabilites

    void set_read_error_probability(int, double, std::vector<double>&, int);                            // Function to generate the probability distribution of x indels occuring in a read
    static bool init_set(int, std::string&);                                                            // Function to set all the static variables

    void generate_read_with_indel();                                                                    // Function that that ultimately call all subfunctions to generate a read with indel damages
    int get_indel_map();                                                                                // Function to generate an indel distribution map for a read
    int get_balanced_indel_map();                                                                       // Function to generate an indel distribution map that will have deletions <= insertions
    void read_maker(std::string&);                                                                      // Function to make read from template using indel map

    void get_read_quality(std::vector<short>&, int);                                                    // Function to generate a read-specific quality score profile for every read 
    void add_baseCall_error(std::vector<short>&);                                                       // Function to add base call errors to the read based on the quality scores of each read positions
    
    static void generate_paired_reads_with_indel(ART&, ART&, int, int);                                 // Function to generate two paired-reads from the opposite ends of the same DNA fragment
    
    const std::string* get_final_read_sequence();
};



#endif