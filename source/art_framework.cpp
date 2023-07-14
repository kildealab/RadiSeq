#include "art_framework.h"
#include "random_generator.h"
#include "fileio.h"
#include "fastafile_handler.h"

#include <vector>
#include <iostream>
#include <map>
#include <fstream>



//--------------------------------------------------------------------------------------------
// Definition of all static variables that can be accessed by all ART class objects
//--------------------------------------------------------------------------------------------
std::vector<std::map<unsigned int, unsigned short>> ART::read1_quality_distribution_vec;                // Quality distribution vector for read 1. This vector will not change for each cell 
std::vector<std::map<unsigned int, unsigned short>> ART::read2_quality_distribution_vec;                // Quality distribution vector for read 2. This vector will not change for each cell 
double ART::baseCall_error_probability[80];                                                             // Array that will hold 80 values of base call error probabilities
int ART::read_length{0};                                                                                // Read length 
std::string tmp ="";
std::string& ART::chromSegmentSeq{tmp};                                                                 // Sequence of the chromosome segment being processed
int ART::valid_region{0};                                                                               // valid_region = chromSegmentSeq.size()-read_length
std::string ART::read_seq;                                                                              // This will be the read sequence from read template with added indels
//--------------------------------------------------------------------------------------------



// Default constructor
ART::ART(){}
    


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
void ART::set_read_error_probability(int read_length, double error_rate, std::vector<double>& error_probability_vec, int max_errors){
    error_probability_vec.clear();                                                                      // Empty the vector 
    double cdf_cutoff = 0.999999;                                                                       // Threshold value for total probability
    
    if(max_errors == 0) return;                                                                         // If max errors is set to 0 (no errors), then return

    if(error_rate<0.000000000000000000000000000001) return;                                             // If user-specified error rate is < 10^-30, set it 0
    
    // Probability that there are no errors in a read of length read_length for a sequencing error rate given
    double probability = rng::binomial_cdf(0, error_rate, read_length);
    double total_probability = probability;
    
    // Determine the probability that there are atleast i errors in a read of length read_length for the sequencing error_rate
    for(size_t i=1; i<read_length; i++){
        probability = rng::binomial_cdf(i, error_rate, read_length);
        error_probability_vec.push_back(probability);
        if(max_errors>0 && (i>=max_errors)) break;                                                      // Return when user-specified limit of errors reached
        total_probability += probability; 
        if(total_probability >= cdf_cutoff) break;                                                      // If total probability exceeded the threshold, return
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function sets all the static members of the ART class
//--------------------------------------------------------------------------------------------
bool ART::init_set(int read_len, std::string& chromSegSeq){
    read_length = read_len;
    chromSegmentSeq = chromSegSeq;
    valid_region = chromSegmentSeq.size() - read_length;
    if(valid_region<0) return false;                                                                    // If read length is > chromSegmentSeq length, then do not proceed with the chromSegmentSeq
    return true;
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function generates a read with random insertions and/or deletions (indel) errors 
// based on the chrosomsome segment sequence that is being processed.
//--------------------------------------------------------------------------------------------
void ART::generate_read_with_indel(){
    long start_position = static_cast<long>(floor(rng::rand_double(0,1)*valid_region));                 // Randomly selects a starting position on the reference sequence within a valid region 
    int length_changed = get_indel_map();                                                               // Generate an indel map randomly for the read. length_changed = length of insertions - deletions      
    
    // If length_changed is negative, we need to sample more than the read length from the chromSegmentSeq to keep the read length fixed for all generated reads
    if((start_position+read_length-length_changed) > chromSegmentSeq.length()){                         // Check if the generated read map requires a sequence that extends beyond the chromSegmentSeq size 
        length_changed = get_balanced_indel_map();                                                      // If it does, generate another indel map, which will have number of deletions <= number of insertions to avoid it happening
    } 
  
    std::string read_template_seq = chromSegmentSeq.substr(start_position, read_length-length_changed); // Obtain the read template sequence from the chromSegmentSeq 
    
    read_maker(read_template_seq);                                                                      // Generate the read sequence by incorporating the indels into the read template sequence
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function generate the indel map, which describes where in the read will be an indel
// present. This function uses the deletion_probability_vec and the insertion_probability_vec
// to find the number of deletions and insertions a read can have at random. Then that many
// deletions and insertions are placed in the indel map at randomly obtained locations that
// corresponds to the read positions. Deletions are indicated with '-' and and insertion gets
// a random nitrogen base inserted at the location.
//--------------------------------------------------------------------------------------------
int ART::get_indel_map(){
    indel_map.clear();
    int insertion_length{0};                                                                            // Variable to hold the size of the insertion made
    int deletion_length{0};                                                                             // Variable to hold the size of the deletion made
    
    // Processing deletions first
    for(int i=deletion_probability_vec.size(); i>0; i--){                                               // The size of the deletion_probability_vec is the max number of errors a read can have,if specified, else is the read length
        if(deletion_probability_vec[i-1]>=rng::rand_double(0,1)){                                       // If probability of having a deletion of size X is greater than a random probability value generated, then
            deletion_length = i;                                                                        // Find size X
            for(int j=i; j>0;){                                                                         // Make X deletions at random in the read and make a map of the deletions made in a read
                int deletion_position = floor((read_length-1)*rng::rand_double(0,1));                   // Randomly generate a deletion location in read. End positions are invalid for deletion. Therefore floor function is used to ignore 'read_length-1'
                if(deletion_position == 0) continue;                                                    // Similarly, position '0' is also ignored
                if(indel_map.count(deletion_position) == 0){                                            // map.count(pos) will return number of elements in the map with the key 'pos'. Therefore, make a deletion only if the position is not already mutated
                    indel_map[deletion_position] = '-';                                                 // Indicate deletion with '-' in the indel map
                    j--;
                }
            }
            break;                                                                                      // Stop processing deletions if deletions of size X is made and mapped once
        }
    }
    // Processing insertions 
    for(int i=insertion_probability_vec.size(); i>0; i--){                                              // The size of the insertion_probability_vec is the max number of errors a read can have,if specified, else is the read length
        if((read_length-deletion_length) < i) continue;                                                 // Ensure that there is enough unchanged position for mutation after introducing deletions. Continue if not
        if(insertion_probability_vec[i-1]>=rng::rand_double(0,1)){                                      // If probability of having an insertion of size X is greater than a random probability value generated, then
            insertion_length = i;
            for(int j=i; j>0;){
                int insertion_position = round((read_length-1)*rng::rand_double(0,1));                  // Randomly generate an insertion location in read in the interval [0,read_length-1].
                if(indel_map.count(insertion_position) == 0){                                           // Add an insertion mutation if this location does not already have a mutation
                    int base = rng::rand_int(1,4);                                                      // Randomly determine the nitrogen base to add to the insertion_position
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
    return (insertion_length-deletion_length);
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function also generate an indel map, like the function 'get_indel_map()'. However, 
// this function forces the number of deletions in a read to be less than or equal to the 
// number of insertions. 
//--------------------------------------------------------------------------------------------
int ART::get_balanced_indel_map(){
    indel_map.clear();
    int insertion_length{0};                                                                            // Variable to hold the size of the insertion made
    int deletion_length{0};                                                                             // Variable to hold the size of the deletion made
    

    // Processing insertions first
    for(int i=insertion_probability_vec.size(); i>0; i--){                                              // The size of the insertion_probability_vec is the max number of errors a read can have,if specified, else is the read length
        if(insertion_probability_vec[i-1]>=rng::rand_double(0,1)){                                      // If probability of having an insertion of size X is greater than a random probability value generated, then
            insertion_length = i;
            for(int j=i; j>0;){
                int insertion_position = round((read_length-1)*rng::rand_double(0,1));                  // Randomly generate an insertion location in read in the interval [0,read_length-1].
                if(indel_map.count(insertion_position) == 0){                                           // Add an insertion mutation if this location does not already have a mutation
                    int base = rng::rand_int(1,4);                                                      // Randomly determine the nitrogen base to add to the insertion_position
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
    for(int i=insertion_length; i>0; i--){                                                              // Make sure the maximum number of deletions is <= the number of insertions
        if((read_length-insertion_length) < i) continue;                                                // Ensure that there is enough unchanged position for mutation after introducing insertions. Continue if not
        if(deletion_probability_vec[i-1]>=rng::rand_double(0,1)){                                       // If probability of having a deletion of size X is greater than a random probability value generated, then
            deletion_length = i;                                                                        // Find size X
            for(int j=i; j>0;){                                                                         // Make X deletions at random in the read and make a map of the deletions made in a read
                int deletion_position = floor((read_length-1)*rng::rand_double(0,1));                   // Randomly generate a deletion location in read. End positions are invalid for deletion. Therefore floor function is used to ignore 'read_length-1'
                if(deletion_position == 0) continue;                                                    // Similarly, position '0' is also ignored
                if(indel_map.count(deletion_position) == 0){                                            // map.count(pos) will return number of elements in the map with the key 'pos'. Therefore, make a deletion only if the position is already mutated
                    indel_map[deletion_position] = '-';                                                 // Indicate deletion with '-' in the indel map
                    j--;
                }
            }
            break;                                                                                      // Stop processing deletions if deletions of size X is made and mapped once
        }
    }
    return (insertion_length-deletion_length);
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function converts a read template sequence obtained from a chromosome segment to a read 
// sequence by incorporating indels in specified locations as per the indel map generated. 
// It first checks if there are any indels, and if not, simply copies the templace sequence. 
// Otherwise, it loops through the templace sequence and adds characters to the read sequence 
// based on the indel map.
//--------------------------------------------------------------------------------------------
void ART::read_maker(std::string& read_template_seq){
    read_seq.clear();
    if(indel_map.size() == 0){                                                                          // If indel map is empty, i.e, there are no insertions and deletions, then
        read_seq = read_template_seq;                                                                   // Make read sequence the same as the chromsome sequence without alterations
        return;
    }
    int k{0};
    for(int i=0; i<read_template_seq.size();){
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
// each and every read based on the read_quality_distribution vectors
//--------------------------------------------------------------------------------------------
void ART::get_read_quality(std::vector<short>& read_quality_vec, int read_number){
    std::vector<std::map<unsigned int, unsigned short>>* quality_distribution_ptr = nullptr;            // Temporary vector pointer of the vector that will hold the quality distribution maps
    if(read_number == 1){                                                                               // If it is read 1, then get the quality distribution vector for read 1 
        quality_distribution_ptr = &read1_quality_distribution_vec;
    }else if(read_number == 2){                                                                         // If it is read 2, then get the quality distribution vector for read 2 
        quality_distribution_ptr = &read2_quality_distribution_vec;
    }
    std::vector<std::map<unsigned int, unsigned short>>& quality_distribution = *quality_distribution_ptr;
    
    unsigned int cumCC;                                                                                 // Temporary integer to hold the map key
    std::map<unsigned int, unsigned short>::iterator it;
    for(int i=0; i<read_length; i++){                                                                   // For each position of the read
        cumCC = static_cast<int>(ceil(rng::rand_double(0,1)*1000000));                                  // Randomly get a map key
        it = quality_distribution[i].lower_bound(cumCC);
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
void ART::add_baseCall_error(std::vector<short>& read_quality_vec){
    for(int i=0; i<read_quality_vec.size(); i++){                                                       // read_quality_vec has same size as the read length
        if(read_seq[i]=='N'){                                                                           // If the base in a read is 'N', then change the quality score to 1 (low)
            read_quality_vec[i]= static_cast<short>(1);
            continue;
        }
        // Sustitution of the existing base with an error
        if(rng::rand_double(0,1) < baseCall_error_probability[read_quality_vec[i]]){                    // If the base call error probability is greater than a random probability, then proceed. 
            char achar = read_seq[i];                                                                   // Get the current base at that location
            while(read_seq[i] == achar){                                                                // Substitute a base that is not the original base the read had
                int base = rng::rand_int(1,4);                                                          // Randomly determine the nitrogen base to add to the base call error position
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
void ART::generate_paired_reads_with_indel(ART& read_1, ART& read_2, int mean_frag_length, int std_dev_frag_length){
    // Determine the fragment length to be simulated from the chromSegmentSequence. 
    int fragment_length;                                                                                // Variable to hold the fragment length that needs to be simulated
    if(mean_frag_length - (2*std_dev_frag_length) > chromSegmentSeq.length()){                          // If chromSegment length < mean - 2*stdDev, then set fragment length to be chromSegmentSeq length
	   fragment_length = chromSegmentSeq.length();
    }else{                                                                                              // Else randomly sample a fragment length from the normal distribution of specific mean and std_dev user-specified
	    fragment_length = static_cast<int>(round(rng::normal_distribution(mean_frag_length, std_dev_frag_length)));
        while (fragment_length<read_length || fragment_length>chromSegmentSeq.length()){                // Sample again and again until the generated fragment length satisfies the required conditions
		    fragment_length = static_cast<int>(round(rng::normal_distribution(mean_frag_length, std_dev_frag_length)));
	    }   
    }

    long start_position_1 = static_cast<long>(floor((chromSegmentSeq.length()-fragment_length)*rng::rand_double(0,1))); // Starting position for read 1
    long start_position_2 = chromSegmentSeq.length()-fragment_length-start_position_1;                  // This would be the starting position for read 2. Must be on the reverse-complementary strand

    int length_changed_1 = read_1.get_indel_map();                                                      // Generate an indel map randomly for read 1. length_changed = length of insertions - deletions
    int length_changed_2 = read_2.get_indel_map();                                                      // Generate an indel map randomly for read 2

    // If length_changed is negative, we need to sample more than the read length from the chromSegmentSeq to keep the read length fixed for all generated reads
    if((start_position_1+read_length-length_changed_1) > chromSegmentSeq.length()){                     // Check if the generated read map requires a sequence that extends beyond the chromSegmentSeq size 
        length_changed_1 = read_1.get_balanced_indel_map();                                             // If it does, generate another indel map, which will have number of deletions <= number of insertions to avoid it happening
    } 
    if((start_position_2+read_length-length_changed_2) > chromSegmentSeq.length()){                     // Check if the generated read map appropriate for read 2 
        length_changed_2 = read_2.get_balanced_indel_map();                                             // Re-generate the indel map for read 2 if needed
    } 
    
    std::string chromSegmentSeq_revComp;                                                                // Temporary string to hold the reverse complementary strand with same damages to make the paired-red from opposite end of the fragment
    getReverseComplementarySeq(chromSegmentSeq, chromSegmentSeq_revComp);                               // Make the reverse complementary sequence strand

    std::string read1_template_seq = chromSegmentSeq.substr(start_position_1, read_length-length_changed_1);         // Obtain the read template sequence from the chromSegmentSeq for read 1
    std::string read2_template_seq = chromSegmentSeq_revComp.substr(start_position_2, read_length-length_changed_2); // Obtain the read template sequence from the chromSegmentSeq for read 2
    
    read_1.read_maker(read1_template_seq);                                                              // Generate the read sequence by incorporating the indels into the template sequence of read 1
    read_2.read_maker(read2_template_seq);                                                              // Generate the read sequence by incorporating the indels into the template sequence of read 2
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This is a function that can return the value of read_seq variable. If this function is 
// called at the end of all read operations, it will return the final read sequence that has
// indel errors and base call errors according to the quality scores. 
//--------------------------------------------------------------------------------------------
const std::string* ART::get_final_read_sequence(){
    return(&read_seq);
}