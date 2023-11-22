#include "fastafile_handler.h"
#include "fileio.h"
#include "support_functions.h"

#include <fstream>
#include <sys/mman.h>


//--------------------------------------------------------------------------------------------
// This function will read the reference genome file provided and generate an undamaged cell
// genome template that will be further used to create damaged cell genomes. The template will
// have forward and complementary strand sequences. If the cell is diploid, then it will even 
// have two copies of each strand. They will get IDs: copy1_chr1 and copy2_chr1 repectively for 
// two copies. Function also returns the total length of the reference sequence in bp.
//--------------------------------------------------------------------------------------------
/* long buildUndamagedGenomeTemplate(const std::string& filepath, int nChrms, int chrmMapping, const std::string* ref_seqPath){
    std::ofstream templateFile(filepath+"/Undamaged_cell.fa");                                            // File to store forward strand sequence. If you modify the name here, change in in other function below as well
    long seqLength{0};                                                                                    // Find the total length of the genome in bp from ref_seq file    
    if (templateFile.is_open()) {
        std::ifstream refseqfile(*ref_seqPath);
        int TotalChrms = nChrms;                                                                          // Holds the total number of chromosome till the end
        std::string chromSeq_ID;                                                                          // String to hold the chromosome sequence ID. This value is not used in this function but needed to pass to the getNextChromSeq functon
        std::string chromSeq;                                                                             // String to store the forward chromseq from reference seq file
        std::string revComp_chromSeq;                                                                     // String to store the reverse complementary sequence for respective chrom-seq
        int chrmCount{0};                                                                                 // Seperate counter to set the seq ID. This value will not be same as nChrms if diploid
        const int batchSize = 4000;                                                                       // Define a batch size for writing data to the output file. These much data will be stored in cache before writing it on the file
        std::vector<std::string> batch_buffer;                                                            // Create a buffer for storing output data.
            
        switch(chrmMapping){                                                                              // Decide how to write the sequences to file depending on the chromosome mapping
            case 0:                                                                                       // If cell is haploid (chromosome mapping type 0)
                while(getNextChromSeq(refseqfile, chromSeq, chromSeq_ID) && nChrms>0){
                    uppercaseString(chromSeq);                                                            // Change lowercase -> Uppercase
                    getReverseComplementarySeq(chromSeq, revComp_chromSeq);                               // Generate the reverse complementary sequence of the chrom sequence
                    chrmCount++;
                    if(nChrms>2){                                                                         // Write according to the mapping 1 pattern
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"a\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"b\n"+revComp_chromSeq+"\n");
                        seqLength+=chromSeq.size();
                        nChrms--; 
                    }else{                                                                                // No need to have copies of X, Y chromosomes
                        batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"a\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"b\n"+revComp_chromSeq+"\n");
                        seqLength+=chromSeq.size();
                        nChrms--;
                    }
                    if (batch_buffer.size() >= batchSize) {                                               // Check if the batch buffer is full, and write it to the file if needed.
                        writeBatchToFile(batch_buffer, templateFile);
                    }
                }
                break;
            case 1:                                                                                       // If the cell is diploid with chromosome mapping type 1 (1,1,2,2,....22,22,X,Y)
                while(getNextChromSeq(refseqfile, chromSeq, chromSeq_ID) && nChrms>0){
                    uppercaseString(chromSeq);                                                            // Change lowercase -> Uppercase
                    getReverseComplementarySeq(chromSeq, revComp_chromSeq);                               // Generate the reverse complementary sequence of the chrom sequence
                    chrmCount++;
                    if(nChrms>2){                                                                         // Until all autosomes are done, 
                        for(int i=0; i<2; i++){                                                           // Write according to the mapping 1 pattern
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"a_copy"+std::to_string(i+1)+"\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"b_copy"+std::to_string(i+1)+"\n"+revComp_chromSeq+"\n");
                        seqLength+=chromSeq.size();
                        nChrms--;
                        } 
                    }else{                                                                                // No need to have copies of X, Y chromosomes
                        batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"a\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"b\n"+revComp_chromSeq+"\n");
                        seqLength+=chromSeq.size();
                        nChrms--;
                    }
                    if (batch_buffer.size() >= batchSize) {                                               // Check if the batch buffer is full, and write it to the file if needed.
                        writeBatchToFile(batch_buffer, templateFile);
                    }
                }
                break;
            case 2:                                                                                       // If the cell is diploid with chromosome mapping type 2 (1,2.....22,1,2.....22,X,Y)
                for(int i=0; i<2; i++){
                    refseqfile.clear();                                                                   // Clear any error flags
                    refseqfile.seekg(0, std::ios::beg);                                                   // Set the position to the beginning of the file
                    chrmCount = 0; 
                    while(getNextChromSeq(refseqfile, chromSeq, chromSeq_ID) && nChrms>0){
                        uppercaseString(chromSeq);                                                        // Change lowercase -> Uppercase
                        getReverseComplementarySeq(chromSeq, revComp_chromSeq);                           // Generate the reverse complementary sequence of the chrom sequence
                        chrmCount++;
                        if(chrmCount<int(TotalChrms/2)){                                                  // Write the autosomes once in the for loop
                            batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"a_copy"+std::to_string(i+1)+"\n"+chromSeq+"\n");
                            batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"b_copy"+std::to_string(i+1)+"\n"+revComp_chromSeq+"\n");
                            seqLength+=chromSeq.size();
                            nChrms--;
                        }
                        if (chrmCount>=int(TotalChrms/2) && i==0){                                        // Skip the sex chromosomes in the first for loop and write autosomes again 
                            break;
                        }else if(chrmCount>=int(TotalChrms/2) && i>0){                                    // In the second loop, write the sex chromosomes as well
                            batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"a\n"+chromSeq+"\n");
                            batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"b\n"+revComp_chromSeq+"\n");
                            seqLength+=chromSeq.size();
                            nChrms--;
                        }
                        if (batch_buffer.size() >= batchSize) {                                           // Check if the batch buffer is full, and write it to the file if needed.
                            writeBatchToFile(batch_buffer, templateFile);
                        }
                    }
                } 
                break;
        }
        if(nChrms!=0){                                                                                    // If fewer sequences than expected was found in the reference file, then error and exit
            std::cerr<<"\n ERROR: Only "<<chrmCount<<" chromosomse seqences were found in "<<*ref_seqPath<<"\n";
            std::cerr<<" A sequence file with "<<TotalChrms<<" chromosomes arranged in 1,2,3.....X,Y fashion is expected \n";
            exit(EXIT_FAILURE);
        }
        writeBatchToFile(batch_buffer, templateFile);                                                     // If there are unwritten data in batch buffer, write that too when the loop ends
        refseqfile.close();
    }
    templateFile.close();
    return seqLength;
} */
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function will process the reference genome memory-map data provided and generate an 
// undamaged cell genome template that will be further used to create damaged cell genomes. 
// The template will have forward and complementary strand sequences. If the cell is diploid, 
// then it will even have two copies of each strand. They will get IDs: copy1_chr1 and copy2_chr1
// repectively for two copies. Function also returns the total length of the reference sequence
// in bp. Function will also populate a memory map of the file for easy handling.
//--------------------------------------------------------------------------------------------
long buildUndamagedGenomeTemplate_MM(char* templateFileMapping, std::size_t templateFileSize, int nChrms, int chrmMapping, const std::string* ref_seqPath){
    char* position_in_MM = templateFileMapping;                                                           // Pointer to the current position in the MM as we write. Starts with the pointer to the beginning
    long seqLength{0};                                                                                    // Find the total length of the genome in bp from ref_seq file  
    if (position_in_MM == nullptr){return 0;}                                                             // Return if the template file memory map is not valid
    
    int TotalChrms = nChrms;                                                                              // Holds the total number of chromosome till the end
    std::string chromSeq_ID;                                                                              // String to hold the chromosome sequence ID. This value is not used in this function but needed to pass to the getNextChromSeq functon
    std::string chromSeq;                                                                                 // String to store the forward chromseq from reference seq file
    std::string revComp_chromSeq;                                                                         // String to store the reverse complementary sequence for respective chrom-seq
    int chrmCount{0};                                                                                     // Seperate counter to set the seq ID. This value will not be same as nChrms if diploid
    
    size_t position{0};                                                                                   // Temporary variable to hold the last-read position in the memory map of the reference file
    size_t refFileSize;                                                                                   // A variable to hold the file size of the reference genome, during memory-mapping
    void* refFileMM = generateInputFileMemoryMap(*ref_seqPath, refFileSize);                              // Create the memory-map of the reference genome file
    const char* refSeqData = static_cast<char*>(refFileMM);                                               // Casting the memory-map void pointer to a const char pointer for further processing

    const int batchSize{10};                                                                              // Define a batch size for writing data to the output memory-mapped file. These much data (buffer vector elements) will be stored in cache before writing it on the file
    std::vector<std::string> batch_buffer;                                                                // Create a buffer for storing output data.

    switch(chrmMapping){                                                                                  // Decide how to write the sequences to file depending on the chromosome mapping
        case 0:                                                                                           // If cell is haploid (chromosome mapping type 0)
            while(getNextChromSeq_MM(refSeqData, refFileSize, position, chromSeq, chromSeq_ID) && nChrms>0){
                uppercaseString(chromSeq);                                                                // Change lowercase -> Uppercase
                getReverseComplementarySeq(chromSeq, revComp_chromSeq);                                   // Generate the reverse complementary sequence of the chrom sequence
                chrmCount++;
                if(nChrms>2){                                                                             // Write according to the mapping 1 pattern
                    batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"a\n"+chromSeq+"\n");
                    batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"b\n"+revComp_chromSeq+"\n");
                    seqLength+=chromSeq.size();                                                           // Increment the reference seq length
                    nChrms--; 
                }else{                                                                                    // No need to have copies of X, Y chromosomes
                    batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"a\n"+chromSeq+"\n");
                    batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"b\n"+revComp_chromSeq+"\n");
                    seqLength+=chromSeq.size();
                    nChrms--;
                }
                if (batch_buffer.size() >= batchSize) {                                                   // Check if the batch buffer is full, and write it to the memory-mapped file if needed.
                    writeBatchToMMFile(batch_buffer, position_in_MM, templateFileMapping, templateFileSize);
                }    
            }
            break;
        case 1:                                                                                           // If the cell is diploid with chromosome mapping type 1 (1,1,2,2,....22,22,X,Y)
            while(getNextChromSeq_MM(refSeqData, refFileSize, position, chromSeq, chromSeq_ID) && nChrms>0){
                uppercaseString(chromSeq);                                                                // Change lowercase -> Uppercase
                getReverseComplementarySeq(chromSeq, revComp_chromSeq);                                   // Generate the reverse complementary sequence of the chrom sequence
                chrmCount++;
                if(nChrms>2){                                                                             // Until all autosomes are done, 
                    for(int i=0; i<2; i++){                                                               // Write according to the mapping 1 pattern
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"a_copy"+std::to_string(i+1)+"\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"b_copy"+std::to_string(i+1)+"\n"+revComp_chromSeq+"\n");
                        seqLength+=chromSeq.size();
                        nChrms--;
                    } 
                }else{                                                                                    // No need to have copies of X, Y chromosomes
                    batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"a\n"+chromSeq+"\n");
                    batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"b\n"+revComp_chromSeq+"\n");
                    seqLength+=chromSeq.size();
                    nChrms--;
                }
                if (batch_buffer.size() >= batchSize) {                                                   // Check if the batch buffer is full, and write it to the memory-mapped file if needed.
                    writeBatchToMMFile(batch_buffer, position_in_MM, templateFileMapping, templateFileSize);
                }  
            }
            break;
        case 2:                                                                                           // If the cell is diploid with chromosome mapping type 2 (1,2.....22,1,2.....22,X,Y)
            for(int i=0; i<2; i++){
                position = 0;                                                                             // Reset the starting positon of the memory map in each iteration
                chrmCount = 0; 
                while(getNextChromSeq_MM(refSeqData, refFileSize, position, chromSeq, chromSeq_ID) && nChrms>0){
                    uppercaseString(chromSeq);                                                            // Change lowercase -> Uppercase
                    getReverseComplementarySeq(chromSeq, revComp_chromSeq);                               // Generate the reverse complementary sequence of the chrom sequence
                    chrmCount++;
                    if(chrmCount<int(TotalChrms/2)){                                                      // Write the autosomes once in the for loop
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"a_copy"+std::to_string(i+1)+"\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chr"+std::to_string(chrmCount)+"b_copy"+std::to_string(i+1)+"\n"+revComp_chromSeq+"\n");
                        seqLength+=chromSeq.size();
                        nChrms--;
                    }
                    if (chrmCount>=int((TotalChrms/2)-1) && i==0){                                        // Skip the sex chromosomes in the first for loop and write autosomes again 
                        break;
                    }else if(chrmCount>=int(TotalChrms/2) && i>0){                                        // In the second loop, write the sex chromosomes as well
                        batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"a\n"+chromSeq+"\n");
                        batch_buffer.push_back(">chrX/Y_"+std::to_string(nChrms)+"b\n"+revComp_chromSeq+"\n");
                        seqLength+=chromSeq.size();
                        nChrms--;
                    }
                    if (batch_buffer.size() >= batchSize) {                                               // Check if the batch buffer is full, and write it to the memory-mapped file if needed.
                        writeBatchToMMFile(batch_buffer, position_in_MM, templateFileMapping, templateFileSize);
                    }
                }
            } 
            break;
    }
    if(nChrms!=0){                                                                                        // If fewer sequences than expected was found in the reference file, then error and exit
        std::cerr<<"\n ERROR: Only "<<chrmCount<<" chromosomse seqences were found in the reference file\n";
        std::cerr<<" A sequence file with "<<TotalChrms<<" chromosomes arranged in 1,2,3.....X,Y fashion is expected \n";
        exit(EXIT_FAILURE);
    } 
    munmap(refFileMM, refFileSize);                                                                       // Unmap the reference genome file to avoid memory-leaks
    writeBatchToMMFile(batch_buffer, position_in_MM, templateFileMapping, templateFileSize);              // If there are unwritten data in batch buffer, write that too when the loop ends                 
    return seqLength;
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function will make a reverse complementary sequence for the chromSeq referenced. The
// generated sequence will be in the direction 3'-5' if the original was 5'-3'. This new sequence
// will be stored in the reference string object called revComp_chromSeq
//--------------------------------------------------------------------------------------------
void getReverseComplementarySeq(const std::string& chromSeq, std::string& revComp_chromSeq){
    revComp_chromSeq.clear();                                                                             // Clear the string if there is previous seq
    revComp_chromSeq.resize(chromSeq.size());                                                             // Resize as with the forward sequence
    for(size_t i=0; i<chromSeq.size(); i++){
        size_t k=chromSeq.size()-i-1;                                                                        // Indexing is reversed to get reverse strand
        switch(chromSeq[i]){
            case 'A':                                                                                     // Generate complementary base values accordingly
                revComp_chromSeq[k]='T'; break;
            case 'C':
                revComp_chromSeq[k]='G'; break;
            case 'G':
                revComp_chromSeq[k]='C'; break;
            case 'T':
                revComp_chromSeq[k]='A'; break;
            default:                                                                                      // Any unrecognized nitrogen base including 'N' will be set to 'N'
                revComp_chromSeq[k]='N';
        }
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function will take the strand breaks, base damages and DNA breakpoint information obtained
// from SDD files and generate a damaged cell genome and store it into a fasta file. For every cell
// a fasta file will be made. The undamaged genome template previously generated will be used to 
// create the damaged cell genome. Base damage locations will be changed to an 'N' and DNA break
// points will be used to break the genome into segments. This function will return an integer
// value that is the total number of lines written in that particular fasta file for the damaged cell
//--------------------------------------------------------------------------------------------
/* int buildDamagedCellGenome(NGSsdd& SDDdata, const std::string& outputPath, const std::string& fileName) {
    std::ofstream outputFile(outputPath + fileName);
    if (!outputFile.is_open()) {return 0;}                                                                      // Return if the output file couldn't be opened

    std::ifstream genomeTemplate(outputPath + "/Undamaged_cell.fa");                                            // Undamaged genome template file that was previously built
    if (!genomeTemplate.is_open()) {return 0;}                                                                  // Return if the genome template file couldn't be opened

    int num_of_lines_written{0};                                                                                // This is a variable to hold the total number of lines written to each damage cell file
    std::string chromID_A;                                                                                      // Strings to hold the chrom IDs and Sequences from the file
    std::string chromSeq_A;
    std::string chromID_B;
    std::string chromSeq_B;
    long seqStartIndex = 0;                                                                                     // Variable that will hold the starting index for each chromsome sequence in the genome
    std::vector<long>& baseDamageLoc1 = SDDdata.get_basestrand1_damage_loc();                                   // Referencing the damage location vector to a temp vector for ease of handling
    std::vector<long>& baseDamageLoc2 = SDDdata.get_basestrand2_damage_loc();
    std::vector<long>& strandbreakLoc1 = SDDdata.get_backbone1_break_loc();
    std::vector<long>& strandbreakLoc2 = SDDdata.get_backbone2_break_loc();
    int i{0}; int j{0}; int k{0}; int l{0};                                                                     // Counter variables to keep a tab on the elements in the damage location vector that we already processed
    const int batchSize = 4000;                                                                                 // Define a batch size for writing data to the output file. These much data will be stored in cache before writing it on the file
    std::vector<std::string> batch_buffer;                                                                      // Create a buffer for storing output data.
        
    while (readFastaTemplate(genomeTemplate, chromID_A, chromSeq_A, chromID_B, chromSeq_B)) {                   // Get forward, backward sequences and their repective IDs for each chroms, one at a time
        long seq_length = chromSeq_A.size();                                                                        // Sequence length is same for both A and B 

        // Introudce Base damages to the forward strand sequence
        while (i<baseDamageLoc1.size()) {                  
            if (baseDamageLoc1[i]<=(seqStartIndex + seq_length)) {                                              // If damage location is in the chromsome that is currently prcessing, then
                chromSeq_A[baseDamageLoc1[i] - seqStartIndex - 1] = 'N';                                        // Replace the nitrogen base at the damage location with 'N'. 1 is subtracted because the chrom_seq index starts from zero
                i++;                                                                                            // Increment the counter to know where the processing should start next in the damage location vector
            }else {break;}                                                                                      // Break the loop as soon as the damage location is more than the final base location of the chromosome. Works because the list is sorted. 
        }

        // Introudce Base damages to the reverse complementary strand sequence
        while (j<baseDamageLoc2.size()) {
            if (baseDamageLoc2[j]<=(seqStartIndex + seq_length)) {
                chromSeq_B[seq_length - (baseDamageLoc2[j] - seqStartIndex)] = 'N';                             // Indexing is different because the strand is not only complementary, but also revesrsed already. Index = chrom_end - (index in Genome - chrom_start)
                j++;
            } else {break;}
        }

        // Breaking chromosomes into segments, wherever there is a SSB
                
        // Processing the forward strand first 
        long A_segment_start{0};                                                                                // Temporary variale that wiil hold the starting location of each DNA segment in forward strand and update its value when a new segment is found
        long A_segment_length;
        int A_segment_index{1};                                                                                 // Temporary index variable to name each DNA segment in forward strand 
        if(strandbreakLoc1.size() !=0){                                                                         // If there are strand breaks present in the forward strand, then
            while(k<strandbreakLoc1.size()){
                if(strandbreakLoc1[k]<=(seqStartIndex+seq_length)){                                             // If the strand break is within the chromosome that is being processed, then split sequencence into segments at each break point
                    batch_buffer.push_back(chromID_A+"_segment_"+std::to_string(A_segment_index)+"\n");
                    A_segment_length = (strandbreakLoc1[k]-seqStartIndex) - A_segment_start;
                    batch_buffer.push_back(chromSeq_A.substr(A_segment_start, A_segment_length)+"\n");                     // substr(a,b): get b charcters starting from index a 
                    A_segment_start = strandbreakLoc1[k]-seqStartIndex;
                    A_segment_index++; k++;
                    num_of_lines_written += 2;
                }else{                                                                                          // If the strand break is not in the chromosome Or if it is the final segment in a chromosome
                    if(A_segment_index!=1){                                                                     // if it is the last segment, write ID with segment tag
                        batch_buffer.push_back(chromID_A+"_segment_"+std::to_string(A_segment_index)+"\n");
                        batch_buffer.push_back(chromSeq_A.substr(A_segment_start)+"\n");
                        num_of_lines_written += 2;
                    }else{                                                                                      // If the strand break is not in chromosome, write without the segmenet tag
                        batch_buffer.push_back(chromID_A+"\n");    
                        batch_buffer.push_back(chromSeq_A.substr(A_segment_start)+"\n");
                        num_of_lines_written += 2;
                    }
                    break;
                }
            }    
        }else{                                                                                                  // If there are no DNA breaks in the forward strand, then just copy the sequence as is without segmenting
            batch_buffer.push_back(chromID_A+"\n"+chromSeq_A+"\n");
            num_of_lines_written += 2;
        }
        
        // Processing the reverse complementary strand next 
        long B_segment_end{seq_length};                                                                         // Temporary variable holding the end location of each segment in reverse strand
        long B_segment_length;
        int B_segment_index{1};                                                                                 // Temporary index variable to name each DNA segment of reverse strand
        if(strandbreakLoc2.size()!=0){                                                                          // If there are strand breaks present in reverse strand, then
            while(l<strandbreakLoc2.size()){                                                            
                if(strandbreakLoc2[l]<=(seqStartIndex+seq_length)){                                             // If the strand break is within the chromosome that is being processed, then split sequencence into segments at each break point
                    batch_buffer.push_back(chromID_B+"_segment_"+std::to_string(B_segment_index)+"\n");
                    B_segment_length = B_segment_end - (seq_length-(strandbreakLoc2[l]-seqStartIndex));
                    batch_buffer.push_back(chromSeq_B.substr((seq_length-(strandbreakLoc2[l]-seqStartIndex)),B_segment_length)+"\n");
                    B_segment_end = seq_length-(strandbreakLoc2[l]-seqStartIndex);                            
                    B_segment_index++; l++;
                    num_of_lines_written += 2;
                }else{                                                                                          // If strand break is not in the chromosome Or if it is the final segment in a chromosome
                    if(B_segment_index!=1){                                                                     // if it is the last segment, write ID with segment tag
                        batch_buffer.push_back(chromID_B+"_segment_"+std::to_string(B_segment_index)+"\n");
                        batch_buffer.push_back(chromSeq_B.substr(0,B_segment_end)+"\n");
                        num_of_lines_written += 2;
                    }else{                                                                                      // If the strand break is not in chromosome, write without the segmenet tag
                        batch_buffer.push_back(chromID_B+"\n");
                        batch_buffer.push_back(chromSeq_B.substr(0,B_segment_end)+"\n");
                        num_of_lines_written += 2;
                    }
                    break;
                }
            }
        }else{                                                                                                  // If there are no strand breaks in the reverse strand, then just copy the sequence as is without segmenting
            batch_buffer.push_back(chromID_B+"\n"+chromSeq_B+"\n");
            num_of_lines_written += 2;
        }

        if (batch_buffer.size() >= batchSize) {                                                                 // Check if the batch buffer is full, and write it to the file if needed.
            writeBatchToFile(batch_buffer, outputFile);
        }
*/
        /* // Breaking chromosomes into segments, wherever there is a DSB 
        long A_segment_start{0};                                                                                // Temporary variale that wiil hold the starting location of each DNA segment in forward strand and update its value when a new segment is found
        long A_segment_length;
        long B_segment_end{seq_length};                                                                         // Temporary variable holding the end location of each segment in reverse strand
        long B_segment_length;
        int segment_index{1};                                                                                   // Temporary index variable to name each DNA segment  
        if(DNAbreaks.size()!=0){                                                                                // If there are DNA breaks present, then
            while(i<DNAbreaks.size()){                                                            
                if(DNAbreaks[i]<=(seqStartIndex+seq_length)){                                                   // If the DNA break is within the chromosome that is being processed, then split sequencence into segments at each break point
                    outputFile<<chromID_A<<"_segment_"<<segment_index<<"\n";
                    A_segment_length = (DNAbreaks[i]-seqStartIndex) - A_segment_start;
                    outputFile<<chromSeq_A.substr(A_segment_start, A_segment_length)<<"\n";                     // substr(a,b): get b charcters starting from index a 
                    A_segment_start = DNAbreaks[i]-seqStartIndex;

                    outputFile<<chromID_B<<"_segment_"<<segment_index<<"\n";
                    B_segment_length = B_segment_end - (seq_length-(DNAbreaks[i]-seqStartIndex));
                    outputFile<<chromSeq_B.substr((seq_length-(DNAbreaks[i]-seqStartIndex)),B_segment_length)<<"\n";
                    B_segment_end = seq_length-(DNAbreaks[i]-seqStartIndex);
                    
                    segment_index++;
                    DNAbreaks.erase(DNAbreaks.begin() + i);
                }else{                                                                                          // If DNA break is not in the chromosome Or if it is the final segment in a chromosome
                    if(segment_index!=1){                                                                       // if it is the last segment, write ID with segment tag
                        outputFile<<chromID_A<<"_segment_"<<segment_index<<"\n";
                        outputFile<<chromSeq_A.substr(A_segment_start)<<"\n";
                        outputFile<<chromID_B<<"_segment_"<<segment_index<<"\n";
                        outputFile<<chromSeq_B.substr(0,B_segment_end)<<"\n";
                    }else{                                                                                      // If the DNA break is not in chromosome, write without the segmenet tag
                        outputFile<<chromID_A<<"\n";    
                        outputFile<<chromSeq_A.substr(A_segment_start)<<"\n";
                        outputFile<<chromID_B<<"\n";
                        outputFile<<chromSeq_B.substr(0,B_segment_end)<<"\n";
                    }
                    break;
                }
            }
        }else{                                                                                                  // If there are no DNA breaks in the genome, then just copy the sequence as is without segmenting
            outputFile<<chromID_A<<"\n"<<chromSeq_A<<"\n"<<chromID_B<<"\n"<<chromSeq_B<<"\n";
        } */
/*
        seqStartIndex += seq_length;
    }

    writeBatchToFile(batch_buffer, outputFile);                                                                 // If there are unwritten data in batch buffer, write that too when the loop ends

    genomeTemplate.close();
    outputFile.close();
    return (num_of_lines_written);
} */
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function is almost identical to buildDamagedCellGenome function. The only difference is
// that this function uses the memory map of the undamaged genome template file to build the
// damaged genomes instead of using the file itself. If this function is called multiple times,
// this approach reduces the I/O operations on a large file hence improving performance. 
//--------------------------------------------------------------------------------------------
int buildDamagedCellGenome_from_MM(NGSsdd& SDDdata, const std::string& outputPath, const std::string& fileName, char* genomeTemplate_data, size_t templateSize, int groupTID){
    size_t outFileSize = templateSize+10000;                                                                    // Temporary variable to hold the size of the current file being processed. This get dynamically changed if needed. Starts with undamagedfile size +10kb extra
    std::string outFilePath = outputPath+fileName;                                                              // Path to the output fasta file that we want to write
    char* outFileMapping = createMemoryMappedFile(outFilePath,outFileSize);                                     // Create a memory mapped output file for each damaged cell genome
    char* position_in_MM = outFileMapping;                                                                      // Pointer to the current position in the MM as we write. Starts with the pointer to the beginning
    if (position_in_MM == nullptr){return 0;}                                                                   // Return if the template file memory map is not valid

    int num_of_lines_written{0};                                                                                // This is a variable to hold the total number of lines written to each damage cell file
    std::string chromID_A;                                                                                      // Strings to hold the chrom IDs and Sequences from the file
    std::string chromSeq_A;
    std::string chromID_B;
    std::string chromSeq_B;
    long seqStartIndex = 0;                                                                                     // Variable that will hold the starting index for each chromsome sequence in the genome
    std::vector<long>& baseDamageLoc1 = SDDdata.get_basestrand1_damage_loc(groupTID);                           // Referencing the damage location vector to a temp vector for ease of handling
    std::vector<long>& baseDamageLoc2 = SDDdata.get_basestrand2_damage_loc(groupTID);
    std::vector<long>& strandbreakLoc1 = SDDdata.get_backbone1_break_loc(groupTID);
    std::vector<long>& strandbreakLoc2 = SDDdata.get_backbone2_break_loc(groupTID);
    size_t i{0}; size_t j{0}; size_t k{0}; size_t l{0};                                                         // Counter variables to keep a tab on the elements in the damage location vector that we already processed
    const int batchSize = 100;                                                                                  // Define a batch size for writing data to the output file. These much data will be stored in cache before writing it on the file
    std::vector<std::string> batch_buffer;                                                                      // Create a buffer for storing output data.
    size_t position = 0;                                                                                        // Temporary variable to hold the last read position in the memory map
    
    while (readFastaMemoryMap(genomeTemplate_data, templateSize, position, chromID_A, chromSeq_A, chromID_B, chromSeq_B)){// Get forward, backward sequences and their repective IDs for each chroms, one at a time
        long seq_length = chromSeq_A.size();                                                                    // Sequence length is same for both A and B 
        
        // Introudce Base damages to the forward strand sequence
        while (i<baseDamageLoc1.size()) {                  
            if (baseDamageLoc1[i]<=(seqStartIndex + seq_length)) {                                              // If damage location is in the chromsome that is currently prcessing, then
                chromSeq_A[baseDamageLoc1[i] - seqStartIndex - 1] = 'N';                                        // Replace the nitrogen base at the damage location with 'N'. 1 is subtracted because the chrom_seq index starts from zero
                i++;                                                                                            // Increment the counter to know where the processing should start next in the damage location vector
            }else {break;}                                                                                      // Break the loop as soon as the damage location is more than the final base location of the chromosome. Works because the list is sorted. 
        }

        // Introudce Base damages to the reverse complementary strand sequence
        while (j<baseDamageLoc2.size()) {
            if (baseDamageLoc2[j]<=(seqStartIndex + seq_length)) {
                chromSeq_B[seq_length - (baseDamageLoc2[j] - seqStartIndex)] = 'N';                             // Indexing is different because the strand is not only complementary, but also revesrsed already. Index = chrom_end - (index in Genome - chrom_start)
                j++;
            } else {break;}
        }

        // Breaking chromosomes into segments, wherever there is a SSB
                
        // Processing the forward strand first 
        long A_segment_start{0};                                                                                // Temporary variale that wiil hold the starting location of each DNA segment in forward strand and update its value when a new segment is found
        long A_segment_length;
        int A_segment_index{1};                                                                                 // Temporary index variable to name each DNA segment in forward strand 
        if(k<strandbreakLoc1.size()){                                                                           // If there are unprocessed strand breaks present in the forward strand, then
            while(k<strandbreakLoc1.size()){
                if(strandbreakLoc1[k]<=(seqStartIndex+seq_length)){                                             // If the strand break is within the chromosome that is being processed, then split sequencence into segments at each break point
                    batch_buffer.push_back(chromID_A+"_segment_"+std::to_string(A_segment_index)+"\n");
                    A_segment_length = (strandbreakLoc1[k]-seqStartIndex) - A_segment_start;
                    batch_buffer.push_back(chromSeq_A.substr(A_segment_start, A_segment_length)+"\n");          // substr(a,b): get b charcters starting from index a 
                    A_segment_start = strandbreakLoc1[k]-seqStartIndex;
                    A_segment_index++; k++;
                    num_of_lines_written += 2;
                }else{                                                                                          // If the strand break is not in the chromosome Or if it is the final segment in a chromosome
                    if(A_segment_index!=1){                                                                     // if it is the last segment, write ID with segment tag
                        batch_buffer.push_back(chromID_A+"_segment_"+std::to_string(A_segment_index)+"\n");
                        batch_buffer.push_back(chromSeq_A.substr(A_segment_start)+"\n");
                        num_of_lines_written += 2;
                    }else{                                                                                      // If the strand break is not in chromosome, write without the segmenet tag
                        batch_buffer.push_back(chromID_A+"\n");    
                        batch_buffer.push_back(chromSeq_A.substr(A_segment_start)+"\n");
                        num_of_lines_written += 2;
                    }
                    break;
                }
            }    
        }else{                                                                                                  // If there are no DNA breaks in the forward strand, then just copy the sequence as is without segmenting
            batch_buffer.push_back(chromID_A+"\n"+chromSeq_A+"\n");
            num_of_lines_written += 2;
        }
        
        // Processing the reverse complementary strand next 
        long B_segment_end{seq_length};                                                                         // Temporary variable holding the end location of each segment in reverse strand
        long B_segment_length;
        int B_segment_index{1};                                                                                 // Temporary index variable to name each DNA segment of reverse strand
        if(l<strandbreakLoc2.size()){                                                                           // If there are unprocessed strand breaks present in reverse strand, then
            while(l<strandbreakLoc2.size()){                                                            
                if(strandbreakLoc2[l]<=(seqStartIndex+seq_length)){                                             // If the strand break is within the chromosome that is being processed, then split sequencence into segments at each break point
                    batch_buffer.push_back(chromID_B+"_segment_"+std::to_string(B_segment_index)+"\n");
                    B_segment_length = B_segment_end - (seq_length-(strandbreakLoc2[l]-seqStartIndex));
                    batch_buffer.push_back(chromSeq_B.substr((seq_length-(strandbreakLoc2[l]-seqStartIndex)),B_segment_length)+"\n");
                    B_segment_end = seq_length-(strandbreakLoc2[l]-seqStartIndex);                            
                    B_segment_index++; l++;
                    num_of_lines_written += 2;
                }else{                                                                                          // If strand break is not in the chromosome Or if it is the final segment in a chromosome
                    if(B_segment_index!=1){                                                                     // if it is the last segment, write ID with segment tag
                        batch_buffer.push_back(chromID_B+"_segment_"+std::to_string(B_segment_index)+"\n");
                        batch_buffer.push_back(chromSeq_B.substr(0,B_segment_end)+"\n");
                        num_of_lines_written += 2;
                    }else{                                                                                      // If the strand break is not in chromosome, write without the segmenet tag
                        batch_buffer.push_back(chromID_B+"\n");
                        batch_buffer.push_back(chromSeq_B.substr(0,B_segment_end)+"\n");
                        num_of_lines_written += 2;
                    }
                    break;
                }
            }
        }else{                                                                                                  // If there are no strand breaks in the reverse strand, then just copy the sequence as is without segmenting
            batch_buffer.push_back(chromID_B+"\n"+chromSeq_B+"\n");
            num_of_lines_written += 2;
        }

        if (batch_buffer.size() >= batchSize) {                                                                 // Check if the batch buffer is full, and write it to the file if needed.
 //           #pragma omp critical
            {
                writeBatchToMMFile(batch_buffer, position_in_MM, outFileMapping, outFileSize, outFilePath);
            }
        }
        /* // Breaking chromosomes into segments, wherever there is a DSB 
        long A_segment_start{0};                                                                                // Temporary variale that wiil hold the starting location of each DNA segment in forward strand and update its value when a new segment is found
        long A_segment_length;
        long B_segment_end{seq_length};                                                                         // Temporary variable holding the end location of each segment in reverse strand
        long B_segment_length;
        int segment_index{1};                                                                                   // Temporary index variable to name each DNA segment  
        if(DNAbreaks.size()!=0){                                                                                // If there are DNA breaks present, then
            while(i<DNAbreaks.size()){                                                            
                if(DNAbreaks[i]<=(seqStartIndex+seq_length)){                                                   // If the DNA break is within the chromosome that is being processed, then split sequencence into segments at each break point
                    outputFile<<chromID_A<<"_segment_"<<segment_index<<"\n";
                    A_segment_length = (DNAbreaks[i]-seqStartIndex) - A_segment_start;
                    outputFile<<chromSeq_A.substr(A_segment_start, A_segment_length)<<"\n";                     // substr(a,b): get b charcters starting from index a 
                    A_segment_start = DNAbreaks[i]-seqStartIndex;

                    outputFile<<chromID_B<<"_segment_"<<segment_index<<"\n";
                    B_segment_length = B_segment_end - (seq_length-(DNAbreaks[i]-seqStartIndex));
                    outputFile<<chromSeq_B.substr((seq_length-(DNAbreaks[i]-seqStartIndex)),B_segment_length)<<"\n";
                    B_segment_end = seq_length-(DNAbreaks[i]-seqStartIndex);
                    
                    segment_index++;
                    DNAbreaks.erase(DNAbreaks.begin() + i);
                }else{                                                                                          // If DNA break is not in the chromosome Or if it is the final segment in a chromosome
                    if(segment_index!=1){                                                                       // if it is the last segment, write ID with segment tag
                        outputFile<<chromID_A<<"_segment_"<<segment_index<<"\n";
                        outputFile<<chromSeq_A.substr(A_segment_start)<<"\n";
                        outputFile<<chromID_B<<"_segment_"<<segment_index<<"\n";
                        outputFile<<chromSeq_B.substr(0,B_segment_end)<<"\n";
                    }else{                                                                                      // If the DNA break is not in chromosome, write without the segmenet tag
                        outputFile<<chromID_A<<"\n";    
                        outputFile<<chromSeq_A.substr(A_segment_start)<<"\n";
                        outputFile<<chromID_B<<"\n";
                        outputFile<<chromSeq_B.substr(0,B_segment_end)<<"\n";
                    }
                    break;
                }
            }
        }else{                                                                                                  // If there are no DNA breaks in the genome, then just copy the sequence as is without segmenting
            outputFile<<chromID_A<<"\n"<<chromSeq_A<<"\n"<<chromID_B<<"\n"<<chromSeq_B<<"\n";
        } */

        seqStartIndex += seq_length;
    }
//    #pragma omp critical
    {
        writeBatchToMMFile(batch_buffer, position_in_MM, outFileMapping, outFileSize, outFilePath);                 // If there are unwritten data in batch buffer, write that too when the loop ends
    }
    if (msync(outFileMapping, outFileSize, MS_SYNC) == -1) {                                                    // After writing your data to the memory-mapped file using mmap, before closing the file, call msync to flush/sync the changes.
        perror("\nERROR: Failed to synchronize memory-mapped data to file.\n");
    }
    munmap(outFileMapping, outFileSize);                                                                        // Unmap the memory-map to avoid memory leaks after use
    return (num_of_lines_written);
}
//--------------------------------------------------------------------------------------------