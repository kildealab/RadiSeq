#include "fastafile_handler.h"
#include "fileio.h"
#include "support_functions.h"

#include <fstream>


//--------------------------------------------------------------------------------------------
// This function will read the reference genome file provided and generate an undamaged cell
// genome template that will be further used to create damaged cell genomes. The template will
// have forward and complementary strand sequences. If the cell is diploid, then it will even 
// have two copies of each strand. They will get IDs: copy1_chr1 and copy2_chr1 repectively for 
// two copies. Forward and reverse strand sequences will be stored in seperate fasta files
// because we will need to create 2 fastq files depending on the strand. Function also returns
// the total length of the reference sequence in bp.
//--------------------------------------------------------------------------------------------
long buildUndamagedGenomeTemplate(const std::string& filepath, int nChrms, int chrmMapping, const std::string* ref_seqPath){
    std::ofstream templateFile(filepath+"/Undamaged_cell.fa");                                            // File to store forward strand sequence. If you modify the name here, change in in other function below as well
    long seqLength{0};                                                                                    // Find the total length of the genome in bp from ref_seq file    
    if (templateFile.is_open()) {
        std::ifstream refseqfile(*ref_seqPath);
        int TotalChrms = nChrms;                                                                          // Holds the total number of chromosome till the end
        std::string chromSeq_ID;                                                                          // String to hold the chromosome sequence ID. This value is not used in this function but needed to pass to the getNextChromSeq functon
        std::string chromSeq;                                                                             // String to store the forward chromseq from reference seq file
        std::string revComp_chromSeq;                                                                     // String to store the reverse complementary sequence for respective chrom-seq
        int chrmCount{0};                                                                                 // Seperate counter to set the seq ID. This value will not be same as nChrms if diploid

        switch(chrmMapping){                                                                              // Decide how to write the sequences to file depending on the chromosome mapping
            case 0:                                                                                       // If cell is haploid (chromosome mapping type 0)
                while(getNextChromSeq(refseqfile, chromSeq, chromSeq_ID) && nChrms>0){
                    uppercaseString(chromSeq);                                                            // Change lowercase -> Uppercase
                    getReverseComplementarySeq(chromSeq, revComp_chromSeq);                               // Generate the reverse complementary sequence of the chrom sequence
                    chrmCount++;
                    if(nChrms>2){                                                                         // Write according to the mapping 1 pattern
                        templateFile<<">chr"+std::to_string(chrmCount)+"a\n"<<chromSeq<<"\n";
                        templateFile<<">chr"+std::to_string(chrmCount)+"b\n"<<revComp_chromSeq<<"\n";
                        seqLength+=chromSeq.size();
                        nChrms--; 
                    }else{                                                                                // No need to have copies of X, Y chromosomes
                        templateFile<<">chrX/Y_"<<nChrms<<"a\n"<<chromSeq<<"\n";
                        templateFile<<">chrX/Y__"<<nChrms<<"b\n"<<revComp_chromSeq<<"\n";
                        seqLength+=chromSeq.size();
                        nChrms--;
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
                        templateFile<<">chr"+std::to_string(chrmCount)+"a_copy"+std::to_string(i+1)+"\n"<<chromSeq<<"\n";
                        templateFile<<">chr"+std::to_string(chrmCount)+"b_copy"+std::to_string(i+1)+"\n"<<revComp_chromSeq<<"\n";
                        seqLength+=chromSeq.size();
                        nChrms--;
                        } 
                    }else{                                                                                // No need to have copies of X, Y chromosomes
                        templateFile<<">chrX/Y_"<<nChrms<<"a\n"<<chromSeq<<"\n";
                        templateFile<<">chrX/Y__"<<nChrms<<"b\n"<<revComp_chromSeq<<"\n";
                        seqLength+=chromSeq.size();
                        nChrms--;
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
                            templateFile<<">chr"+std::to_string(chrmCount)+"a_copy"+std::to_string(i+1)+"\n"<<chromSeq<<"\n";
                            templateFile<<">chr"+std::to_string(chrmCount)+"b_copy"+std::to_string(i+1)+"\n"<<revComp_chromSeq<<"\n";
                            seqLength+=chromSeq.size();
                            nChrms--;
                        }
                        if (chrmCount>=int(TotalChrms/2) && i==0){                                        // Skip the sex chromosomes in the first for loop and write autosomes again 
                            break;
                        }else if(chrmCount>=int(TotalChrms/2) && i>0){                                    // In the second loop, write the sex chromosomes as well
                            templateFile<<">chrX/Y_"<<nChrms<<"a\n"<<chromSeq<<"\n";
                            templateFile<<">chrX/Y__"<<nChrms<<"b\n"<<revComp_chromSeq<<"\n";
                            seqLength+=chromSeq.size();
                            nChrms--;
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
        refseqfile.close();
    }
    templateFile.close();
    return seqLength;
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function will make a reverse complementary sequence for the chromSeq referenced. The
// generated sequence will be in the direction 3'-5' if the original was 5'-3'. This new sequence
// will be stored in the reference string object called revComp_chromSeq
//--------------------------------------------------------------------------------------------
void getReverseComplementarySeq(std::string& chromSeq, std::string& revComp_chromSeq){
    revComp_chromSeq.clear();                                                                             // Clear the string if there is previous seq
    revComp_chromSeq.resize(chromSeq.size());                                                             // Resize as with the forward sequence
    for(int i=0; i<chromSeq.size(); i++){
        int k=chromSeq.size()-i-1;                                                                        // Indexing is reversed to get reverse strand
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
// points will be used to break the genome into segments. 
//--------------------------------------------------------------------------------------------
void buildDamagedCellGenome(NGSsdd& SDDdata, const std::string& outputPath, const std::string& fileName){
    std::ofstream outputFile(outputPath+fileName);
    if(outputFile.is_open()){
        std::ifstream genomeTemplate(outputPath+"/Undamaged_cell.fa");                                          // Undamaged genome template file that was previously built
        std::string chromID_A;                                                                                  // Strings to hold the chrom IDs and Sequences from the file
        std::string chromSeq_A;
        std::string chromID_B;
        std::string chromSeq_B;
        long seqStartIndex{0};                                                                                  // Variable that will hold the starting index for each chromsome sequence in the genome
        std::vector<long>& baseDamageLoc1 = SDDdata.get_basestrand1_damage_loc();                               // Referencing the damage location vector to a temp vector for ease of handling
        std::vector<long>& baseDamageLoc2 = SDDdata.get_basestrand2_damage_loc();
        std::vector<long>& strandbreakLoc1 = SDDdata.get_backbone1_break_loc();
        std::vector<long>& strandbreakLoc2 = SDDdata.get_backbone2_break_loc();
        //std::vector<long>& DNAbreaks = SDDdata.get_DNA_breakPoints();

        if(genomeTemplate.is_open()){
            while(readFastaTemplate(genomeTemplate, chromID_A, chromSeq_A, chromID_B, chromSeq_B)){             // Get forward, backward sequences and their repective IDs for each chroms, one at a time
                long seq_length = chromSeq_A.size();                                                            // Sequence length is same for both A and B 
                int i{0};

                // Introudce Base damages to the forward strand sequence
                while(i<baseDamageLoc1.size()){
                    if(baseDamageLoc1[i]<=(seqStartIndex+seq_length)){                                          // If damage location is in the chromsome that is currently prcessing, then
                        chromSeq_A[baseDamageLoc1[i]-seqStartIndex-1] = 'N';                                    // Replace the nitrogen base at the damage location with 'N'. 1 is subtracted because the chrom_seq index starts from zero
                        baseDamageLoc1.erase(baseDamageLoc1.begin() + i);                                       // Remove processed locations from the original list (referenced list) of damage locations permanantly
                    }else{break;}                                                                               // Break the loop as soon as the damage location is more than the final base location of the chromosome. Works because the list is sorted. 
                }

                // Introudce Base damages to the reverse complementary strand sequence
                while(i<baseDamageLoc2.size()){
                    if(baseDamageLoc2[i]<=(seqStartIndex+seq_length)){
                        chromSeq_B[seq_length-(baseDamageLoc2[i]-seqStartIndex)] = 'N';                         // Indexing is different because the strand is not only complementary, but also revesrsed already. Index = chrom_end - (index in Genome - chrom_start)
                        baseDamageLoc2.erase(baseDamageLoc2.begin() + i);
                    }else{break;}
                }

                // Breaking chromosomes into segments, wherever there is a SSB
                
                // Processing the forward strand first 
                long A_segment_start{0};                                                                        // Temporary variale that wiil hold the starting location of each DNA segment in forward strand and update its value when a new segment is found
                long A_segment_length;
                int A_segment_index{1};                                                                         // Temporary index variable to name each DNA segment in forward strand 
                if(strandbreakLoc1.size() !=0){                                                                 // If there are strand breaks present in the forward strand, then
                    while(i<strandbreakLoc1.size()){
                        if(strandbreakLoc1[i]<=(seqStartIndex+seq_length)){                                     // If the strand break is within the chromosome that is being processed, then split sequencence into segments at each break point
                            outputFile<<chromID_A<<"_segment_"<<A_segment_index<<"\n";
                            A_segment_length = (strandbreakLoc1[i]-seqStartIndex) - A_segment_start;
                            outputFile<<chromSeq_A.substr(A_segment_start, A_segment_length)<<"\n";             // substr(a,b): get b charcters starting from index a 
                            A_segment_start = strandbreakLoc1[i]-seqStartIndex;
                            A_segment_index++;
                            strandbreakLoc1.erase(strandbreakLoc1.begin() + i);
                        }else{                                                                                  // If the strand break is not in the chromosome Or if it is the final segment in a chromosome
                            if(A_segment_index!=1){                                                             // if it is the last segment, write ID with segment tag
                                outputFile<<chromID_A<<"_segment_"<<A_segment_index<<"\n";
                                outputFile<<chromSeq_A.substr(A_segment_start)<<"\n";
                            }else{                                                                              // If the strand break is not in chromosome, write without the segmenet tag
                                outputFile<<chromID_A<<"\n";    
                                outputFile<<chromSeq_A.substr(A_segment_start)<<"\n";
                            }
                            break;
                        }
                    }    
                }else{                                                                                          // If there are no DNA breaks in the forward strand, then just copy the sequence as is without segmenting
                    outputFile<<chromID_A<<"\n"<<chromSeq_A<<"\n";
                }
                
                // Processing the reverse complementary strand next 
                long B_segment_end{seq_length};                                                                 // Temporary variable holding the end location of each segment in reverse strand
                long B_segment_length;
                int B_segment_index{1};                                                                         // Temporary index variable to name each DNA segment of reverse strand
                if(strandbreakLoc2.size()!=0){                                                                  // If there are strand breaks present in reverse strand, then
                    while(i<strandbreakLoc2.size()){                                                            
                        if(strandbreakLoc2[i]<=(seqStartIndex+seq_length)){                                     // If the strand break is within the chromosome that is being processed, then split sequencence into segments at each break point
                            outputFile<<chromID_B<<"_segment_"<<B_segment_index<<"\n";
                            B_segment_length = B_segment_end - (seq_length-(strandbreakLoc2[i]-seqStartIndex));
                            outputFile<<chromSeq_B.substr((seq_length-(strandbreakLoc2[i]-seqStartIndex)),B_segment_length)<<"\n";
                            B_segment_end = seq_length-(strandbreakLoc2[i]-seqStartIndex);                            
                            B_segment_index++;
                            strandbreakLoc2.erase(strandbreakLoc2.begin() + i);
                        }else{                                                                                  // If strand break is not in the chromosome Or if it is the final segment in a chromosome
                            if(B_segment_index!=1){                                                             // if it is the last segment, write ID with segment tag
                                outputFile<<chromID_B<<"_segment_"<<B_segment_index<<"\n";
                                outputFile<<chromSeq_B.substr(0,B_segment_end)<<"\n";
                            }else{                                                                              // If the strand break is not in chromosome, write without the segmenet tag
                                outputFile<<chromID_B<<"\n";
                                outputFile<<chromSeq_B.substr(0,B_segment_end)<<"\n";
                            }
                            break;
                        }
                    }
                }else{                                                                                          // If there are no strand breaks in the reverse strand, then just copy the sequence as is without segmenting
                    outputFile<<chromID_B<<"\n"<<chromSeq_B<<"\n";
                }

                /* // Breaking chromosomes into segments, wherever there is a DSB 
                long A_segment_start{0};                                                                        // Temporary variale that wiil hold the starting location of each DNA segment in forward strand and update its value when a new segment is found
                long A_segment_length;
                long B_segment_end{seq_length};                                                                 // Temporary variable holding the end location of each segment in reverse strand
                long B_segment_length;
                int segment_index{1};                                                                           // Temporary index variable to name each DNA segment  
                if(DNAbreaks.size()!=0){                                                                        // If there are DNA breaks present, then
                    while(i<DNAbreaks.size()){                                                            
                        if(DNAbreaks[i]<=(seqStartIndex+seq_length)){                                           // If the DNA break is within the chromosome that is being processed, then split sequencence into segments at each break point
                            outputFile<<chromID_A<<"_segment_"<<segment_index<<"\n";
                            A_segment_length = (DNAbreaks[i]-seqStartIndex) - A_segment_start;
                            outputFile<<chromSeq_A.substr(A_segment_start, A_segment_length)<<"\n";             // substr(a,b): get b charcters starting from index a 
                            A_segment_start = DNAbreaks[i]-seqStartIndex;

                            outputFile<<chromID_B<<"_segment_"<<segment_index<<"\n";
                            B_segment_length = B_segment_end - (seq_length-(DNAbreaks[i]-seqStartIndex));
                            outputFile<<chromSeq_B.substr((seq_length-(DNAbreaks[i]-seqStartIndex)),B_segment_length)<<"\n";
                            B_segment_end = seq_length-(DNAbreaks[i]-seqStartIndex);
                            
                            segment_index++;
                            DNAbreaks.erase(DNAbreaks.begin() + i);
                        }else{                                                                                  // If DNA break is not in the chromosome Or if it is the final segment in a chromosome
                            if(segment_index!=1){                                                               // if it is the last segment, write ID with segment tag
                                outputFile<<chromID_A<<"_segment_"<<segment_index<<"\n";
                                outputFile<<chromSeq_A.substr(A_segment_start)<<"\n";
                                outputFile<<chromID_B<<"_segment_"<<segment_index<<"\n";
                                outputFile<<chromSeq_B.substr(0,B_segment_end)<<"\n";
                            }else{                                                                              // If the DNA break is not in chromosome, write without the segmenet tag
                                outputFile<<chromID_A<<"\n";    
                                outputFile<<chromSeq_A.substr(A_segment_start)<<"\n";
                                outputFile<<chromID_B<<"\n";
                                outputFile<<chromSeq_B.substr(0,B_segment_end)<<"\n";
                            }
                            break;
                        }
                    }
                }else{                                                                                          // If there are no DNA breaks in the genome, then just copy the sequence as is without segmenting
                    outputFile<<chromID_A<<"\n"<<chromSeq_A<<"\n"<<chromID_B<<"\n"<<chromSeq_B<<"\n";
                } */

                seqStartIndex += seq_length;
            }
        }
        genomeTemplate.close();
    }
    outputFile.close();
}
//--------------------------------------------------------------------------------------------
