#include "sddfile_handler.h"
#include "fileio.h"
#include "support_functions.h"
#include "random_generator.h"

#include <algorithm>
#include <cmath>
#include <omp.h>

// Default constructor
NGSsdd::NGSsdd(){
}



//--------------------------------------------------------------------------------------------
// This function processes all the SDD files given by the user
//--------------------------------------------------------------------------------------------
void NGSsdd::process_sddfile(NGSsdd& SDDdata, NGSParameters& parameter){

    // function returns a pointer to a vector holding SDD file paths. Therefore, de-reference the pointer to get the actual vector object
    const std::vector<std::string>& sdd_paths = parameter.get_sddfile_path();                          
    set_num_of_SDDs(sdd_paths.size());                                                                  // number of SDD files passed
    std::vector<int> exposure_count_vector;                                                             // vector to hold number of exposures from all SDDs
    exposureLine_vec.resize(get_num_of_SDDs());                                                         // vector to hold vectors of line numbers that corresponds to new exposures. It will store one vector for every SDD file
    
    for (int i=0; i<get_num_of_SDDs(); i++){
        bool readFullHeader = (i == 0);                                                                 // True for the first iteration only
        readSDDfileHeader(&sdd_paths[i], SDDdata, readFullHeader);                                      // full header will be read only for the first SDD file, otherwise only the dose is read
        exposure_count_vector.push_back(countExposuresSDD(&sdd_paths[i], exposureLine_vec[i]));         // count exposure in each SDD one by one and also store the line numbers of each new exposure
    }
    set_num_of_exposures(*std::min_element(exposure_count_vector.begin(), exposure_count_vector.end()));// find the least (common) number of exposures in all SDD(s) to merge 
    // If there are unequal number of exposure data entries (runs) in SDDs that there is to combine, print warning
    if(!(std::adjacent_find(exposure_count_vector.begin(), exposure_count_vector.end(), std::not_equal_to<int>()) == exposure_count_vector.end())){
        std::cerr<<"\n WARNING: The number of exposures in multiple SDD files you have provided to merge damages are not the same.\n"
        <<" RadiSeq will only merge exposures that are present in all SDD files and ignore the rest of the data.\n"
        <<" ---- Found "<<get_num_of_exposures()<<" exposure data to merge from SDDs ----- \n";
    }
    if(parameter.get_num_of_cells_to_sequence()<get_num_of_exposures()){                                // Check if the number of cells we want to sequence is less than the number of cells we have in the SDD file
        set_num_of_damagedCells_toBuild(parameter.get_num_of_cells_to_sequence());                      // if true, build only the required number of damaged cell genomes
    }else{
        set_num_of_damagedCells_toBuild(get_num_of_exposures());                                        // otherwise build all the damaged cell genomes in the SDD file
    }
    if(parameter.get_adjust_damages_with_actual_dose()){                                                // if adjust damages with actual dose flag is ON, read actual dose values
        // function returns a pointer to a vector holding actual dose file paths. Therefore, de-reference the pointer to get the actual vector object
        const std::vector<std::string>& actualdose_paths = *parameter.get_actual_dosefile_path();       
        for(int i=0; i<get_num_of_SDDs(); i++){
            actual_dose_delivered.push_back(readActualDosefile(&actualdose_paths[i], get_num_of_exposures()));
        }
    }
    lines_read_sdd.resize(get_num_of_SDDs());                                                           // one value each for every SDD 
    for (size_t i = 0; i < lines_read_sdd.size(); i++){                                                 // initiating the vector with empty values
        lines_read_sdd[i] = std::streampos(std::streamoff(0));                                          // streamoff(0) is setting the position to be at the beginning of the stream
    }

}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// set_ functions to set SDD data variable values
//--------------------------------------------------------------------------------------------
void NGSsdd::set_num_of_SDDs(int sddValue){
    num_of_SDDs = sddValue;
}
void NGSsdd::set_num_of_exposures(int sddValue){
    num_of_exposures = sddValue;
}
void NGSsdd::set_num_of_damagedCells_toBuild(int sddValue){
    num_of_damagedCells_toBuild = sddValue;
}
void NGSsdd::set_expected_dose_gy(std::string* sddField){
    std::vector<std::string> temp_vec;                                                                  // Temporary vector to hold stringToVec function return
    stringToVec(',', sddField, temp_vec);                                                               // Converts comma-seperated sddField to temp_vec
    expected_dose_gy.push_back(std::stold(temp_vec[1]));                                                // Obtain the dose value and convert string to double
}
void NGSsdd::set_chrom_size_bp(std::string* sddField){
    std::vector<std::string> temp_vec;
    stringToVec(',', sddField, temp_vec);
    for(int i=0; i<std::stoi(temp_vec[1]); i++){                                                        // Iterate through each chromosome sizes. temp_vec[1] is the number of chroms
        if (temp_vec[i+2].back() == ';'){                                                               // Check if the ';' character is at the end. Last entry might have it.
            temp_vec[i+2].erase(temp_vec[i+2].size() - 1);                                              // Remove the ';' character
        }
        chrom_size_bp.push_back(std::stold(temp_vec[i+2])*1e6);                                         // Convert chrom size from Mbp to bp before storing
    }
}
void NGSsdd::set_cell_ploidy_and_chrom_mappping(std::string* sddField){
    std::vector<std::string> temp_vec;
    stringToVec(',', sddField, temp_vec);
    num_chroms = std::stoi(temp_vec[1]);
    if (temp_vec[2]==temp_vec[3]){                                                                      // It is diploid (same len chrom copy is there) and also the chrom mapping is 1,1,2,2......22,22,X,Y
        chrom_mapping_type = 1;
    }else if (temp_vec[2]==temp_vec[(num_chroms/2)+1]){                                                 // It is diploid and the chrom mapping is 1,2,......22,1,2,......22,X,Y
        chrom_mapping_type = 2;
    }else{                                                                                              // Haploid
        chrom_mapping_type = 0;
    }
    set_chrom_end_loc();
}
void NGSsdd::set_chrom_end_loc(){
    long chrom_end = 0;
    for(size_t i=0; i<chrom_size_bp.size(); i++){
        chrom_end += chrom_size_bp[i];
        chrom_end_loc.push_back(chrom_end);                                                             // Store chrom ends as a vector for later use. Includes 0 and end of genome as well.
    }
}
void NGSsdd::set_lines_read_sdd(std::streampos sddPos, int sddfilenumber){
    lines_read_sdd[sddfilenumber] = sddPos;                                                             // Update the position in SDD file
}
void NGSsdd::init_set_data_holders(int nGroupThreads){
    original_num_damages.clear();
    original_num_damages.resize(nGroupThreads);
    
    damage_start_loc_vec.clear();
    damage_start_loc_vec.resize(nGroupThreads);
    
    chrom_ID_vec.clear();
    chrom_ID_vec.resize(nGroupThreads);
    
    temp_backbone1_breaks_vec.clear();
    temp_backbone1_breaks_vec.resize(nGroupThreads);
    
    temp_backbone2_breaks_vec.clear();
    temp_backbone2_breaks_vec.resize(nGroupThreads);
    
    temp_basestrand1_damages_vec.clear();
    temp_basestrand1_damages_vec.resize(nGroupThreads);
    
    temp_basestrand2_damages_vec.clear();
    temp_basestrand2_damages_vec.resize(nGroupThreads);
    
    backbone1_break_loc_vec.clear();
    backbone1_break_loc_vec.resize(nGroupThreads);
    
    backbone2_break_loc_vec.clear();
    backbone2_break_loc_vec.resize(nGroupThreads);
    
    basestrand1_damage_loc_vec.clear();
    basestrand1_damage_loc_vec.resize(nGroupThreads);
    
    basestrand2_damage_loc_vec.clear();
    basestrand2_damage_loc_vec.resize(nGroupThreads);
}
void NGSsdd::init_set_workThread_data_holders(int groupTID, int nWorkerThreads){
    damage_start_loc_vec[groupTID].resize(nWorkerThreads);
    chrom_ID_vec[groupTID].resize(nWorkerThreads);
    temp_backbone1_breaks_vec[groupTID].resize(nWorkerThreads);
    temp_backbone2_breaks_vec[groupTID].resize(nWorkerThreads);
    temp_basestrand1_damages_vec[groupTID].resize(nWorkerThreads);
    temp_basestrand2_damages_vec[groupTID].resize(nWorkerThreads);
}
void NGSsdd::set_all_sdd_data_fields(std::string* sddEntry, int groupTID, int workerTID){               // This function process sdd damage entry line as a whole and passes each field to respective set functions
    std::vector<std::string> temp_vec;
    stringToVec(';', sddEntry, temp_vec);
    set_sddData_field3(&temp_vec[2], groupTID, workerTID);
    set_sddData_field4(&temp_vec[3], groupTID, workerTID);
    set_sddData_field7(&temp_vec[6], groupTID, workerTID);
}
void NGSsdd::set_sddData_field3(std::string* sddField3, int groupTID, int workerTID){
    std::vector<std::string> temp_vec;
    stringToVec(',', sddField3, temp_vec);                                                              // temp_vec[1] is chromosome ID, temp_vec[2] is chromatid ID. 
    chrom_ID_vec[groupTID][workerTID] = std::stoi(temp_vec[1]);                                         // Store chromosome ID. Chrom ID = 1 for first chromosome. Assuming there is only one chromatid per chromosome for now.  
}
void NGSsdd::set_sddData_field4(std::string* sddField4, int groupTID, int workerTID){
    damage_start_loc_vec[groupTID][workerTID] = std::stol(*sddField4);
}
void NGSsdd::set_sddData_field7(std::string* sddField7, int groupTID, int workerTID){
    long damage_location{0};                                                                            // Temporary variable will hold the damage location for each damage block in each chromosome
    int chrom_ID = chrom_ID_vec[groupTID][workerTID];
    long damage_start_loc = damage_start_loc_vec[groupTID][workerTID];
    std::vector<std::string> temp_vec1;                                                                 // This temporary vector will hold all damages in a field7. 1 element = 1 damage
    stringToVec('/', sddField7, temp_vec1);                                                             // Isolate individual damage descriptions 
    for(size_t i=0; i<temp_vec1.size(); i++){                                                           // Iterate through each damage in a damage block
        std::vector<int> temp_vec2;                                                                     // This temporary vector will hold one damage at a time. 3 elements corresponding to 3 comma seperated damage description
        stringToVec(',', &temp_vec1[i], temp_vec2);                                                     
        if (temp_vec2[2] == 0){continue;}                                                               // Skip the damage if it is a non-damage interaction
        else{
            #pragma omp atomic                                                                          // Ensure that the following counter will be updated by only one thread at a time
            original_num_damages[groupTID]++;                                                           // Count total number of damages in one exposure in one SDD if it is a valid damage
            switch(temp_vec2[0]){                                                                       // temp_vec2[0] indicate where the damage is: base vs backbone
                case 1:                                                                                 // Backbone 1
                    damage_location = chrom_end_loc[chrom_ID-1]+damage_start_loc+(temp_vec2[1]-1);      // damage location = end_loc_of previous chrom+ start_loc_of damage block wrt chrom + (damage location w.r.t damage block - 1). 
                    temp_backbone1_breaks_vec[groupTID][workerTID].push_back(damage_location);          // 1 is subtracted because damage location w.r.t damage block starts from 1 not 0.
                    break;
                case 2:                                                                                 // Base strand 1
                    
                    damage_location = chrom_end_loc[chrom_ID-1]+damage_start_loc+(temp_vec2[1]-1);
                    temp_basestrand1_damages_vec[groupTID][workerTID].push_back(damage_location);
                    break;
                case 3:                                                                                 // Base strand 2
                    damage_location = chrom_end_loc[chrom_ID-1]+damage_start_loc+(temp_vec2[1]-1);
                    temp_basestrand2_damages_vec[groupTID][workerTID].push_back(damage_location);
                    break;
                case 4:                                                                                 // Backbone 2
                    damage_location = chrom_end_loc[chrom_ID-1]+damage_start_loc+(temp_vec2[1]-1);
                    temp_backbone2_breaks_vec[groupTID][workerTID].push_back(damage_location);
                    break;
            }
        }
    }
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// get_ functions to return the SDD data variable values
//--------------------------------------------------------------------------------------------
int NGSsdd::get_num_of_SDDs(){
    return(num_of_SDDs);
} 
int NGSsdd::get_num_of_exposures(){
    return(num_of_exposures);
} 
int NGSsdd::get_num_of_damagedCells_toBuild(){
    return(num_of_damagedCells_toBuild);
}
const std::vector<std::vector<double>>* NGSsdd::get_actual_dose_delivered(){
    return(&actual_dose_delivered);
}
double NGSsdd::get_expected_dose_gy(int sdd_index){
    return(expected_dose_gy[sdd_index]);
}
const std::vector<long>* NGSsdd::get_chrom_size_bp(){
    return(&chrom_size_bp);
}
int NGSsdd::get_chrom_mapping(){
    return(chrom_mapping_type);
}
int NGSsdd::get_num_chrom(){
    return(num_chroms);
}
const std::vector<long>* NGSsdd::get_chrom_end_loc(){
    return(&chrom_end_loc);
}
long NGSsdd::get_sdd_genome_length(){
    return(chrom_end_loc.back());                                                                       // Return the last element of the vector
}
std::vector<std::streampos>* NGSsdd::get_lines_read_sdd(){
    return(&lines_read_sdd);
}
int NGSsdd::get_original_num_damages(int groupTID){
    return(original_num_damages[groupTID]);
}
std::streampos NGSsdd::get_exposureLines(int fileNumber, int exposureNumber){
    int vec_size = static_cast<int>(exposureLine_vec[fileNumber].size());
    std::streampos value;
    if(exposureNumber<vec_size){                                                                        // Since the index starts from 0 - vec_size-1, if the index given is within the limit, all ok
        value = exposureLine_vec[fileNumber][exposureNumber];
    }else{                                                                                              // If the index is out of bounds, then get the file end position
        value = exposureLine_vec[fileNumber][vec_size-1];     
    }
    return(value);                                                                                      // Return the appropriate element in the fileNumber vector
}
//std::vector<long>& NGSsdd::get_DNA_breakPoints(){
//    return(DNA_breakPoints);
//}
std::vector<long>& NGSsdd::get_backbone1_break_loc(int groupTID){
    return(backbone1_break_loc_vec[groupTID]);
}
std::vector<long>& NGSsdd::get_backbone2_break_loc(int groupTID){
    return(backbone2_break_loc_vec[groupTID]);
}
std::vector<long>& NGSsdd::get_basestrand1_damage_loc(int groupTID){
    return(basestrand1_damage_loc_vec[groupTID]);
}
std::vector<long>& NGSsdd::get_basestrand2_damage_loc(int groupTID){
    return(basestrand2_damage_loc_vec[groupTID]);
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// Reset functions to reset already set variables to re-use them
//--------------------------------------------------------------------------------------------
void NGSsdd::reset_workThread_data_holders(int groupTID){
    original_num_damages[groupTID] = 0;                                                                 // resetting the damage number counter to zero
    damage_start_loc_vec[groupTID].clear();                                                             // empty the damage start location vector of the group
    chrom_ID_vec[groupTID].clear();                                                                     // empty the chromosome ID vector of the group
}
void NGSsdd::reset_temporary_damage_vecs(int groupTID){                                                 // empty all the temporary vectors before the next SDD
    temp_backbone1_breaks_vec[groupTID].clear();
    temp_backbone2_breaks_vec[groupTID].clear();
    temp_basestrand1_damages_vec[groupTID].clear();
    temp_basestrand2_damages_vec[groupTID].clear();
}
void NGSsdd::reset_permanent_damage_vecs(int groupTID){                                                 // empty all the permanent vectors before the next exposure
    backbone1_break_loc_vec[groupTID].clear();
    backbone2_break_loc_vec[groupTID].clear();
    basestrand1_damage_loc_vec[groupTID].clear();
    basestrand2_damage_loc_vec[groupTID].clear();
    //DNA_breakPoints.clear();
}
//--------------------------------------------------------------------------------------------




//--------------------------------------------------------------------------------------------
// This function will read the group thread vector that has multiple worker thread vectors as
// elements (i.e, [ [],[],[] ] format) and contatenate the worker thread vectors into one to 
// get a single vector for the group (i.e, [ [] ] format)
//--------------------------------------------------------------------------------------------
void NGSsdd::merge_workerThread_vectors(int groupTID){
    
    // Process backbone 1 breaks
    std::vector<long> mergedVector_bb1;                                                                 // Temporary vector to hold the concatenated vector
    for (const auto& workVector_bb1 : temp_backbone1_breaks_vec[groupTID]){                             // Iterate over each work thread vector in the group thread
        mergedVector_bb1.insert(mergedVector_bb1.end(), workVector_bb1.begin(), workVector_bb1.end());  // Concatenate the elements of each workVector into the mergedVector
    }
    temp_backbone1_breaks_vec[groupTID].clear();                                                        // Empty the group vector and remove individual thread vectors
    temp_backbone1_breaks_vec[groupTID].push_back(mergedVector_bb1);                                    //  Push the concatenated vector into the group vector

    // Process backbone 2 breaks
    std::vector<long> mergedVector_bb2;                                                                 // Temporary vector to hold the concatenated vector
    for (const auto& workVector_bb2 : temp_backbone2_breaks_vec[groupTID]){                             // Iterate over each work thread vector in the group thread
        mergedVector_bb2.insert(mergedVector_bb2.end(), workVector_bb2.begin(), workVector_bb2.end());  // Concatenate the elements of each workVector into the mergedVector
    }
    temp_backbone2_breaks_vec[groupTID].clear();                                                        // Empty the group vector and remove individual thread vectors
    temp_backbone2_breaks_vec[groupTID].push_back(mergedVector_bb2);                                    //  Push the concatenated vector into the group vector

    // Process basestrand 1 breaks
    std::vector<long> mergedVector_bs1;                                                                 // Temporary vector to hold the concatenated vector
    for (const auto& workVector_bs1 : temp_basestrand1_damages_vec[groupTID]){                          // Iterate over each work thread vector in the group thread
        mergedVector_bs1.insert(mergedVector_bs1.end(), workVector_bs1.begin(), workVector_bs1.end());  // Concatenate the elements of each workVector into the mergedVector
    }
    temp_basestrand1_damages_vec[groupTID].clear();                                                     // Empty the group vector and remove individual thread vectors
    temp_basestrand1_damages_vec[groupTID].push_back(mergedVector_bs1);                                 //  Push the concatenated vector into the group vector

    // Process basestrand 2 breaks
    std::vector<long> mergedVector_bs2;                                                                 // Temporary vector to hold the concatenated vector
    for (const auto& workVector_bs2 : temp_basestrand2_damages_vec[groupTID]){                          // Iterate over each work thread vector in the group thread
        mergedVector_bs2.insert(mergedVector_bs2.end(), workVector_bs2.begin(), workVector_bs2.end());  // Concatenate the elements of each workVector into the mergedVector
    }
    temp_basestrand2_damages_vec[groupTID].clear();                                                     // Empty the group vector and remove individual thread vectors
    temp_basestrand2_damages_vec[groupTID].push_back(mergedVector_bs2);                                 //  Push the concatenated vector into the group vector
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This function will determine how many damages should be removed from the total number of 
// damages obtained in one exposure data of 1 SDD, using the actual dose delivered, if needed. 
// Then this function will also remove that many damages randomly from the temporary vector 
// holding the damage locations before pushing back the data to the final damage location vectors. 
// Final damages vectors will also be sorted.
//--------------------------------------------------------------------------------------------
void NGSsdd::adjust_damages_data(int exposure_index, int sdd_index, bool adjust_dose_flag, int groupTID){
    int num_damages_to_remove{0};                                                                       // temporary counter variable
    
    // If adjust_damage_with_actula_dose flag is ON
    if(adjust_dose_flag){
        const std::vector<std::vector<double>>& actual_dose = *get_actual_dose_delivered();
        num_damages_to_remove = std::round(get_original_num_damages(groupTID)*(1 - (get_expected_dose_gy(sdd_index)/(actual_dose[sdd_index][exposure_index]))));
        // damages to remove = total damages(1-(expected dose/actual dose));
    }

    while(num_damages_to_remove>0){
        //remove damage randomly from one temp vector at random
        switch (rng::rand_int(1,4,groupTID)){
            case 1:
                if(temp_backbone1_breaks_vec[groupTID][0].empty()){continue;}
                removeAnElement(temp_backbone1_breaks_vec[groupTID][0]);
                break;
            case 2: 
                if(temp_backbone2_breaks_vec[groupTID][0].empty()){continue;}
                removeAnElement(temp_backbone2_breaks_vec[groupTID][0]);
                break;
            case 3:
                if(temp_basestrand1_damages_vec[groupTID][0].empty()){continue;}
                removeAnElement(temp_basestrand1_damages_vec[groupTID][0]);
                break;
            case 4:
                if(temp_basestrand2_damages_vec[groupTID][0].empty()){continue;}
                removeAnElement(temp_basestrand2_damages_vec[groupTID][0]);
                break;
        }
        num_damages_to_remove--;
    }
    
    // finally push-back one element at a time to the final damage location vectors
    // destination.insert(destination.end(), source.begin(), source.end());
    backbone1_break_loc_vec[groupTID].insert(backbone1_break_loc_vec[groupTID].end(), temp_backbone1_breaks_vec[groupTID][0].begin(), temp_backbone1_breaks_vec[groupTID][0].end());
    backbone2_break_loc_vec[groupTID].insert(backbone2_break_loc_vec[groupTID].end(), temp_backbone2_breaks_vec[groupTID][0].begin(), temp_backbone2_breaks_vec[groupTID][0].end());
    basestrand1_damage_loc_vec[groupTID].insert(basestrand1_damage_loc_vec[groupTID].end(), temp_basestrand1_damages_vec[groupTID][0].begin(), temp_basestrand1_damages_vec[groupTID][0].end());
    basestrand2_damage_loc_vec[groupTID].insert(basestrand2_damage_loc_vec[groupTID].end(), temp_basestrand2_damages_vec[groupTID][0].begin(), temp_basestrand2_damages_vec[groupTID][0].end());
    
    // Sort the damage locations stored, in ascending order and remove duplicate damages coming from multiple SDDs.
    sortNremoveDuplicates_inVector(backbone1_break_loc_vec[groupTID]);
    sortNremoveDuplicates_inVector(backbone2_break_loc_vec[groupTID]);
    sortNremoveDuplicates_inVector(basestrand1_damage_loc_vec[groupTID]);
    sortNremoveDuplicates_inVector(basestrand2_damage_loc_vec[groupTID]);
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// Function to find the double-strand break (DSBs) locations, which are calculated as the mid
// point of two SSBs on opposite backbones and if they are within the DSB threshold provided, 
// and if they are on a single chromosome. Then remove backbone breaks contributed to the DSB
// from the list of backbone breaks. 
//--------------------------------------------------------------------------------------------
/* void NGSsdd::find_DNA_breakPoints(int DSBthreshold){
    std::vector<long>::iterator site1 = backbone1_break_loc.begin();
	std::vector<long>::iterator site2 = backbone2_break_loc.begin();
    // Proceed until have completely processed SSBs in either backbone
    while (site1 != backbone1_break_loc.end() && site2 != backbone2_break_loc.end()){
        int siteDiff = *site1 - *site2;                                                                 // separation in number of bp
		bool isDSB{0};                                                                                  // initiating with zero
        
        if(abs(siteDiff) <= DSBthreshold){                                                              // if the seperation is within threshold, then
            isDSB = is_nums_in_same_interval(*get_chrom_end_loc(), *site1, *site2);                     // check if these damages are on the same chromosome
        }

        if(isDSB){
            DNA_breakPoints.push_back(*site2+(std::round(siteDiff/2)));                                 // siteDiff will be -ve if site2>site1
            site1 = backbone1_break_loc.erase(site1);                                                   // remove the site if it's part of a DSB
            site2 = backbone2_break_loc.erase(site2);                                                   // next site will be assigned to site2 after erasing the old site2
        }else{
            if (siteDiff>0){site2++;}                                                                   // if site1 is after site2; update site2
            else if (siteDiff<0){site1++;}                                                              // if site1 is before site2; update site1
        }
    }
    sortNremoveDuplicates_inVector(DNA_breakPoints);                                                    // sort and remove duplicates from DNA breakpoints if there is any
} */
//--------------------------------------------------------------------------------------------