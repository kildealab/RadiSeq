#ifndef FASTAFILE_HANDLER_H
#define FASTAFILE_HANDLER_H

#include <string>

#include "sddfile_handler.h"

long buildUndamagedGenomeTemplate(const std::string&, int, int, const std::string*);
void getReverseComplementarySeq(std::string&, std::string&);
void buildDamagedCellGenome(NGSsdd&, const std::string&, const std::string&);

#endif