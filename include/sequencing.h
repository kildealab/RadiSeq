#ifndef SEQUENCING_H
#define SEQUENCING_H

#include "parameter_handler.h"

#include <iostream>
#include <string>
#include <vector>

void single_cell_sequencing(NGSParameters&, const std::vector<std::string>&, long);
void bulk_cell_sequencing(NGSParameters&, const std::vector<std::string>&, long);

#endif