#ifndef LIBRARY_MMAZZARINI_UTILITIES_H
#define LIBRARY_MMAZZARINI_UTILITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "library_MMazzarini_structs.h"

std::string myString(int iFile);

void reset_structures(struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, std::vector<struct h5_particle> *pp0, std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2,  std::vector<struct h5_particle> *pp3,std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp, struct array_hd *ap);

bool acompare(struct h5_particle ptL, struct h5_particle ptR);

void read_particle_numbers(struct header_h5 *hp);


#endif
