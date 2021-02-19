#ifndef LIBRARY_MMAZZARINI_ARRAYS_H
#define LIBRARY_MMAZZARINI_ARRAYS_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <hdf5.h>
#include "H5Cpp.h"
#include "library_MMazzarini_structs.h"

using namespace H5;

void kernel_general(struct header_h5 *hp, struct array_hd *ap);

void array_kernel(struct array_hd *ap);

void make_xyz_arrays(std::string myHdf5File, struct mass_center *dcms, std::vector<struct h5_particle> *pt, struct header_h5 *hp,struct array_hd *ap,int *pType);

#endif




