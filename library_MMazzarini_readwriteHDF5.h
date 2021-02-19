#ifndef LIBRARY_MMAZZARINI_RDWRHDF5_H
#define LIBRARY_MMAZZARINI_RDWRHDF5_H

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

void write_Hdf5_postproc(std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2, std::vector<struct h5_particle> *pp3, std::vector<struct h5_particle> *pp4,std::vector<struct h5_particle> *pp5,struct header_h5 *hp, struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, struct mass_center *dcms_Baryon, struct mass_center *dcms_Dark);

void create_HDF5_arrays_file(struct header_h5 *hp,struct array_hd *ap);

void read_HDF5_file(std::vector<struct h5_particle> *pt0, std::vector<struct h5_particle> *pt1, std::vector<struct h5_particle> *pt2, std::vector<struct h5_particle> *pt3, std::vector<struct h5_particle> *pt4, std::vector<struct h5_particle> *pt5,struct header_h5 *hp);

void myFunc_Assign_H5parts(std::vector<struct h5_particle> *pt0, std::vector<struct h5_particle> *pt1, std::vector<struct h5_particle> *pt2, std::vector<struct h5_particle> *pt3, std::vector<struct h5_particle> *pt4, std::vector<struct h5_particle> *pt5, int *grindex, int *start, int *end, float myCoord[][3], float myVels[][3], float myAcc[][3],  float myMass[],  float myPot[], int myIDs[]);

void read_sbox_bin_MW(std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2, std::vector<struct h5_particle> *pp3, struct header_h5 *hp, int *myTime);

void read_sbox_bin_sat(std::vector<struct h5_particle> *ppsat, struct header_h5 *hp, int *myTime, int *sat_numb);

#endif
