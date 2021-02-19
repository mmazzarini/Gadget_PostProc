#ifndef LIBRARY_MMAZZARINI_MASS_DISTRIBUTION_H
#define LIBRARY_MMAZZARINI_MASS_DISTRIBUTION_H

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
//#include "externalModule_spline.h"

void bound_mass_calculate(std::vector<struct mass_center> *dcms, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp);

void integrate_MW_mass(struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2, std::vector<struct h5_particle> *pp3, struct header_h5 *hp);

void tidal_radius_calculate(struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, struct header_h5 *hp);

void calculate_TidalMass_Fraction(std::vector<struct mass_center> *dcms_Sat, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp,std::string filenumb);

void density_center_relocate(struct mass_center *dcms, std::vector<struct h5_particle> *pp, int *nPart);

void angular_momentum_MW_calculate(struct mass_center *dcms, std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2,  std::vector<struct h5_particle> *pp3,struct header_h5 *hp);

void angular_momentum_corrections(struct mass_center *dcms, std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2,  std::vector<struct h5_particle> *pp3, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5,struct header_h5 *hp);

void coordinates_rotate(float *lonNodeCont, float *diskIncCont,std::vector<struct h5_particle> *pp, int *nPart);

void mass_distribution_flag(std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp, int *thisFile);

void mass_distribution_disk(std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp);

void mass_distribution_bulge(std::vector<struct h5_particle> *pp2, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp);

void mass_distribution_halo(std::vector<struct h5_particle> *pp3, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, struct header_h5 *hp);

void calc_tensor_inertia(std::vector<struct h5_particle> *pp3, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp,struct mass_center *dcms,std::vector<struct mass_center> *dcms_Sat, std::string filenumb);

void calc_En_Ang_25kpc(std::vector<struct h5_particle> *pp3, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp, std::string myStringHalo, std::string myStringDark, std::string myStringStar);

void calc_Toomre_Space(std::vector<struct h5_particle> *pp3, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp, std::string myStringHalo, std::string myStringDark, std::string myStringStar);

void calc_Energy_Lang_DcmsSat(struct mass_center *dcms, std::vector<struct mass_center> *dcms_Sat, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp);

#endif
