#ifndef LIBRARY_MMAZZARINI_DENSCALC_H
#define LIBRARY_MMAZZARINI_DENSCALC_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include "library_MMazzarini_structs.h"


void read_densities(struct mass_center *dcms,std::vector<struct mass_center> *dcms_Sat, struct header_h5 *hp, int *myTime);

void density_center_calc_main(struct mass_center *dcms, std::vector<struct gas_particle> *gp, std::vector<struct dark_particle> *dp, std::vector<struct star_particle> *sp,struct header *hp);

void dcms_initialize(struct mass_center *cm);

void dcms_update(struct mass_center *ocm, struct mass_center *ncm);

void dcms_calculate(struct mass_center *ncm,struct mass_center *ocm,std::vector<struct gas_particle> *gp, std::vector<struct dark_particle> *dp, std::vector<struct star_particle> *sp,struct header *hp,float *maxRad,int *nparts,const int flag);

void dcms_determine(struct mass_center *dcms, struct mass_center *ncm);

void density_center_write(std::string myTipsyFile,struct mass_center *dcms);

void density_center_trunc(struct mass_center *dcms, std::vector<struct gas_particle> *gp, std::vector<struct dark_particle> *dp, std::vector<struct star_particle> *sp, struct header *hp,float *virRad);

void center_mass_increment_gas(struct mass_center *ncm, struct gas_particle *gpp);

void center_mass_increment_dark(struct mass_center *ncm, struct dark_particle *dpp);

void center_mass_increment_star(struct mass_center *ncm, struct star_particle *spp);

void density_center_calc_MW(struct mass_center *dcms, std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2,  std::vector<struct h5_particle> *pp3,struct header_h5 *hp);

void density_center_calc_Sat(std::vector<struct mass_center> *dcms_sat, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5,struct header_h5 *hp);

void dcms_calculate_MW(struct mass_center *ncm,struct mass_center *ocm, std::vector<struct h5_particle> *pp1, std::vector<struct h5_particle> *pp2,  std::vector<struct h5_particle> *pp3,struct header_h5 *hp,float *maxRad,int *nparts,const int flag);

void dcms_calculate_Sat(struct mass_center *ncm,struct mass_center *ocm, std::vector<struct h5_particle> *pp4, std::vector<struct h5_particle> *pp5, struct header_h5 *hp,float *maxRad,int *nparts,const int flag,int *ind_sat, int *initPartDark, int *finPartDark, int *initPartStar, int *finPartStar);

void center_mass_increment_h5Particle(struct mass_center *cm, struct h5_particle *pp);

void density_center_calc_SatelliteDark(std::vector<struct mass_center> *dcms_sat, std::vector<struct h5_particle> *pp5,struct header_h5 *hp);

void density_center_calc_SatelliteBaryon(struct mass_center *dcms_Baryon, std::vector<struct h5_particle> *pp4, struct header_h5 *hp);

void dcms_calculate_SatelliteBaryons(struct mass_center *ncm,struct mass_center *ocm, std::vector<struct h5_particle> *pp4, struct header_h5 *hp,float *maxRad,int *nparts,const int flag);

void dcms_calculate_SatelliteDark(struct mass_center *ncm,struct mass_center *ocm, std::vector<struct h5_particle> *pp5, struct header_h5 *hp,float *maxRad,int *nparts,const int flag);

#endif
