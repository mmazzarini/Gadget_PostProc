#ifndef LIBRARY_MMAZZARINI_STRUCTS_H
#define LIBRARY_MMAZZARINI_STRUCTS_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>


struct mass_center{
  
  float mass;
  float pos[3];
  float vel[3];
  float angM[3];

};

//TIPSY Header

/*

struct header{

  double time;
  int    nbodies;
  int    nsph;
  int    ndark;
  int    nstar;
  int    ndim;
  
};

*/

//HDF5 header

struct header_h5{

  std::string myHDF5Root;
  std::string AquariusClass;
  std::string myHDF5File;
  std::string myHDF5Maps;
  std::string myHDF5Post;
  int first;
  int last;
  int step;
  float time;
  int sfr;
  int cool;
  int feed;
  float redshift;
  float boxSize;
  float omega0;
  float omegaL;
  float Hubble;
  int dPrec;	
  int nFile;
  float mass_table[6];
  int nPart_here_table[6];
  int nPart_tot_table[6];
  int nPart_high_table[6];
  int nbodies;
  int wantBaryons;
  int sat_numb;	
  int nsph; // Tags gas
  int nbody0; // 0 -> 5 tags all types of particles, independent from gas/no_gas
  int nbody1;
  int nbody2;
  int nbody3;
  int nbody4;
  int nbody5;
  std::vector<int> nbodySat = std::vector<int>(0);
  std::vector<float> softenings = std::vector<float>(0);

  // Tidal/bound quantities	

  float MW_concentration;
  float tot_halo_mass;	
  float mass_out_of_242_p4;
  float mass_out_of_242_p5;

  //std::vector<float> tot_halo_mass = std::vector<float>(0);
  std::vector<float> sat_massBound = std::vector<float>(0); //Describes the fracion of the initial satellite mass that is now still bound
  std::vector<float> sat_massBound_p4 = std::vector<float>(0);
  std::vector<float> sat_massBound_p5 = std::vector<float>(0);
  std::vector<float> sat_totMass = std::vector<float>(0); //Describes the initial satellite mass	
  std::vector<float> sat_totMass_p4 = std::vector<float>(0);
  std::vector<float> sat_totMass_p5 = std::vector<float>(0);
  std::vector<float> sat_tidRadius = std::vector<float>(0); //Tidal radius of satellite
  std::vector<float> MW_integratedMass = std::vector<float>(0); //MW mass within dens. center of MW and of satellite	   
  std::vector<float> sat_tidalMass = std::vector<float>(0);	
  std::vector<float> sat_tidalMass_p4 = std::vector<float>(0);
  std::vector<float> sat_tidalMass_p5 = std::vector<float>(0);
  std::vector<float> sat_ecc = std::vector<float>(0);
  std::vector<float> sat_ax = std::vector<float>(0);
  std::vector<float> sat_outTidMass = std::vector<float>(0);
  std::vector<float> sat_outTidMass_p4 = std::vector<float>(0);
  std::vector<float> sat_outTidMass_p5 = std::vector<float>(0);
  std::vector<int> myNumbDark = std::vector<int>(0);
  std::vector<int> myNumbStar = std::vector<int>(0);
 
  // We dont need the following 9 


  /*	
  float sat_massDisk;
  float sat_massDisk_p4;
  float sat_massDisk_p5;
  float sat_massBulge;
  float sat_massBulge_p4;
  float sat_massBulge_p5;
  float sat_massHalo;
  float sat_massHalo_p4;
  float sat_massHalo_p5;
  */	

  //Arrays
	//Disk

  float inertiaTensorDark[3][3][50];
  float inertiaTensorDark_Unbound[3][3][50];
  float inertiaTensorStar[3][3][50];
  float inertiaTensorStar_Unbound[3][3][50];
  float inertiaTensorHalo[3][3][50];
  float inertiaTensorHalo_Unbound[3][3][50];
  float geomTensorDark[3][3][50];
  float geomTensorDark_Unbound[3][3][50];
  float geomTensorStar[3][3][50];
  float geomTensorStar_Unbound[3][3][50];
  float geomTensorHalo[3][3][50];
  float geomTensorHalo_Unbound[3][3][50];
  float disk_Vtan_p1[40];
  float disk_Zed_p1[40];
  float disk_Vzed_p1[40];
  float disk_Mass_p1[40];
  float disk_SigVzed_p1[40];
  float disk_SigZed_p1[40];
  //float disk_vertical_total_p1[8];
  //float disk_vertical_zone_p1[4][20];
	//Bulge
  //float bulge_vels_p2[10];
  float bulge_Vrad_p2[10];
  float bulge_Vtheta_p2[10];
  float bulge_Vphi_p2[10];
  float bulge_Mass_p2[10];
  float bulge_Vphi_p4[10];
  float bulge_Vtheta_p4[10];
  float bulge_Vrad_p4[10];
  float bulge_Vphi_p5[10];
  float bulge_Vrad_p5[10];
  float bulge_Vtheta_p5[10];
  float bulge_Sigrad_p2[10];
  float bulge_Sigtheta_p2[10];
  float bulge_Sigphi_p2[10];
  float bulge_Sigma_p2[10];
  float disk_Mass_p4[40];	
  float bulge_Mass_p4[10];
  float halo_Mass_p4_tot;
  float disk_Mass_p5[40];
  float bulge_Mass_p5[10];
  float halo_Mass_p5_tot;

  //last step analysis	
  float bulge_vels_p4[10];
  float disk_Vtan_p4[40];
  float disk_Vzed_p4[40];
  float disk_Zed_p4[40];
  float bulge_vels_p5[10];
  float disk_Vtan_p5[40];
  float disk_Vzed_p5[40];
  float disk_Zed_p5[40];

  //Eccentricity and halo
  float halo_ecc_p3[40];
  float halo_ecc_p4[40];
  float halo_ecc_p5[40];	
  float smj_ax_p3[40];
  float smj_ax_p4[40];
  float smj_ax_p5[40];
  float halo_Mass_p3[40];
  float halo_Mass_hiRes_p3[50];
  float halo_Mass_hiRes2_p3[250];
  float halo_Mass_hiRes0_p3[25];
  float halo_Vrad_p3[40];
  float	halo_Vtheta_p3[40];
  float	halo_Vphi_p3[40];
  float	halo_Vrad_p4[40];
  float	halo_Vtheta_p4[40];
  float	halo_Vphi_p4[40];
  float	halo_Vrad_p5[40];
  float	halo_Vtheta_p5[40];
  float	halo_Vphi_p5[40];
  float	halo_SigVrad_p3[40];
  float	halo_SigVtheta_p3[40];
  float	halo_SigVphi_p3[40];
  float	halo_SigVrad_p4[40];
  float	halo_SigVtheta_p4[40];
  float	halo_SigVphi_p4[40];
  float	halo_SigVrad_p5[40];
  float	halo_SigVtheta_p5[40];
  float	halo_SigVphi_p5[40];	
  std::vector<float> halo_Mass_p4 = std::vector<float>(0);
  std::vector<float> halo_Mass_p5 = std::vector<float>(0);
  std::vector<float> halo_Mass_hiRes0_p4 = std::vector<float>(0);
  std::vector<float> halo_Mass_hiRes0_p5 = std::vector<float>(0);
  std::vector<float> halo_Mass_hiRes_p4 = std::vector<float>(0);
  std::vector<float> halo_Mass_hiRes_p5 = std::vector<float>(0);
  std::vector<float> halo_Mass_hiRes2_p4 = std::vector<float>(0);
  std::vector<float> halo_Mass_hiRes2_p5 = std::vector<float>(0);

};

	//Halo	
  /*
  float halo_Vrad_p3[20];		
  float halo_Vtheta_p3[20];
  float halo_Vphi_p3[20];	
  float halo_Mass_p3[20];
  float halo_Sigrad_p3[20];
  float halo_Sigtheta_p3[20];
  float halo_Sigphi_p3[20];
  float halo_vels_p3[20];
  */

	//Part4
  //float disk_Vtan_p4[40];


  /*
  float disk_vertical_total_p4[8];
  float disk_vertical_zone_p4[4][20];
  float bulge_vels_p4[10];
  float bulge_Vrad_p4[10];
  float bulge_Vtheta_p4[10];
  float bulge_Vphi_p4[10];
  float bulge_Sigrad_p4[10];
  float bulge_Sigtheta_p4[10];
  float bulge_Sigphi_p4[10];
  float halo_Vrad_p4[20];		
  float halo_Vtheta_p4[20];
  float halo_Vphi_p4[20];	
  float halo_Mass_p4[20];
  */
	

  /*
  float halo_Sigrad_p4[20];
  float halo_Sigtheta_p4[20];
  float halo_Sigphi_p4[20];
  float halo_vels_p4[20];
  */
	//Part5
 
  /*
  float disk_vertical_total_p5[8];
  float disk_vertical_zone_p5[4][20];
  float bulge_vels_p5[10];	
  float bulge_Vrad_p5[10];
  float bulge_Vtheta_p5[10];
  float bulge_Vphi_p5[10];
  */

  
  /*
  float bulge_Sigrad_p5[10];
  float bulge_Sigtheta_p5[10];
  float bulge_Sigphi_p5[10];
  float halo_Vrad_p5[20];		
  float halo_Vtheta_p5[20];
  float halo_Vphi_p5[20];	
  float halo_Mass_p5[20];
  float halo_Mass_p5_tot;	
  float halo_Sigrad_p5[20];
  float halo_Sigtheta_p5[20];
  float halo_Sigphi_p5[20];
  float halo_vels_p5[20];
  */




//TIPSY particles

struct gas_particle{

  float mass;
  float pos[3];
  float vel[3];
  float acc[3];	
  float rho;
  float temp;
  float hsmooth;
  float metals;
  float phi;
  int   idf;
};

struct dark_particle{

  float mass;
  float pos[3];
  float vel[3];
  float acc[3];	
  float eps;
  float phi;
  int   idf;

};

struct star_particle{

  float mass;
  float pos[3];
  float vel[3];
  float acc[3];	
  float metals;
  float tform;
  float eps;
  float phi;
  int   idf;
};

//Define a generic structure for a particle to appear in the HDF5 routines. This will allow me to use particles as gas, dark or star depending on the simulation setting only

struct h5_particle{

  float mass;
  float pos[3];
  float vel[3];
  float acc[3];	
  float rho;  //In case the particle is sph
  float temp; // " "
  float metals;//In case the particle is star
  float tform;// " "
  float eps;// " Softening"
  float phi;
  int   idf;
  int flag_bound;
  float potKin; //Sum of potential and kinetic energy of a particle
  int flag_zone;		

};

struct array_hd{


  int sat_numb;	
  float xMin_rel;
  float yMin_rel;	
  float zMin_rel;
  float xMax_rel;
  float yMax_rel;
  float zMax_rel;		
  std::vector<float> res = std::vector<float>(0);
  int nX[6];
  int nY[6];
  int nZ[6]; 
	 
};

#endif																										      	
